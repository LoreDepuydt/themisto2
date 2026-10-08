"""Runs `themisto2 merge` on the benchmark datasets and reports time and memory per phase.

Usage: python3 run.py [--low-ram] [--repeat N] <tag> [dataset ...]     (after `cargo build --release`)

Builds the two input indexes of each dataset once (k = 31, 8 threads, colors from the fasta headers),
then merges them with 8 threads while sampling the resident memory (RSS) of the process every 20 ms.
Every log line is timestamped when it arrives (the log itself has only whole seconds), and the
samples are matched to the phases below. Logs go to data/<dataset>/log_<tag>.txt.

Phases (each runs until the next one starts):
  load          loading the two input indexes
  interleaving  computing the merge interleaving of the two SBWTs
  sbwt merge    merging the SBWTs, removing redundant dummies, building select support
    cleanup     of which the removal of redundant dummy nodes (a sub-phase, not a separate row in time)
  (dummy marking: the time spent in "Marking dummy nodes" of sbwt::dbg, summed over the DBGs of the
   merged index (in dbg) and of the two inputs (in colors 1))
  lcs           building the LCS array of the merged SBWT
  dbg           initializing the DBG of the merged SBWT (the DBGs of the inputs are initialized in colors 1)
  colors 1/2/3  the three phases of merging the colors (key k-mers, fingerprints, color set storage)
  write         writing the merged index

For each phase the average number of busy cores (CPU time / wall time) is reported too, and after
the table the stretches of at least 3 s in which at most 1.5 cores were busy, with the log message
that was current when they began.

A dataset is a directory data/<name>, or a path relative to this directory if the name contains a
'/' (e.g. ../EColi/ecoli100). Its inputs are built from reads1.fa and reads2.fa (with
--seq-colors), or from fof1.txt and fof2.txt (with --file-colors) if those exist. If it has a
file merge_args.txt, its contents are added to the merge command (e.g. --merge-shared-colors).

With MEMORY_MAX set (e.g. MEMORY_MAX=10G), every build and merge runs under that memory limit
(systemd-run --user --scope -p MemoryMax=...), so that a process that grows too large is killed
instead of the whole machine running out of memory.

With --low-ram, the merge uses its low-memory mode. With --repeat N, each merge runs N times and the
median of each number is reported. The binary is target/release/themisto2, or $THEMISTO2_BIN if set.
The script refuses to run target/release/themisto2 if it is older than the sources of themisto2 or of
the local sbwt checkout it is built against.
"""
import datetime, os, re, subprocess, sys, threading, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '../..')
SBWT_SRC = os.path.join(ROOT, '../externalSoftware/SBWT/sbwt-rs-cli/api/src')  # Path dependency in Cargo.toml
DEFAULT_BIN = os.path.join(ROOT, 'target/release/themisto2')
BIN = os.environ.get('THEMISTO2_BIN') or DEFAULT_BIN
K, THREADS = 31, 8
DATASETS = ['reads150', 'reads40', 'kmerset']

# (phase name, prefix of the log message that starts it)
PHASES = [
    ('load', 'Merging '),  # "Merging <index 1> and <index 2> into ...", the first such message
    ('interleaving', 'Computing the sbwt merge plan'),
    ('sbwt merge', 'Merging SBWTs'),
    ('lcs', 'Building the LCS array for the merged SBWT'),
    ('dbg', 'Initializing DBG for the merged SBWT'),
    ('colors 1', '=== Phase 1/3'),
    ('colors 2', '=== PHASE 2/3'),
    ('colors 3', '=== PHASE 3/3'),
    ('write', 'Serializing merged index'),
]

def check_binary_is_current():
    if BIN != DEFAULT_BIN:
        return  # A chosen binary may be older on purpose
    sources = [os.path.join(d, f) for src in (os.path.join(ROOT, 'src'), SBWT_SRC)
               for d, _, files in os.walk(src) for f in files if f.endswith('.rs')]
    newest = max(sources, key=os.path.getmtime)
    if os.path.getmtime(BIN) < os.path.getmtime(newest):
        sys.exit(f'{BIN} is older than {newest}: run cargo build --release in the themisto2 root')

def sbwt_version():
    """The commit of the local sbwt checkout, with a + if it has uncommitted changes."""
    git = ['git', '-C', SBWT_SRC]
    commit = subprocess.run(git + ['rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()
    branch = subprocess.run(git + ['branch', '--show-current'], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(git + ['status', '--porcelain', '--', '.'], capture_output=True, text=True).stdout.strip()
    return f'{branch} {commit}{"+" if dirty else ""}'

def limited(cmd):
    """The command, under the memory limit MEMORY_MAX if that is set."""
    limit = os.environ.get('MEMORY_MAX')
    if not limit:
        return cmd
    return ['systemd-run', '--user', '--scope', '-q', '-p', f'MemoryMax={limit}', '-p', 'MemorySwapMax=0'] + cmd

def dataset_dir(name):
    return name if '/' in name else f'data/{name}'

def build_inputs(d):
    os.makedirs(f'{d}/tmp', exist_ok=True)
    for i in (1, 2):
        if not os.path.exists(f'{d}/idx{i}.thm2'):
            if os.path.exists(f'{d}/fof{i}.txt'):
                colors = ['--file-colors', f'{d}/fof{i}.txt']
            else:
                colors = ['--seq-colors', f'{d}/reads{i}.fa', '--seq-colors-by-name']
            with open(f'{d}/build{i}.log', 'w') as log:
                subprocess.run(limited([BIN, 'build'] + colors + ['-o', f'{d}/idx{i}.thm2', '-k', str(K),
                                '-t', str(THREADS), '--temp-dir', f'{d}/tmp']),
                               check=True, stderr=log, stdout=log)
    with open(f'{d}/inputs.txt', 'w') as f:
        f.write(f'{d}/idx1.thm2\n{d}/idx2.thm2\n')

def run_merge(d, tag, extra_args):
    log_path = f'{d}/log_{tag}.txt'
    merge_args = open(f'{d}/merge_args.txt').read().split() if os.path.exists(f'{d}/merge_args.txt') else []
    p = subprocess.Popen(limited([BIN, 'merge', '--index-file-list', f'{d}/inputs.txt', '--temp-dir', f'{d}/tmp',
                          '-o', f'{d}/merged.thm2', '-t', str(THREADS)] + merge_args + extra_args),
                         stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    messages = []  # (arrival time, message)
    def read_log():
        with open(log_path, 'w') as log:
            for line in p.stderr:
                now = datetime.datetime.now()
                log.write(f'{now.isoformat(timespec="milliseconds")} {line}')
                m = re.match(r'\[\S+ \w+\s+[^\]]+\] (.*)', line)
                messages.append((now, m.group(1) if m else line.rstrip('\n')))
    reader = threading.Thread(target=read_log)
    reader.start()
    samples = []
    while p.poll() is None:
        try:
            status = open(f'/proc/{p.pid}/status').read()
            stat = open(f'/proc/{p.pid}/stat').read().rsplit(')', 1)[1].split()
            cpu_seconds = (int(stat[11]) + int(stat[12])) / os.sysconf('SC_CLK_TCK')  # utime + stime
            samples.append((datetime.datetime.now(), int(re.search(r'VmRSS:\s+(\d+)', status).group(1)), cpu_seconds))
        except (OSError, AttributeError, IndexError):
            pass
        time.sleep(0.02)
    reader.join()
    end = datetime.datetime.now()
    assert p.returncode == 0, f'merge failed, see {log_path}'
    return samples, messages, end

def first(messages, prefix, after=None):
    return next(t for t, msg in messages if msg.startswith(prefix) and (after is None or t >= after))

def peak(samples, start, end):
    return max((r for t, r, _ in samples if start <= t < end), default=0) / 1024

def cores_busy(samples, start, end):
    """Average number of busy cores between start and end (CPU time / wall time)."""
    inside = [(t, c) for t, _, c in samples if start <= t < end]
    if len(inside) < 2:
        return 0.0
    return (inside[-1][1] - inside[0][1]) / (inside[-1][0] - inside[0][0]).total_seconds()

def single_threaded_stretches(samples, messages, min_seconds=3.0, max_cores=1.5):
    """Stretches of at least min_seconds in which at most max_cores cores were busy on average,
    measured in 1-second windows, with the log message that was current when each stretch began."""
    windows = []  # (start, end, cores busy)
    i = 0
    while i < len(samples):
        j = i
        while j < len(samples) and (samples[j][0] - samples[i][0]).total_seconds() < 1.0:
            j += 1
        if j >= len(samples):
            break
        windows.append((samples[i][0], samples[j][0], (samples[j][2] - samples[i][2]) / (samples[j][0] - samples[i][0]).total_seconds()))
        i = j
    stretches, current = [], None
    for start, end, busy in windows:
        if busy <= max_cores:
            current = [current[0] if current else start, end]
        else:
            if current and (current[1] - current[0]).total_seconds() >= min_seconds:
                stretches.append(tuple(current))
            current = None
    if current and (current[1] - current[0]).total_seconds() >= min_seconds:
        stretches.append(tuple(current))
    result = []
    for start, end in stretches:
        message = max(((t, m) for t, m in messages if t <= start), default=(start, '(start)'))[1]
        result.append(((end - start).total_seconds(), message))
    return result

def measure(samples, messages, end):
    """[merged nodes (M), removed (M), total time, total peak, then time and peak per phase, cleanup time and peak]."""
    starts = []
    for _, prefix in PHASES:  # Each phase starts after the previous one
        starts.append(first(messages, prefix, starts[-1] if starts else None))
    starts.append(end)
    sets = next((int(re.search(r'(\d+)$', m).group(1)) for _, m in messages if m.startswith('Number of sets in merged SBWT')), 0)
    removed = next((int(re.search(r'removing (\d+)', m).group(1)) for _, m in messages if m.startswith('[strip] Filtering')), 0)
    values = [sets / 1e6, removed / 1e6, (end - starts[0]).total_seconds(), max(r for _, r, _ in samples) / 1024]
    for a, b in zip(starts, starts[1:]):
        values += [(b - a).total_seconds(), peak(samples, a, b), cores_busy(samples, a, b)]
    cleanup = [t for t, m in messages if m.startswith('[merge] Removing redundant dummy nodes')]
    if cleanup:
        cleanup_end = next(t for t, m in messages if t > cleanup[0] and m.startswith('Building the subset rank structure'))
        values += [(cleanup_end - cleanup[0]).total_seconds(), peak(samples, cleanup[0], cleanup_end)]
    else:
        values += [0.0, 0.0]
    dummy_marking = sum((b[0] - a[0]).total_seconds() for a, b in zip(messages, messages[1:])
                        if a[1].startswith('Marking dummy nodes'))
    return values + [dummy_marking]

def report(name, values):
    sets, removed, total_time, total_peak = values[:4]
    print(f'\n#### {name}: {sets:.0f}M merged nodes, {removed:.1f}M removed, {total_time:.1f} s, peak {total_peak:.0f} MB\n')
    print('| phase | time | peak | cores busy |')
    print('|---|---|---|---|')
    for i, (phase, _) in enumerate(PHASES):
        print(f'| {phase} | {values[4 + 3*i]:.2f} s | {values[5 + 3*i]:.0f} MB | {values[6 + 3*i]:.1f} |')
    print(f'| (cleanup, part of sbwt merge) | {values[-3]:.2f} s | {values[-2]:.0f} MB |')
    print(f'| (dummy marking, part of dbg and colors 1) | {values[-1]:.2f} s | |', flush=True)

if __name__ == '__main__':
    args = sys.argv[1:]
    extra_args = ['--low-ram-mode'] if '--low-ram' in args else []
    args = [a for a in args if a != '--low-ram']
    repeat = 1
    if '--repeat' in args:
        i = args.index('--repeat')
        repeat = int(args[i + 1])
        args = args[:i] + args[i + 2:]
    tag, names = args[0], args[1:] or DATASETS
    os.chdir(HERE)
    check_binary_is_current()
    print(f'themisto2 merge, sbwt {sbwt_version()}, {THREADS} threads{", low-RAM mode" if extra_args else ""}, '
          f'median of {repeat}')
    for name in names:
        d = dataset_dir(name)
        build_inputs(d)
        runs, stretches = [], None
        for i in range(repeat):
            samples, messages, end = run_merge(d, f'{tag}_{i}' if repeat > 1 else tag, extra_args)
            runs.append(measure(samples, messages, end))
            stretches = stretches or single_threaded_stretches(samples, messages)
        report(name, [sorted(column)[len(column) // 2] for column in zip(*runs)])
        if stretches:
            print('\nStretches of at least 3 s with at most 1.5 cores busy (first run):\n')
            for seconds, message in stretches:
                print(f'- {seconds:.0f} s, from: {message[:100]}')
