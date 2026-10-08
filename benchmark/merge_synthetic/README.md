# Merge benchmark on synthetic data

Measures the time and memory of each phase of `themisto2 merge`, on generated inputs with few and
with many dummy nodes. Each change to the merge is measured here before and after, and listed below
with the dataset that shows why it is useful. The SBWT part of the merge has its own, more detailed
benchmark in sbwt-rs-cli (`benchmarks/merge_cleanup`, on the branch `benchmark-tracking`).

## Running

```sh
cargo build --release
cd benchmark/merge_synthetic
python3 gen.py                   # writes data/<dataset>/reads{1,2}.fa, deterministic (about 2 GB)
python3 run.py <tag>             # builds the inputs once, merges, prints a table per dataset
python3 run.py --low-ram <tag>   # the same, with --low-ram-mode
python3 run.py --repeat 3 <tag>  # each merge three times, reports the median of each number
```

`data/` is not committed. Building the six input indexes takes about 25 minutes. `run.py` runs
`target/release/themisto2` and refuses to if that binary is older than the sources of themisto2 or
of the local sbwt checkout (`../externalSoftware/SBWT/sbwt-rs-cli`, see `Cargo.toml`), whose branch
and commit it prints. `THEMISTO2_BIN=<path>` runs another binary instead. The phases are explained
at the top of `run.py`.

## Datasets

The same reads as the sbwt benchmark (same generator and seeds), with colors: every read gets one
of 16 colors of its input, so the merged index has 32 colors. themisto2 also indexes the reverse
complements, so the indexes have about twice as many nodes as in the sbwt benchmark.

| dataset | input | merged nodes | removed dummies |
|---|---|---|---|
| reads150 | 4M reads of 150 bp | 412M | 33.2M |
| reads40 | 6M reads of 40 bp | 193M | 54.8M |
| kmerset | 12M single k-mers | 469M | 350.5M |

## Results

Intel Core Ultra 7 265H, 30 GB RAM, 8 threads.

### Baseline (themisto2 5365720, sbwt benchmark-tracking 44f5d56)

One run each. Time and peak resident memory per phase:

| phase | reads150 | reads40 | kmerset |
|---|---|---|---|
| load | 0.9 s, 881 MB | 0.3 s, 323 MB | 0.6 s, 676 MB |
| interleaving | 72.2 s, 2483 MB | 30.1 s, 998 MB | 48.0 s, 2172 MB |
| sbwt merge | 8.2 s, 1674 MB | 3.7 s, 766 MB | 8.7 s, 1665 MB |
| lcs | 35.2 s, 2261 MB | 13.6 s, 964 MB | 12.6 s, 1603 MB |
| dbg | 13.0 s, 1945 MB | 46.1 s, 849 MB | 61.0 s, 1500 MB |
| colors 1 | 243.4 s, 2174 MB | 128.2 s, 918 MB | 431.2 s, 1630 MB |
| colors 2 | 30.0 s, **2646 MB** | 3.6 s, **1325 MB** | 3.7 s, 2107 MB |
| colors 3 | 3.4 s, 2172 MB | 0.5 s, 958 MB | 1.0 s, 1586 MB |
| write | 0.5 s, 849 MB | 0.3 s, 586 MB | 0.2 s, 634 MB |
| **total** | **406.5 s, 2646 MB** | **226.5 s, 1325 MB** | **567.1 s, 2172 MB** |
| of which dummy marking | 50.4 s (12%) | 125.8 s (56%) | 407.4 s (72%) |
| of which removing redundant dummies | 5.0 s | 2.5 s | 6.9 s |

- **Dummy marking** is the largest single cost on inputs with many dummy nodes: 56% of the merge on
  `reads40` and 72% on `kmerset`. Every `Dbg::new` computes the dummy marks with
  `SbwtIndex::compute_dummy_node_marks`, a depth-first search from the root on one thread, and the
  merge creates three DBGs: of the merged index (in *dbg*) and of both inputs (in *colors 1*).
- **colors 1** is the slowest phase on every dataset. Besides the dummy marking of the two input
  DBGs, it iterates the unitigs of three DBGs; on `reads150` each iteration takes about 67 s.
- **Memory**: the peak is in *colors 2* on `reads150` and `reads40`, and in the interleaving on
  `kmerset`. The SBWT merge, including the removal of redundant dummies, is not the peak on any
  dataset, and takes 2 to 4% of the time.

The second run of the baseline binary gave totals within 2 to 4% of the first (416.7, 234.4 and
582.8 s) and peaks within 3 MB, so single runs are used below.

### 1. Parallel dummy marking (sbwt parallel-dummy-marks, kept as an option on its own branch)

`SbwtIndex::compute_dummy_node_marks` found the dummy nodes with a depth-first search from the
root, on one thread. It now walks the dummy nodes level by level from the root (the node at
distance d has k - d leading dollars, so there are at most k levels), with each level a bit vector
that is processed in parallel; `Dbg::new` runs it with its `n_threads`. This needs 2n bits more
memory while it runs (the current and the next level), where the depth-first search needed only a
stack of at most about 3k entries.

| phase | reads150 | reads40 | kmerset |
|---|---|---|---|
| dbg | 13.0 → 2.9 s | 46.1 → 1.0 s | 61.0 → 1.0 s |
| colors 1 | 243.4 → 203.8 s | 128.2 → 50.1 s | 431.2 → 91.6 s |
| of which dummy marking (all three DBGs) | 50.4 → 4.8 s | 125.8 → 1.6 s | 407.4 → 4.1 s |
| **total** | **406.5 → 357.7 s** | **226.5 → 104.5 s** | **567.1 → 168.1 s** |
| peak | 2646 → 2649 MB | 1325 → 1375 MB | 2172 → 2191 MB |

The dummy marking is 10 to 100 times faster, more than the 8 threads explain. Probably the access
pattern matters as well: each level is visited in colex order, so the set lookups and rank queries
read memory roughly in order, where the depth-first search jumped along one path at a time (not
measured separately). The whole merge takes 12, 54 and 70% less time. The peak grows by 50 MB on `reads40`
and 19 MB on `kmerset` (it is in *colors 2*, after the marking, so most likely memory the allocator
kept), and on `kmerset` the peak moves from the interleaving to *colors 2*. The merged indexes of
`reads40` and `kmerset` are byte for byte identical to those of the baseline.

The higher peaks are memory that glibc keeps after it is freed, not memory in use. glibc returns a
freed block to the OS only if it is above a threshold that grows with the blocks it has seen, up to
32 MB. The bit vectors of the levels are about 24 MB on `reads40` (kept) and 51 to 58 MB on
`reads150` and `kmerset` (returned). With a fixed threshold (`MALLOC_MMAP_THRESHOLD_=1048576`),
both versions on `reads40`:

| phase | baseline | parallel | difference |
|---|---|---|---|
| dbg (marking of the merged index) | 652 MB | 685 MB | +33 MB |
| colors 1 (marking of the inputs) | 734 MB | 747 MB | +13 MB |
| colors 2 (peak) | 1113 MB | 1113 MB | 0 |
| write | 244 MB | 244 MB | 0 |
| total time | 237.0 s | 106.8 s | |

So the real cost is the 2n bits while the marks are computed (48 MB for the merged index of
`reads40`), in phases that are far below the peak of the merge. It would only raise the peak of a
program in which `Dbg::new` runs while little else is in memory; `themisto2 build` and `export` are
not measured yet. If that matters, the levels could be lists of positions while they are small.

Side finding: with the fixed threshold, the peak of the baseline itself drops from 1325 to 1113 MB
(16%) at the same running time, so the merge loses about 200 MB to memory that glibc keeps.

Remaining: on `reads150`, *colors 1* still takes 204 of 358 s, mostly iterating the unitigs of the
three DBGs (about 67 s each).

### 2. Fixed glibc mmap threshold

glibc returns a freed block to the operating system only if it is larger than a threshold that by
default grows with the largest block freed so far, up to 32 MiB. Many structures of the merge are
tens of MiB, so after a while they stay in the heap of the process after they are freed, and later
phases run on top of them. `main` now fixes the threshold at 1 MiB with `mallopt` (Linux with glibc
only). Measured with `MALLOC_MMAP_THRESHOLD_=1048576`, which has the same effect, on top of version 1:

| | reads150 | reads40 | kmerset |
|---|---|---|---|
| peak | 2649 → 2405 MB (−9%) | 1375 → 1113 MB (−19%) | 2191 → 1930 MB (−12%) |
| total time | 357.7 → 367.0 s | 104.5 → 104.5 s | 168.1 → 172.9 s |

Every phase after loading is lower (e.g. *write* on `reads150`: 851 → 593 MB), and the times are
within the noise of 2 to 4%. With the `mallopt` call in the code instead of the environment
variable, `reads40` gives the same result (peak 1113 MB, 103.6 s) and a byte for byte identical
merged index.

### 3. Faster unitig walk (sbwt faster-unitig-walk, on top of the baseline)

`Dbg::walk_unitig_from` looked up the (k-1)-suffix group of each node four times, and computed
the in-degree of the next node with a select (`inverse_lf_step`) that leads back to the group it
came from. It now looks up the group once per node and computes the in-degree of the next node
from that group: one rank per node instead of a rank and a select. Measured against the baseline
(old dummy marking, no fixed mmap threshold), one run (the unit tests check that the walks are
identical; the merged indexes were not compared here):

| | baseline (two runs) | faster walk |
|---|---|---|
| reads150, total | 406.5 / 416.7 s | 367.3 s (−10%) |
| reads150, the three unitig walks | 67 / 68 / 67 s | 49 / 57 / 57 s |
| reads40, total | 226.5 / 234.4 s | 232.0 s |
| kmerset, total | 567.1 / 582.8 s | 576.7 s |

On `reads40` and `kmerset` *colors 1* is dominated by the dummy marking and by the callbacks of
the walk in themisto2, and the unitigs are short, so the difference is within the noise there.
Peaks do not change.

A first version created the node vector with `vec![v]` (room for one node) and was 55% slower on
`reads40`: `perf` showed the time in `realloc` and in waiting for the allocator lock, from growing
the vector of every short unitig on 8 threads. With room for 4 nodes it is faster on both. The
callbacks of themisto2 in *colors 1* also allocate several vectors per unitig, which may be worth
a look for the same reason.

### Current state: 1 + 3 (sbwt parallel-dummy-marks c7b53ac)

Both sbwt changes together, without the fixed mmap threshold. `run.py` now also reports the
average number of busy cores per phase (8 threads), and the stretches of at least 3 s with at most
1.5 cores busy.

| phase | reads150 | reads40 | kmerset |
|---|---|---|---|
| interleaving | 72.6 s, 6.2 cores | 30.0 s, 6.2 cores | 48.0 s, 6.2 cores |
| sbwt merge | 8.1 s, 5.6 cores | 3.8 s, 5.9 cores | 9.1 s, 6.4 cores |
| lcs | 34.7 s, 6.3 cores | 13.7 s, 6.3 cores | 12.6 s, 6.3 cores |
| dbg | 2.9 s, 1.9 cores | 1.2 s, 3.7 cores | 1.0 s, 4.7 cores |
| colors 1 | 162.8 s, 7.6 cores | 49.3 s, 7.5 cores | 100.8 s, 7.3 cores |
| colors 2 | 29.8 s, 7.3 cores | 3.6 s, 5.0 cores | 3.6 s, 3.9 cores |
| **total** | **315.3 s** (baseline 407) | **102.7 s** (baseline 227) | **177.0 s** (baseline 567) |
| peak | 2647 MB | 1366 MB | 2157 MB |

The merge keeps most of the 8 cores busy. The only single-threaded stretches of 3 s or more are
3 to 4 s in the unitig iterations (the end of the acyclic part, where the last long unitigs are
walked, and the cyclic part, which runs on one thread), two per dataset at most.

Keep in mind that these datasets are synthetic and small: one random genome without repeats,
random colors per read (so the color set changes every few k-mers, unlike real genomes), and
inputs dealt so that many dummies become redundant. They are good for finding mechanisms and
worst cases, but priorities should be checked on real data.

### 4. Positions in the merged index from the interleaving (themisto2)

In *colors 1*, the callback of the unitig walk over each input index needs the position in the
merged index of the first and last k-mer of every colored subunitig. It built the string of every
unitig (`push_unitig_string`, which reconstructs the first k-mer with k selects) and searched the
two k-mers in the merged index. A `perf` profile of `reads40` (built with frame pointers) showed
this string as the largest single cost: 13.7% of the whole merge, against 0.6% for the searches.
The position now comes from the interleaving: input k-mer i is at interleaving position
select(s1, i) (s2 for the second input), minus the number of removed positions before it. This
needs a copy of s1 or s2 with select support and the removed positions with rank support while
*colors 1* runs.

On top of the current state (1 + 3), one run each; the merged index of `reads40` is byte for byte
identical:

| | reads150 | reads40 | kmerset |
|---|---|---|---|
| walk + callback per input index | 55.9 / 53.5 → 45.2 / 44.9 s | 20.1 / 20.1 → 8.2 / 8.2 s | 43.3 / 49.5 → 4.8 / 4.8 s |
| colors 1 | 162.8 → 145.7 s | 49.3 → 26.3 s | 100.8 → 18.1 s |
| **total** | **315.3 → 300.9 s** (−5%) | **102.7 → 80.4 s** (−22%) | **177.0 → 96.3 s** (−46%) |
| peak of colors 1 | 2175 → 2327 MB | 952 → 1034 MB | 1621 → 1815 MB |
| peak of the run | 2647 → 2664 MB | 1366 → 1382 MB | 2157 → 2183 MB |

The gain is largest where unitigs are short (the string costs k selects per unitig, however short
the unitig). The extra memory is below the peak of the run, which stays in *colors 2*. On
`kmerset`, *colors 1* now keeps only 3.9 cores busy: of its 18 s, two stretches of 3 to 4 s run on
one thread, in the cyclic part of the unitig iteration.

## Real data: E. coli assemblies

`../EColi/setup.py 100` extracts the first 200 of the 3682 E. coli assemblies
(`coli3682_dataset.tar.gz`, not committed) and writes two `--file-colors` lists of 100 genomes
each. Every genome gets a random color c1..c10 (seed 1), with the same 10 names in both indexes,
and the merge uses `--merge-shared-colors`. Run with
`MEMORY_MAX=10G python3 run.py <tag> ../EColi/ecoli100`.

The merged index has 109M nodes, of which 0.1% dummies (14k redundant ones are removed), 4.6M
unitigs (about 23 k-mers each, because of the variation between the genomes), 10 shared colors and
1023 distinct color sets. Like the merge of two human haplotypes (assemblies, few dummies), it
spends most of its time in *colors 1*, and the dummy marking does not matter (0.3 s).

Baseline (themisto2 5365720, sbwt benchmark-tracking 44f5d56) against the current state (sbwt
1 + 3, themisto2 change 4), one run each; the merged indexes are byte for byte identical:

| phase | baseline | current |
|---|---|---|
| interleaving | 19.7 s, 6.2 cores | 19.6 s, 6.2 cores |
| sbwt merge | 2.1 s | 2.0 s |
| lcs | 11.1 s | 11.2 s |
| colors 1 | 61.9 s | 29.2 s |
| of which the three unitig walks | 13.3 / 25.7 / 22.4 s | 11.3 / 8.5 / 8.4 s |
| colors 2 + 3 | 3.0 s | 3.0 s |
| **total** | **98.8 s** | **65.7 s** (−33%) |
| peak (in colors 2) | 771 MB | 841 MB |

With a fixed mmap threshold (`MALLOC_MMAP_THRESHOLD_=1048576`, see change 2) both peaks are 674 MB:
the 70 MB difference is memory that glibc keeps, not memory in use. The real extra memory of change
4 is 40 MB in *colors 1* (532 → 572 MB), below the peak. The fixed threshold itself lowers the
peak by 13% (baseline) and 20% (current).

Now the interleaving (30%) and the LCS array (17%) are the largest parts after *colors 1* (44%).
