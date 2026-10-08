"""Generates the benchmark datasets into data/<name>/reads{1,2}.fa. Deterministic (fixed seeds).

The same three datasets as the merge cleanup benchmark of sbwt-rs-cli (benchmarks/merge_cleanup),
with colors added. Every dataset is two read sets from one random genome. Reads are dealt to the two
inputs so that neighbouring reads end up in different inputs: a read start is a dummy chain in its
own index that the other index often makes redundant.

Each read gets one of N_COLORS colors per input, written as the fasta header and used with
`themisto2 build --seq-colors --seq-colors-by-name`. The colors of the two inputs have different
names (a1..a16 and b1..b16), so the merged index has 32 colors.

  reads150  4M reads of 150 bp from a 200 Mbp genome. Few dummies (12% of the merged nodes).
  reads40   6M reads of 40 bp from a 50 Mbp genome. Fragmented: 64% of the merged nodes are dummies.
  kmerset   Single k-mers (31 bp reads) at 80% of the positions of a 15 Mbp genome; the k-mer at
            position i goes to input i % 2. Worst case for the cleanup of the SBWT merge.
"""
import os, random, sys

N_COLORS = 16

def write(name, reads):
    os.makedirs(f'data/{name}', exist_ok=True)
    outs = [open(f'data/{name}/reads{i}.fa', 'w') for i in (1, 2)]
    counts = [0, 0]
    for inp, seq in reads:
        color = counts[inp] % N_COLORS + 1
        counts[inp] += 1
        outs[inp].write(f'>{"ab"[inp]}{color}\n{seq}\n')
    for o in outs: o.close()

def random_reads(name, genome_len, n_reads, read_len, seed):
    random.seed(seed)
    genome = ''.join(random.choices('ACGT', k=genome_len))
    def gen():
        for r in range(n_reads):
            p = random.randrange(genome_len - read_len)
            yield r % 2, genome[p:p + read_len]
    write(name, gen())

def kmer_set(name, genome_len, density, k, seed):
    random.seed(seed)
    genome = ''.join(random.choices('ACGT', k=genome_len))
    def gen():
        for i in range(genome_len - k + 1):
            if random.random() < density:
                yield i % 2, genome[i:i + k]
    write(name, gen())

DATASETS = {
    'reads150': lambda: random_reads('reads150', 200_000_000, 4_000_000, 150, 1),
    'reads40':  lambda: random_reads('reads40', 50_000_000, 6_000_000, 40, 1),
    'kmerset':  lambda: kmer_set('kmerset', 15_000_000, 0.8, 31, 1),
}

if __name__ == '__main__':
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    for name in sys.argv[1:] or DATASETS:
        DATASETS[name]()
