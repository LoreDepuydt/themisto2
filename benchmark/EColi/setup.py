"""Extracts a subset of the 3682 E. coli assemblies and writes the --file-colors lists of two
input indexes for the merge benchmark. Deterministic (fixed seed).

Usage: python3 setup.py <n_per_index>     (from this directory, next to coli3682_dataset.tar.gz)

Takes the first 2 * n_per_index genomes of the archive (in archive order): the first half for
index 1, the second half for index 2. Every genome gets a random color c1..c10, from the same
10 names in both indexes, so that `themisto2 merge --merge-shared-colors` merges equally named
colors. Writes ecoli<n>/fof1.txt, ecoli<n>/fof2.txt and ecoli<n>/merge_args.txt.
"""
import os, random, subprocess, sys

n = int(sys.argv[1])
os.chdir(os.path.dirname(os.path.abspath(__file__)))
names = subprocess.run(['tar', '-tzf', 'coli3682_dataset.tar.gz'], capture_output=True, text=True, check=True).stdout.split()
genomes = [x for x in names if x.endswith('.fna')][:2 * n]
missing = [g for g in genomes if not os.path.exists(g)]
if missing:
    subprocess.run(['tar', '-xzf', 'coli3682_dataset.tar.gz'] + missing, check=True)

random.seed(1)
d = f'ecoli{n}'
os.makedirs(d, exist_ok=True)
for i, part in enumerate((genomes[:n], genomes[n:]), start=1):
    with open(f'{d}/fof{i}.txt', 'w') as f:
        for g in part:
            f.write(f'{os.path.abspath(g)}\tc{random.randint(1, 10)}\n')
with open(f'{d}/merge_args.txt', 'w') as f:
    f.write('--merge-shared-colors\n')
print(f'{d}: {n} + {n} genomes, {sum(os.path.getsize(g) for g in genomes) / 1e9:.2f} GB of sequence')
