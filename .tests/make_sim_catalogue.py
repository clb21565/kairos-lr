"""Simulated catalogue for integration tests: hosts with their own genes sharing one gene cassette."""
import random
import sys

rng = random.Random(7)
STOPS = {"TAA", "TAG", "TGA"}
CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS]
COMP = str.maketrans("ACGT", "TGCA")


def gene(n_codons=250):
    return "ATG" + "".join(rng.choice(CODONS) for _ in range(n_codons)) + "TAA"


def spacer(n=150):
    return "".join(rng.choice("ACGT") for _ in range(n))


def block(n_genes):
    return "".join(spacer() + gene(rng.randint(150, 350)) for _ in range(n_genes)) + spacer()


cassette = block(4)                       # the shared segment, ~4 kb
hosts = {
    "asm1_host1": block(6) + cassette + block(6),
    "asm2_host2": block(5) + cassette + block(7),
    "asm3_host3": block(8) + cassette.translate(COMP)[::-1] + block(4),
    "asm1_host4": block(10),               # no shared genes
}
with open(sys.argv[1], "w") as out:
    for name, seq in hosts.items():
        out.write(f">{name} simulated\n{seq}\n")
print(f"cassette {len(cassette)} bp", file=sys.stderr)
