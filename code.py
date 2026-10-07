import os

# Table standard du code génétique (ADN -> Acide Aminé)
CODON_TABLE = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGA": "*",
    "TGC": "C",
    "TGT": "C",
    "TGG": "W",
}


def lire_fasta(chemin):
    with open(chemin, "r") as f:
        lignes = f.readlines()
    return "".join(
        line.strip().upper() for line in lignes if not line.startswith(">")
    )


def traduire_dna(dna_seq):
    prot = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i : i + 3]
        prot.append(CODON_TABLE.get(codon, "X"))
    return "".join(prot)


# Chargement des données depuis le dossier 'donnees'
# Modifiez le chemin pour pointer vers "ressources/donnees"
dna_seq = lire_fasta(os.path.join("ressources", "donnees", "sequence.fasta"))
prot_seq = lire_fasta(os.path.join("ressources", "donnees", "geneX.fasta"))

motif = prot_seq[:6]
print(f"Motif recherché : {motif}\n")

for frame in [1, 2, 3]:
    dna_sub = dna_seq[frame - 1 :]
    translated = traduire_dna(dna_sub)
    if motif in translated:
        pos = translated.find(motif)
        print(f"✅ TROUVÉ ! La protéine X est dans le CADRE +{frame}")
        print(f"   Position du premier acide aminé (M) : index {pos}")