#!/usr/bin/python3

from Bio import AlignIO
from collections import defaultdict

alignment = AlignIO.read("proteins_aln.fasta", "fasta")

# Define your capripoxvirus samples
capripox_ids = ["KX894508_proteins", "NC_004002_proteins", "NC_004003_proteins"]  # Replace with your actual sequence IDs

# Create lists of sequences for each group
capripox_seqs = [seq for seq in alignment if seq.id in capripox_ids]
non_capripox_seqs = [seq for seq in alignment if seq.id not in capripox_ids]

type1_sites = []  # Conserved in CPPV, different in others
type2_sites = []  # Present in CPPV, gap in all others
for i in range(alignment.get_alignment_length()):
    capripox_residues = [seq[i] for seq in capripox_seqs]
    non_capripox_residues = [seq[i] for seq in non_capripox_seqs]
    if len(set(capripox_residues)) == 1:
        capri_aa = capripox_residues[0]
        if capri_aa == '-': 
            continue
    else:
        continue
    if all(res == '-' for res in non_capripox_residues):  
        type2_sites.append((i + 1, capri_aa))
    elif all(res != '-' for res in non_capripox_residues) and all(res != capri_aa for res in non_capripox_residues): 
        non_capri_residues = set(non_capripox_residues)
        type1_sites.append((i + 1, capri_aa, non_capri_residues))

print("1: Positions where capripoxviruses have a conserved amino acid, different from all others")
for pos, capri_aa, other_aa in type1_sites:
    print(f"Position {pos}: Capripox = {capri_aa}, Others = {other_aa}")

print("2: Positions where capripoxviruses have a conserved amino acid, all others have gaps")
for pos, capri_aa in type2_sites:
    print(f"Position {pos}: Capripox = {capri_aa}, Others = gap")

# Output summary statistics
print(f"\n\nTotal 1 sites: {len(type1_sites)}")
print(f"Total 2 sites: {len(type2_sites)}")