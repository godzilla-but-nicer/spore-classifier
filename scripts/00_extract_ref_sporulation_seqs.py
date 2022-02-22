import pandas as pd
import re
from Bio import SeqIO, SeqRecord

# species names as they appear in uniprot-sprot
SPROT_B_SUBTILIS = "Bacillus subtilis (strain 168)"
SPROT_C_DIFFICILE = "Clostridioides difficile"

# load the list of sporulation genes for each
baci_df = pd.read_csv('data/reference_sequences/subtiwiki-sporulation-genes.csv')
clos_df = pd.read_csv('data/reference_sequences/fimlaid-sporulation-genes.csv')

# get just the gene names for each species in a set
baci_genes = set(baci_df['locus tag'].values)
clos_genes = set(clos_df['Name'].values)

# make some empty lists to grow of our reference genes
baci_sprot = []
clos_sprot = []

for idx, gene in enumerate(SeqIO.parse('data/reference_sequences/uniprot_sprot_bacteria.fasta', format="fasta")):
    gene_name_given = re.search(r"GN=([a-zA-Z0-9_]+) PE=", gene.description)
    if gene_name_given:
        gene_name = re.findall(r"GN=([a-zA-Z0-9_]+) PE=", gene.description)[0]
        if (SPROT_B_SUBTILIS in gene.description) and (gene_name in baci_genes):
            baci_sprot.append(gene)
            #print(f'B. sub {gene_name} found!')
        if (SPROT_C_DIFFICILE in gene.description) and (gene_name in clos_genes):
            clos_sprot.append(gene)
            print(f"C. diff {gene_name} found!")

baci_sprot.extend(clos_sprot)
SeqIO.write(baci_sprot, 'data/reference_sequences/sporulation_sprot.fasta', format='fasta')
