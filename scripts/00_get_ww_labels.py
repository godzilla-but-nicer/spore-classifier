import pandas as pd
import re
from Bio import SeqIO
from glob import glob

# get the list of genomes and labels
genomes = pd.read_csv('data/labelled/ww_genomes.txt')
label_sheet = pd.read_csv('data/labelled/weller_and_wu_sporulators.csv')

# for each genome we need to try and get the taxa name from the faa (have to try a few)
# and then look that up in the excel file. then we will add the label to the list
labels = []
taxon_names = []
genome_ids = []
for gid in genomes['genome']:
    potential_names = []

    # grab unique names from the fasta
    for idx, gene in enumerate(SeqIO.parse(f"data/labelled/ww_proteins/{gid}.faa", format='fasta')):
        # pull out the name in brackets
        if re.search(r'\[([a-zA-Z0-9_ ]+)\]', gene.description):
            pot_name = re.findall(r'\[([a-zA-Z0-9_ ]+)\]', gene.description)[0]

            # if its a new name add it to the list
            if pot_name not in potential_names:
                potential_names.append(pot_name)

    # now we need to check all of the unique names against the excel file
    for name in potential_names:
        # if we find a match we can grab the label and be done
        if name in label_sheet['Taxon Name'].values:
            print('found one')
            raw_label = label_sheet[label_sheet['Taxon Name'] == name]['Spore-Forming?'].values[0]
            # we want integers for our labels
            labels.append(int(raw_label == 'Y'))
            taxon_names.append(name)
            genome_ids.append(gid)
            break

out_df = pd.DataFrame({'genome': genome_ids,
                       'taxon_name': taxon_names,
                       'label': labels})
out_df.to_csv('data/labelled/ww_genomes_labeled.csv')


