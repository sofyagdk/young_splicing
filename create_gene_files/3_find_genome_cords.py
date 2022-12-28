import pandas as pd
import sys
sys.path.append('/home/sofya/main/Lab/introns/young_splicing/splib/')
from maf_parser import iterator_over_maf_blocks

gene_table = pd.read_csv('../../data/gene_files/gene_cords_table.csv', sep=' ', header=0)
gene_table['kb_before_start'] = gene_table['start'] - 1000
gene_table['kb_after_end'] = gene_table['end'] + 1000

chr_genes = gene_table[gene_table['chr'] == '1'].reset_index(drop=True)

output_folder = '../../data/gene_files/chr1/'


species_list = ['hg38', 'calJac3', 'otoGar3', 'rheMac8', 'micMur3']
seqs_column = {sp: list() for sp in species_list}
starts_column = {sp: list() for sp in species_list}
ends_column = {sp: list() for sp in species_list}

for index, gene_id in enumerate(chr_genes['gene_id'].values):
    filename = output_folder + gene_id + '.txt'
    sp_to_genome_cords = {sp: dict() for sp in species_list}

    with open(filename, 'r') as fin:
        line = fin.readline()
        while line[:2] != '$2':
            line = fin.readline()

        for block in iterator_over_maf_blocks(fin):
            for maf in block:
                s, e = maf._start, maf._end
                sp, source_seq = maf.sp, maf.source_seq
                if len(sp_to_genome_cords[sp]) == 0:
                    sp_to_genome_cords[sp][source_seq] = (s, e)
                else:
                    if source_seq in sp_to_genome_cords[sp]:
                        cs, ce = sp_to_genome_cords[sp][source_seq]
                        sp_to_genome_cords[sp][source_seq] = (min(s, cs), max(e, ce))

    for sp in species_list:
        seqs_column[sp].append(','.join(sp_to_genome_cords[sp].keys()))
        starts_column[sp].append(','.join(map(str, [region[0] for region in sp_to_genome_cords[sp].values()])))
        ends_column[sp].append(','.join(map(str, [region[1] for region in sp_to_genome_cords[sp].values()])))

for sp in species_list:
    chr_genes[sp + '_source'] = seqs_column[sp]
    chr_genes[sp + '_starts'] = starts_column[sp]
    chr_genes[sp + '_ends'] = ends_column[sp]

chr_genes.to_csv('../../data/gene_files/chr1_gene_cords_table.csv')
