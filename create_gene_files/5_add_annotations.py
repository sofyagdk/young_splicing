import pandas as pd
import sys
sys.path.append('/home/sofya/main/Lab/introns/young_splicing/splib/')
from gff3_parser import GFF3Line, iterator_over_genes_with_annotation


added_genomes = pd.read_csv('added_genomes_with_chr.txt', sep=' ', header=0)
added_genomes['human_chr_name'] = added_genomes['human_chr_name'].astype('str')


annotations = ['Homo_sapiens.GRCh38.104.chr.gff3',
               'Callithrix_jacchus.C_jacchus3.2.1.91.chr.gff3',
               'Macaca_mulatta.Mmul_8.0.1.97.chr.gff3',
               'Microcebus_murinus.Mmur_3.0.104.chr.gff3',
               'Otolemur_garnettii.OtoGar3.104.gff3']

species_list = ['hg38', 'calJac3', 'rheMac8', 'micMur3', 'otoGar3']

chr_names = ['1', '2', '3', '4', '5', '6', '7', '8', '9',
             '10', '11', '12', '13', '14', '15', '16',
             '17', '18', '19', '20', '21', '22', 'X', 'Y']

for chr_name in chr_names:

    chr_genes = pd.read_csv(f'../../data/gene_files/chr{chr_name}_gene_cords_table.csv', header=0, sep=' ')
    for human_gene_id in chr_genes['gene_id']:
        with open(f'../../data/gene_files/chr{chr_name}/' + human_gene_id + '.txt', 'a') as fout:
            fout.write('$4 species annotations\n')

for index, sp in enumerate(species_list):
        print(sp)

        sp_table = added_genomes[added_genomes['sp'] == sp]
        filename = '../../data/annotations/' + annotations[index]

        with open(filename, 'r') as fin:
            for annotation_block in iterator_over_genes_with_annotation(fin, genetypes=['gene', 'ncRNA_gene']):
                gene = annotation_block[0]
                chr_num, start, end = gene.seq_id, gene.start, gene.end

                rows = sp_table[sp_table['chr_num'] == chr_num]
                rows = rows[rows['start'] < end]
                rows = rows[rows['end'] > start]

                for index, row in rows.iterrows():
                    human_gene_id = row['human_gene_id']
                    human_chr_name = row['human_chr_name']

                    if str(human_chr_name) != '3':
                        continue

                    with open(f'../../data/gene_files/chr{human_chr_name}/{human_gene_id}.txt', 'a') as fout:
                        fout.write(f'>>{sp} {chr_num} {start} {end}\n')
                        for gff3_line in annotation_block:
                            fout.write(str(gff3_line))








otogar_alias = pd.read_csv('../../data/genomes/otogar_alias.txt', sep = '\t')
ucsc_to_refseq = {pair[1]:pair[0] for pair in otogar_alias[['refseq', 'ucsc']].values}

sp = 'otoGar3'



sp_table = added_genomes[added_genomes['sp'] == sp]
sp_table['chr_num'] = sp_table['chr_num'].astype(str)
filename = '../../data/annotations/' + annotations[4]

with open(filename, 'r') as fin:

    for annotation_block in iterator_over_genes_with_annotation(fin, genetypes=['gene', 'ncRNA_gene']):
        gene = annotation_block[0]
        chr_num, start, end = gene.seq_id, gene.start, gene.end

        rows = sp_table[sp_table['chr_num'] == str(chr_num)]
        print(chr_num)


        rows = rows[rows['start'] < end]
        rows = rows[rows['end'] > start]

        for index, row in rows.iterrows():
            human_gene_id = row['human_gene_id']
            human_chr_name = row['human_chr_name']

            with open(f'../../data/gene_files/chr{human_chr_name}/{human_gene_id}.txt', 'a') as fout:
                fout.write(f'>>{sp} {chr_num} {start} {end}\n')
                for gff3_line in annotation_block:
                    fout.write(str(gff3_line))
