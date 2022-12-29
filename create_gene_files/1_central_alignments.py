import time
import sys
sys.path.append('/home/sofya/main/Lab/introns/young_splicing/splib/')
from gff3_parser import iterator_over_genes_with_annotation


human_annotation = '../../data/annotations/Homo_sapiens.GRCh38.104.chr.gff3'
chr_names = ['1', '2', '3', '4', '5', '6', '7', '8', '9',
             '10', '11', '12', '13', '14', '15', '16',
             '17', '18', '19', '20', '21', '22', 'X', 'Y']

gene_table = open('../../data/gene_files/gene_cords_table.csv', 'w')
gene_table.write('chr start end gene_id\n')

s = time.time()

with open(human_annotation, 'r') as fin:
    for region in iterator_over_genes_with_annotation(fin):
        gene = region[0]

        output_folder = f'../../data/gene_files/chr{gene.seq_id}/' if gene.seq_id in chr_names \
            else '../../data/gene_files/unplaced/'

        with open(output_folder + gene.attributes['gene_id'] + '.txt', 'w') as fout:
            fout.write(f'$1 human_gene cords=chr{gene.seq_id}:{gene.start}-{gene.end}\n')
            gene_table.write(f'{gene.seq_id} {gene.start} {gene.end} {gene.attributes["gene_id"]}\n')
            for annot in region:
                fout.write(str(annot))

print(f'Success, time={time.time() - s}')
