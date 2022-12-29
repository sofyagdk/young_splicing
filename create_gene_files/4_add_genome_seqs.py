import pandas as pd
import gzip
import sys
sys.path.append('/home/sofya/main/Lab/introns/young_splicing/splib/')
from fasta_parser import iterator_over_fasta


chr_names = ['1', '2', '3', '4', '5', '6', '7', '8', '9',
             '10', '11', '12', '13', '14', '15', '16',
             '17', '18', '19', '20', '21', '22', 'X', 'Y']

species_list = ['hg38', 'calJac3', 'rheMac8', 'micMur3', 'otoGar3']
genomes = ['GCF_000001405.38_GRCh38.p12_genomic.fna.gz',
           'GCF_000004665.1_Callithrix_jacchus-3.2_genomic.fna.gz',
           'GCF_000772875.2_Mmul_8.0.1_genomic.fna.gz',
           'GCF_000165445.2_Mmur_3.0_genomic.fna.gz',
           'GCF_000181295.1_OtoGar3_genomic.fna.gz']

added_genomes = open('added_genomes_with_chr.txt', 'a')
added_genomes.write('human_gene_id sp chr_num start end human_chr_name\n')


def add_genome_sequence(sp, chr_num, seq, start, end, human_gene_id, human_chr_name):
    with open(f'../../data/gene_files/chr{human_chr_name}/{human_gene_id}.txt', 'a') as fout:
        fout.write(f'>{sp} chr{chr_num}:{start}-{end}\n')
        fout.write(seq[start:end] + '\n')
        added_genomes.write(f'{human_gene_id} {sp} {chr_num} {start} {end} {human_chr_name}\n')


for chr_name in chr_names:
    print(chr_name)

    chr_genes = pd.read_csv(f'../../data/gene_files/chr{chr_name}_gene_cords_table.csv', header=0, sep=' ')

    output_folder = f'../../data/gene_files/chr{chr_name}/'

    for human_gene_id in chr_genes['gene_id']:
        with open(f'../../data/gene_files/chr{chr_name}/' + human_gene_id + '.txt', 'a') as fout:
            fout.write('$3 genome sequences\n')

    for index, genome in enumerate(genomes[:4]):
        sp = species_list[index]
        print(sp)
        with gzip.open('../../data/genomes/' + genome) as fastafin:
            for seq_id, info, seq in iterator_over_fasta(fastafin, decode=True):
                if 'unlocalized genomic scaffold' in info:
                    continue
                chr_num = info.split(',')[0].split()[-1]
                print(chr_num)
                rows = chr_genes[chr_genes[sp + '_source'] == 'chr' + str(chr_num)]
                for index, row in rows.iterrows():
                    start, end, human_gene_id = int(row[sp + '_starts']), int(row[sp + '_ends']), str(row['gene_id'])
                    add_genome_sequence(sp, chr_num, seq, start, end, human_gene_id, chr_name)

    #
    otogar_alias = pd.read_csv('../../data/genomes/otogar_alias.txt', sep = '\t')
    refseq_to_ucsc = {pair[0]:pair[1] for pair in otogar_alias[['refseq', 'ucsc']].values}


    with gzip.open('../../data/genomes/' + genomes[4]) as fastafin:
        sp = 'otoGar3'
        print(sp)
        for seq_id, info, seq in iterator_over_fasta(fastafin, decode=True):
            ucsc = refseq_to_ucsc[seq_id]
            print(ucsc)
            rows = chr_genes[chr_genes[sp + '_source'] == ucsc]
            for index, row in rows.iterrows():
                start, end, human_gene_id = int(row[sp + '_starts']), int(row[sp + '_ends']), str(row['gene_id'])
                add_genome_sequence(sp, ucsc, seq, start, end, human_gene_id, chr_name)

added_genomes.close()
