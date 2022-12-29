import pandas as pd
import gzip
import typing as tp
import sys
sys.path.append('/home/sofya/main/Lab/introns/young_splicing/splib/')
from maf_parser import iterator_over_maf_blocks, SMafLine


def find_block_with_position(positions, block_cords):
    block_pointer = 0
    list_of_blocks = list()

    start_from = 0
    end_with = -1

    while positions[start_from] < block_cords[0][0]:
        list_of_blocks.append(None)
        start_from += 1

    while positions[end_with] >= block_cords[-1][1]:
        end_with -= 1

    for pos in positions[start_from:len(positions) + end_with + 1]:
        while pos < block_cords[block_pointer][0] and block_pointer > 0:
            block_pointer -= 1

        while pos >= block_cords[block_pointer][1]:
            block_pointer += 1

        if pos >= block_cords[block_pointer][0]:
            list_of_blocks.append(block_pointer)
        else:
            list_of_blocks.append(block_pointer + 1)

        assert pos < block_cords[block_pointer + 1][0], f'suggested {pos} in {block_cords[block_pointer]} \
        and between {block_cords[block_pointer - 1]} and {block_cords[block_pointer + 1]}. Please check if sorted ' \
                                                        f'lists are given'
        assert pos >= block_cords[block_pointer - 1][1], f'suggested {pos} in {block_cords[block_pointer]}\
        and between {block_cords[block_pointer - 1]} and {block_cords[block_pointer + 1]}. Please check if sorted ' \
                                                         f'lists are given'

    while end_with != -1:
        list_of_blocks.append(None)
        end_with += 1

    return list_of_blocks

gene_table = pd.read_csv('../../data/gene_files/gene_cords_table.csv', sep=' ', header=0)
gene_table['kb_before_start'] = gene_table['start'] - 1000
gene_table['kb_after_end'] = gene_table['end'] + 1000


chr_names = ['1', '2', '3', '4', '5', '6', '7', '8', '9',
             '10', '11', '12', '13', '14', '15', '16',
             '17', '18', '19', '20', '21', '22', 'X', 'Y']


for chr_name in chr_names:
    print(chr_name)
    if chr_name != '3':
        continue

    maf_human_cords: list[tuple[int, int]] = list()  # human_start, human_end
    mafblock_list: list[list[tp.Any]] = list()

    maf_file = gzip.open(f'../../data/mafs/chr{chr_name}.maf.gz', 'r')
    for block in iterator_over_maf_blocks(maf_file, decode=True):
        mafblock_list.append(block)
        maf_human_cords.append((block[0]._start, block[0]._end))
    maf_file.close()

    print('finished reading maf blocks')

    chr_genes = gene_table[gene_table['chr'] == int(chr_name)]

    start_containing_blocks = find_block_with_position(chr_genes['kb_before_start'].values, maf_human_cords)
    end_containing_blocks = find_block_with_position(chr_genes['kb_after_end'].values, maf_human_cords)

    output_folder = f'../../data/gene_files/chr{chr_name}/'

    for index, gene_id in enumerate(chr_genes['gene_id'].values):
        print(index, gene_id)

        with open(output_folder + gene_id + '.txt', 'a') as fout:
            s_block, e_block = start_containing_blocks[index], end_containing_blocks[index]
            if s_block is None or e_block is None:
                fout.write("$2 Error: gene didn't fully get into alignment\n")
                continue

            fout.write("$2 maf alignment blocks\n")
            for block in mafblock_list[max(s_block - 1, 0): min(e_block + 1, len(mafblock_list))]:
                for maf in block:
                    fout.write(str(maf))
                    if maf.i_line is not None:
                        fout.write(maf.i_line)
                fout.write('a\n')
