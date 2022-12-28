import sys
sys.path.append('/home/sofya/main/Lab/introns/young_splicing/splib/')
import time
from collections import defaultdict

from gff3_parser import GFF3Line, iterator_over_genes_with_annotation
from maf_parser import iterator_over_maf_blocks, SMafLine
from fasta_parser import iterator_over_fasta



class GeneFile:
    def __init__(self, filename):
        with open(filename, 'r') as fin:
            lines = fin.readlines()

        self.human_center_annot = list()
        self.maf_blocks = list()
        self.genome_seqs = dict()  # species: chr start end
        self.species_annotations = defaultdict(list)  # species: list[annot_blocks]
        self.alignment = dict()  # species: chr

        block = list()
        current_species = str()

        for line in lines:
            if line[:2] == '$2':
                print(block)
                self.human_center_annot = list(iterator_over_genes_with_annotation(block))[0]
                block = list()
            elif line[:2] == '$3':
                self.maf_blocks = list(iterator_over_maf_blocks(block))
                block = list()
            elif line[:2] == '$4':
                for seq_id, info, seq in iterator_over_fasta(block, decode=False ):
                    self.genome_seqs[seq_id] = [info, seq]
                block = list()
            elif line[:2] == '>>' or line[:2] == '$5':
                if current_species:
                    self.species_annotations[current_species].\
                        append(list(iterator_over_genes_with_annotation(block))[0])
                current_species = line[2:].split()[0]
                block = list()
            elif line[:2] == '$1':
                continue
            else:
                block.append(line)
