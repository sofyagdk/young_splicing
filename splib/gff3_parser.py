import typing as tp

# 1	ensembl	exon	198729	198866	.	+	.	Parent=transcript:ENSMICT00000073141;Name=ENSMICE00000419148;constitutive=0;ensembl_end_phase=0;ensembl_phase=0;exon_id=ENSMICE00000419148;rank=2;version=1

class GFF3Line:
    def __init__(self, line: str):
        separated = line.strip().split('\t')
        self.seq_id: str = separated[0]
        if '.' in self.seq_id:
            self.seq_id = self.seq_id.split('.')[0]
        self.source: str = separated[1]
        self.type: str = separated[2]
        self.start: int = int(separated[3])
        self.end: int = int(separated[4])
        self.score: tp.Optional[float] = None if separated[5] == '.' else float(separated[5])

        self.strand: tp.Optional[int] = None
        if separated[6] == '+':
            self.strand = 1
        elif separated[6] == '-':
            self.strand = -1

        self.phase: tp.Optional[int] = None if separated[7] == '.' else int(separated[7])

        self.attributes = dict()
        for pair in separated[8].split(';'):
            div_pair = pair.split('=')
            self.attributes[div_pair[0]] = div_pair[1]

        self.initial_line = line

    def __repr__(self) -> str:
        return self.initial_line


# gene
# ncRNA_gene
# pseudogene

def iterator_over_genes_with_annotation(line_iterator: tp.Iterator,
                                        genetypes: list = ['gene', 'ncRNA_gene']) -> list:
    # with open(filename, 'r') as fin:

        current_annotation_block: list[GFF3Line] = list()

        for line in line_iterator:
            if line.startswith('###') and len(current_annotation_block) > 0:
                yield current_annotation_block
                current_annotation_block = list()
                continue

            if line[0] == '#':
                continue
            gff3_line = GFF3Line(line)
            if gff3_line.type in genetypes:
                current_annotation_block.append(gff3_line)
            else:
                if len(current_annotation_block) > 0:
                    current_annotation_block.append(gff3_line)
        if len(current_annotation_block) > 0:
            yield current_annotation_block
