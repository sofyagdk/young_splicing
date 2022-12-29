import typing as tp
import gzip


class MafLine:
    def __init__(self, line):
        separated = line.strip().split()
        self.type = separated[0]
        self.maf_line = line
        if self.type != 'a':
            self.sp, self.source_seq = separated[1].split('.')

    def __repr__(self):
        return self.maf_line


class SMafLine(MafLine):

    def __init__(self, line):
        super().__init__(line)
        separated = self.maf_line.strip().split()
        self.start = int(separated[2])
        self.length = int(separated[3])
        self.ori = separated[4]
        self.source_length = int(separated[5])
        self.seq = separated[6]

        self._start = self.start
        self._end = self.start + self.length

        if self.ori == '-':
            self._end = self.source_length - self.start
            self._start = self.source_length - self.start - self.length

        self.leftStatus = None
        self.rightStatus = None
        self.leftCount = None
        self.rightCount = None
        self.i_line = ''

    def add_i_line(self, i_line):
        separated = i_line.strip().split()
        assert [self.sp, self.source_seq] == separated[1].split('.'), 'The passed maf i line is for a different species'
        self.i_line = i_line
        """
        A character that specifies the relationship between the sequence in this block
        and the sequence that appears in the previous block
        
        C -- the sequence before or after is contiguous with this block.
        I -- there are bases between the bases in this block and the one before or after it.
        N -- this is the first sequence from this src chrom or scaffold.
        n -- this is the first sequence from this src chrom or scaffold but it is bridged by another alignment from a 
             different chrom or scaffold.
        M -- there is missing data before or after this block (Ns in the sequence).
        T -- the sequence in this block has been used before in a previous block (likely a tandem duplication)
        """
        self.leftStatus = str(separated[2])
        self.rightStatus = str(separated[4])

        """
        Usually the number of bases in the aligning species between the start of
        this alignment and the end of the previous one.
        """
        self.leftCount = int(separated[3])
        self.rightCount = str(separated[5])


# species_list = ['hg38', 'calJac3', 'otoGar3', 'rheMac8', 'micMur3']




def iterator_over_maf_blocks(fin,
                             species_list=['hg38', 'calJac3', 'otoGar3', 'rheMac8', 'micMur3'],
                             decode: bool = False):

# def iterator_over_maf_blocks(filename='../../data/mafs/chr1.maf.gz',
#                              species_list=['hg38', 'calJac3', 'otoGar3', 'rheMac8', 'micMur3']):
#
#     with gzip.open(filename) as fin:
        current_block: list[SMafLine] = list()
        for enc_line in fin:
            if decode:
                line = enc_line.decode('utf8')
            else:
                line = enc_line
            if line[0] == '#' or line == '\n':
                continue

            if line[0] == '$': #### for gene files!
                break

            if line[0] == 'a':
                if len(current_block) != 0:
                    yield current_block
                current_block = list()
                continue

            sp = line[2:9].split('.')[0]
            if sp not in species_list:
                continue

            if line[0] == 's':
                current_block.append(SMafLine(line))
            elif line[0] == 'i':
                current_block[-1].add_i_line(line)
        if len(current_block) > 0:
            yield current_block
