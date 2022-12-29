
def iterator_over_fasta(fastafin, decode=False):
    seq = str()
    seq_id = str()
    info = str()


    for enc_line in fastafin:
        if decode:
            line = enc_line.decode('utf8')
        else:
            line = enc_line

        if line[0] == '>':
            if seq:
                yield seq_id, info, seq
            seq = str()
            seq_id = line[1:].strip().split()[0]
            info = line[1:].strip()
        else:
            seq += line.strip()
    if seq:
        yield seq_id, info, seq


def read_fasta_to_dict(filename):
    fasta = dict()
    ident, seq = str(), str()
    with open(filename, 'r') as fin:
        for line in fin:
            if line[0] == '>':
                if ident != '':
                    fasta[ident] = seq
                ident = line[1:].strip().split()[0]
                seq = ''
            else:
                seq += line.strip()
        if ident != '':
            fasta[ident] = seq
    return fasta
