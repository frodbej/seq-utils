from itertools import product

# some useful functions to deal with some basics operations on biological sequences

def reverse_complement(seq):
    """Take as input a sequence and return the reverse complement sequence."""

    comp_base = {
        'A':'T',
        'T':'A',
        'C':'G',
        'G':'C'
    }

    return ''.join([comp_base[nt] for nt in seq[::-1]])


def GCCont(seq):
    """Return GC content of a sequence."""

    count, l = 0, len(seq)
    for i in range(0,l):
        if seq[i] in "GC":
            count += 1
    return count/l


def maxHomopolymerLength(seq):
    """Return the length of the longest homopolymer in a sequence."""

    lMaxHom, lCurHom, l = 1, 1, len(seq)
    for i in range(1,l):
        if seq[i] == seq[i-1]:
            lCurHom += 1
        else:
            if lCurHom > lMaxHom:
                lMaxHom = lCurHom
            lCurHom = 1
    return lMaxHom


def degenToSet(degen):
    """Take as input a degenerate sequence and return a list containing the corresponding non-degenerate sequences."""

    conversion_table = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']}

    return list(map("".join, product(*map(conversion_table.get, degen))))


def loadFasta(filename):
    """Open fasta file and return two lists: one for headers and one for sequences"""

    with open(filename, 'r') as f:
        # split at headers
        data = f.read().split(">")
        data.pop(0)
        headers, sequences = [], []
        for sequence in data:
            lines = sequence.split('\n')
            headers.append(lines.pop(0))
            sequences.append(''.join(lines))
    return (headers, sequences)