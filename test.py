dna_seq = 'GCGTTTGACCGCGCTTGGGTGGCCTGGGACCCTGTGGGAGGCTTCCCCGGCGCCGAGAGCCCTGGCTGACGGCTGATGGGGAGGAGCCGGCGGGCGGAGAAGGCCACGGGCTCCCCAGTACCCTCACCTGCGCGGGATCGCTGCGGGAAACCAGGGGGAGCTTCGGCAGGGCCTGCAGAGAGGACAAGCGAAGTTAAGAGCCTAGTGTACTTGCCGCTGGGAGCTGGGCTAGGCCCCCAACCTTTGCCCTGAAGATGCTGGCAGAGCAGGATGTTGTAACGGGAAATGTCAGAAATACTGCAAGCAAACTGAAAACAACCCATCCATGTAGGAAAGAATAACACGG'

protein_seq = 'MGRSRRAEKATGSPVPSPARDRCGKPGGASAGPAERTSEVKSLVYLPLGAGLGPQPLP_'


def translate_DNA(sequence):
    start_codon = 'ATG'
    stop_codons = ('TAA', 'TAG', 'TGA')
    codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'
    }

    start = sequence.find(start_codon)
    # Create a codon list to store codons generated from coding seq.
    codons = []

    for i in range(start, len(sequence), 3):
        codon = sequence[i:i + 3]
        print(codon)
        if (codon in stop_codons):
            print("stop found")
            break
        else:
            codons.append(codon)

      # Finish the loop or write your own code to find the coding sequence which should
      # start at the first start codon and stop at the occurrence of any stop
      # codon.
    # Translate condons to protein seq.
    protein_sequence = ''.join([codontable[codon] for codon in codons])
    return "{0}_".format(protein_sequence)

# print(protein_seq)
print(translate_DNA(dna_seq))
