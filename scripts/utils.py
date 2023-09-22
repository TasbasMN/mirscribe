import sys


def get_size_in_ram(obj):
    size_in_bytes = sys.getsizeof(obj)
    size_in_mb = size_in_bytes / (1024 ** 3)
    print(f"The object uses approximately {size_in_mb} GB of RAM.")




def reverse_complement_rna_to_dna(rna_sequence):
    complement = {'A': 'T', 'U': 'A', 'C': 'G', 'G': 'C'}
    reverse_seq = rna_sequence[::-1]
    return ''.join(complement[base] for base in reverse_seq)

def reverse_complement_dna_to_rna(rna_sequence):
    complement = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_seq = rna_sequence[::-1]
    return ''.join(complement[base] for base in reverse_seq)

def calculate_au_content(sequence):
    au_count = sequence.count(
        'A') + sequence.count('T') + sequence.count('U')
    return None if len(sequence) == 0 else au_count / len(sequence)


