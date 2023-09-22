import sys
import editdistance
import subprocess
import pandas as pd
from sklearn.metrics import (accuracy_score, f1_score, precision_score,
                             recall_score, roc_auc_score)

from scripts.globals import *


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


def find_most_different_string(original_string, column, cache=None):
    if cache is None:
        cache = {}
    max_distance = float('-inf')
    most_different_string = ''
    for sequence in column:
        if sequence in cache:
            distance = cache[sequence]
        else:
            distance = editdistance.eval(original_string, sequence)
            cache[sequence] = distance
        if distance > max_distance:
            max_distance = distance
            most_different_string = sequence
    return most_different_string


def invoke_rnaduplex(long_sequence, short_sequence, energy_range = 5.0,
                     rnaduplex_location = RNADUPLEX_LOCATION):

    input_sequence = f"{long_sequence}\n{short_sequence}".encode()

    rnaduplex_subprocess = subprocess.Popen(
        [rnaduplex_location, "-e", f"{energy_range}", "-s"],
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    output, error = rnaduplex_subprocess.communicate(input=input_sequence)
    rnaduplex_subprocess.wait()

    first_line = output.decode().split("\n")[0].split()

    dot_bracket_long, dot_bracket_short = first_line[0].split("&")
    start_long, end_long = map(int, first_line[1].split(","))
    start_short, end_short = map(int, first_line[3].split(","))
    energy = float(first_line[-1].strip("()"))

# -1's here convert biological coordinates into 0-index coordinates
    return start_long-1, end_long-1, dot_bracket_long, start_short-1, end_short-1, dot_bracket_short, energy

def find_matches_with_rnaduplex(df):

    mrna_starts = []
    mrna_ends = []
    mirna_starts = []
    mirna_ends = []
    mirna_dot_brackets = []
    energies = []

    for _, row in df.iterrows():
        start_long, end_long, dot_bracket_long, start_short, end_short, dot_bracket_short, energy = invoke_rnaduplex(
            row.mrna_sequence, row.mirna_sequence)

        mrna_starts.append(start_long)
        mrna_ends.append(end_long)
        mirna_starts.append(start_short)
        mirna_ends.append(end_short)
        mirna_dot_brackets.append(dot_bracket_short)
        energies.append(energy)

    df = pd.DataFrame({"mrna_start": mrna_starts,
                       "mrna_end": mrna_ends,
                       "pred_energy": energies,
                       "mirna_start": mirna_starts,
                       "mirna_end": mirna_ends,
                       "mirna_dot_bracket_5to3": mirna_dot_brackets,
                       })

    return df


def report_performance(model, X, y):
    # Make predictions on the input data
    y_pred = model.predict(X)

    # Calculate performance metrics
    accuracy = accuracy_score(y, y_pred)
    precision = precision_score(y, y_pred)
    recall = recall_score(y, y_pred)
    f1 = f1_score(y, y_pred)
    roc_auc = roc_auc_score(y, y_pred)

    return {
        'Accuracy': accuracy,
        'Precision': precision,
        'Recall': recall,
        'F1-Score': f1,
        'ROC AUC': roc_auc,
    }
    

def get_nucleotides_in_interval(chrom, start, end):
    # sourcery skip: extract-method
    file_path = f"{GRCH37_DIR}/Homo_sapiens.GRCh37.dna.chromosome.{chrom}.fa"
    with open(file_path, 'r') as file:
        file.readline()
        byte_position = file.tell()
        line_length = len(file.readline().strip())
        start_offset = start - 1
        end_offset = end - 1
        num_start_new_lines = start_offset // line_length
        num_end_new_lines = end_offset // line_length
        start_byte_position = byte_position + start_offset + num_start_new_lines
        end_byte_position = byte_position + end_offset + num_end_new_lines
        file.seek(start_byte_position)

        # Read the nucleotides in the interval
        nucleotides = file.read(end_byte_position - start_byte_position + 1)

    # Remove newlines from the nucleotides
    nucleotides = nucleotides.replace('\n', '')

    return nucleotides

def get_nucleotide_at_position(chrom, position):
    file_path = f"{GRCH37_DIR}/Homo_sapiens.GRCh37.dna.chromosome.{chrom}.fa"
    with open(file_path, 'r') as file:
        file.readline()
        byte_position = file.tell()
        line_length = len(file.readline().strip())
        offset = position - 1
        num_new_lines = offset // line_length
        byte_position = byte_position + offset + num_new_lines
        file.seek(byte_position)

        # Read the nucleotide at the position
        nucleotide = file.read(1)
    return nucleotide