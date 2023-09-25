import csv
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
import xgboost as xgb
from editdistance import eval
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
            distance = eval(original_string, sequence)
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


############### 
# multithreaded

def prepare_jobs_from_df(df):
    
    mirna_df = pd.read_csv(MIRNA_CSV)
    
    wt_sequences = df["sequence"].to_list()
    mutated_sequences = df["mutated_sequence"].to_list()
    identifiers = df["id"].to_list()
    mirna_identifiers = mirna_df["mirna_accession"].to_list()
    mirna_sequences = (mirna_df["sequence"]
                       .tolist())

    wt_pairs = [(i, j) for i in wt_sequences for j in mirna_sequences]
    mutated_pairs = [(k, l)
                     for k in mutated_sequences for l in mirna_sequences]
    id_pairs = [(m, n) for m in identifiers for n in mirna_identifiers]

    wt_jobs_set = {
        wt_pair + id_pair for wt_pair, id_pair in zip(wt_pairs, id_pairs)
    }
    mutated_jobs_set = {
        mutated_pair + id_pair
        for mutated_pair, id_pair in zip(mutated_pairs, id_pairs)
    }

    return wt_jobs_set, mutated_jobs_set


def run_rnaduplex_multithreaded(long_sequence, short_sequence, long_identifier, short_identifier):

    # handle the sequence input
    input_sequence = f"{long_sequence}\n{short_sequence}"

    result = subprocess.run(
        [RNADUPLEX_LOCATION, "-e", "5.0", "-s"],
        input=input_sequence,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,  # Capture stderr
        text=True)

    try:
        lines = result.stdout.split("\n")
        first_line = lines[0].split()
        brackets, long_start_end, _, short_start_end, energy = first_line
        long_bracket, short_bracket = brackets.split("&")
        start_long, end_long = map(int, long_start_end.split(","))
        start_short, end_short = map(int, short_start_end.split(","))
        energy_value = float(energy.strip("()"))

        return start_long - 1, end_long, long_bracket, start_short - 1, end_short, short_bracket, energy_value, long_identifier, short_identifier, long_sequence, short_sequence

    # predictions with positive energy prints out "( 5.0163)" with an extra whitespace in the beginning, so strip function adds another element "(" to the array.
    # we try to capture 6 elements with 5 assignments so the interpreter returns ValueError: too many values to unpack (expected 5). This part handles this issue.
    # predictions with negative energy prints out "(-5.0163)" so it works without problem.
    except ValueError as e:
        lines = result.stdout.split("\n")
        first_line = lines[0].split()

        if first_line == []:
            return 0, 1, ".", 0, 1, ".", 0, long_identifier, short_identifier, long_sequence, short_sequence
        
        brackets, long_start_end, _, short_start_end, _, energy = first_line
        long_bracket, short_bracket = brackets.split("&")
        start_long, end_long = map(int, long_start_end.split(","))
        start_short, end_short = map(int, short_start_end.split(","))
        energy_value = float(energy.strip("()"))

        return start_long - 1, end_long, long_bracket, start_short - 1, end_short, short_bracket, energy_value, long_identifier, short_identifier, long_sequence, short_sequence


def run_jobs_multithreaded(job_list, binary_value, result_file):

    with ProcessPoolExecutor() as executor:
        future_jobs = [executor.submit(
            run_rnaduplex_multithreaded, *job) for job in job_list]

        results = []

        for future in as_completed(future_jobs):
            result = future.result()
            result_with_binary = result + (binary_value,)
            results.append(result_with_binary)

        with open(result_file, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(results)
    
    return results

def handle_target_file(result_file):
    if os.path.isfile(result_file):
        # Get the current timestamp
        timestamp = time.strftime("%Y%m%d%H%M%S")
        # Rename the existing CSV file with a ".bak" suffix and a timestamp
        bak_filename = f"{result_file}_{timestamp}.bak"
        os.rename(result_file, bak_filename)
        

def split_df_to_num_thread_chunks(df):
    chunk_size = (len(df) + NUM_CORES - 1) // NUM_CORES
    return [df.iloc[i:i+chunk_size] for i in range(0, len(df), chunk_size)]


##
# xgb functions


def make_predictions_regressor(df, df_filtered, model_name=XGB_MODEL):

    # Create an empty DMatrix model
    dmatrix_model = xgb.DMatrix(df_filtered)

    # Load the pre-trained model
    model = xgb.Booster()
    model.load_model(model_name)

    # Make predictions on the df_filtered DataFrame
    predictions = model.predict(dmatrix_model)

    # Append the predictions to the "prediction" column in the df DataFrame
    df["prediction"] = predictions
    df["binary_prediction"] = (df["prediction"] > 0.5).astype(int)
    df.sort_values(["id", "is_mutated"], inplace=True)

    # creating wt and mut dfs
    wt = df[df.is_mutated == 0].reset_index(drop=True)
    mut = df[df.is_mutated == 1].reset_index(drop=True)

    # Calculate the difference between wt and mut predictions
    wt['pred_difference'] = mut['prediction'] - wt['prediction']
    wt['pred_difference_binary'] = mut['binary_prediction'] - \
        wt['binary_prediction']

    # Merge the difference values back to the original df DataFrame
    df = df.merge(
        wt[['id', 'pred_difference', 'pred_difference_binary']], on='id', how='left')

    return df

def filter_columns_for_xgb_prediction(df):
    cols_to_keep = [
        "pred_energy",
        "pred_num_basepairs",
        "pred_seed_basepairs",
        "ta_log10",
        "sps_mean",
        "anchor_a",
        "6mer_seed",
        "match_8",
        "6mer_seed_1_mismatch",
        "compensatory_site",
        "supplementary_site",
        "supplementary_site_2",
        "empty_seed",
        "9_consecutive_match_anywhere",
        "mirna_conservation",
        "seed_8mer",
        "seed_7mer_a1",
        "seed_7mer_m8",
        "seed_compensatory",
        "seed_clash_2",
        "seed_clash_3",
        "seed_clash_4",
        "seed_clash_5",
        "mre_au_content",
        "local_au_content"
    ]
    return df[cols_to_keep]


