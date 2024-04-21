import gc, os, argparse, logging, subprocess, csv
from scripts.globals import *

import xgboost as xgb
import pandas as pd
import numpy as np
import psutil

from concurrent.futures import ProcessPoolExecutor, as_completed


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

def parse_cli_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="VCF file to analyze")
    parser.add_argument("-s", "--start", default=0, type=int, help="start index")
    parser.add_argument("-e", "--end", default=-1, type=int, help="end index")
    parser.add_argument("-o", "--output_dir", default="./results", help="Output directory for runtime file, defaults to ./results")
    return parser.parse_args()

def load_vcf_into_df(file_path):
    """
    Load a VCF (Variant Call Format) file into a pandas DataFrame.
    
    Parameters:
    - file_path (str): The full path to the VCF file to be loaded.
    
    Returns:
    - DataFrame: A pandas DataFrame containing the VCF file data, with columns for 'chr', 'pos', 'id', 'ref', and 'alt'.
    
    Raises:
    - FileNotFoundError: If the specified file does not exist.
    - ValueError: If the file is not in the expected VCF format.
    - Exception: For other unforeseen errors during file loading.
    """
    # Ensure the file path ends with '.vcf', indicating it is a VCF file
    if not file_path.endswith('.vcf'):
        raise ValueError("The file specified does not appear to be a VCF file. Please provide a valid .vcf file.")

    try:
        return pd.read_csv(
            file_path,
            sep="\t",
            header=None,
            names=["chr", "pos", "id", "ref", "alt"],
            comment='#',
        )

    except FileNotFoundError as e:
        # Raise a more descriptive error if the file is not found
        raise FileNotFoundError(f"The file at {file_path} was not found.") from e
    
    except pd.errors.EmptyDataError as e:
        # Handle the case where the file is empty or does not contain expected data
        raise ValueError(
            f"The file at {file_path} is empty or not in the expected VCF format."
        ) from e
    except Exception as e:
        # Catch-all for other exceptions, with a generic error message
        print(f"An error occurred while trying to read the file at {file_path}: {e}")
        raise

def augment_id(df):
    """Augment the 'id' column with additional information."""
    df['id'] += '_' + df['chr'].astype(str) + '_' + df['pos'].astype(str) + '_' + df['ref'] + '_' + df['alt']


def validate_ref_nucleotides(df, output_file="invalid_rows.csv"):
    """
    Check if the 'ref' column matches the 'nucleotide_at_position' column, fetched from FASTA, in a DataFrame.
    Write invalid rows to a file and return the valid rows for downstream analysis.

    The 'nucleotide_at_position' column is populated as follows:
    - For rows where the reference length (ref_len) is greater than 1, the get_nucleotides_in_interval function
      is used to fetch the nucleotides in the interval [pos, pos + ref_len - 1] for the given chromosome (chr).
    - For rows where the reference length (ref_len) is 1, the get_nucleotide_at_position function is used to
      fetch the nucleotide at the given position (pos) for the given chromosome (chr).
    - For all other cases, an empty string is assigned.

    Args:
        df (pandas.DataFrame): The input DataFrame.
        output_file (str, optional): The file name to write invalid rows to. Default is "invalid_rows.csv".

    Returns:
        pandas.DataFrame: The DataFrame containing valid rows with matching 'ref' and 'nucleotide_at_position',
                           without the 'nucleotide_at_position' column.
    """
    # Add ref_len and alt_len columns
    df["ref_len"] = df["ref"].apply(len)
    df["alt_len"] = df["alt"].apply(len)

    # For rows where the reference length (ref_len) is greater than 1:
    # - Use the get_nucleotides_in_interval function to fetch the nucleotides
    #   in the interval [pos, pos + ref_len - 1] for the given chromosome (chr)
    df['nuc_at_pos'] = np.where(
        df['ref_len'] > 1,
        df.apply(lambda x: get_nucleotides_in_interval(x['chr'], x['pos'], x["pos"] + x["ref_len"] - 1), axis=1),
        # For rows where the reference length (ref_len) is 1:
        # - Use the get_nucleotide_at_position function to fetch the nucleotide
        #   at the given position (pos) for the given chromosome (chr)
        np.where(
            df['ref_len'] == 1,
            df.apply(lambda x: get_nucleotide_at_position(x['chr'], x['pos']), axis=1),
            # For all other cases, set the value to an empty string
            ""
        )
    )

    ## Check if ref matches nucleotide_at_position
    mask = df['ref'] != df['nuc_at_pos']

    ## Isolate invalid rows
    invalid_rows = df[mask]

    ## Write invalid rows to a file
    if not invalid_rows.empty:
        print(f"Writing invalid rows to {output_file}")
        invalid_rows[["id"]].to_csv(output_file, index=False)

    return df[~mask].drop("nuc_at_pos", axis=1)


def generate_is_mirna_column(df, grch):
    """
    Adds two columns to the input DataFrame:
    'is_mirna': 1 if the mutation falls within a miRNA region, 0 otherwise
    'mirna_accession': the accession number of the miRNA if 'is_mirna' is 1, None otherwise

    Args:
        df (pandas.DataFrame): The input DataFrame containing mutation data
        grch (int): The genome reference coordinate system version (e.g., 37, 38)

    Returns:
        pandas.DataFrame: The input DataFrame with two additional columns ('is_mirna' and 'mirna_accession')
    """
    # Construct the miRNA coordinates file path
    mirna_coords_file = os.path.join(MIRNA_COORDS_DIR, f"grch{grch}_coordinates.csv")

    # Load miRNA coordinates
    coords = pd.read_csv(mirna_coords_file)

    # Initialize new columns
    df['is_mirna'] = 0
    df['mirna_accession'] = None

    # Iterate over each mutation in the mutations dataframe
    for index, row in df.iterrows():
        mutation_chr = row['chr']
        mutation_start = row['pos']

        # Find matching miRNAs
        matching_rnas = coords.loc[(coords['chr'] == mutation_chr) &
                                   (coords['start'] <= mutation_start) &
                                   (coords['end'] >= mutation_start)]

        if not matching_rnas.empty:
            # Update the 'is_mirna' and 'mirna_accession' columns
            df.at[index, 'is_mirna'] = 1
            df.at[index, 'mirna_accession'] = matching_rnas['mirna_accession'].values[0]

    return df
    
def get_upstream_sequence(row, n=30):
    """
    Get the upstream sequence of length n from the given position.

    Args:
        row (pandas.Series): A row from the DataFrame containing the 'chr', 'pos', and 'ref_len' columns.
        n (int, optional): The length of the upstream sequence. Defaults to 30.

    Returns:
        str: The upstream sequence.
    """
    chrom = row['chr']
    pos = row['pos']
    upstream_start = max(1, pos - n)
    upstream_end = pos - 1
    return get_nucleotides_in_interval(chrom, upstream_start, upstream_end)

def get_downstream_sequence(row, n=30):
    """
    Get the downstream sequence of length n from the given position.

    Args:
        row (pandas.Series): A row from the DataFrame containing the 'chr', 'pos', and 'ref_len' columns.
        n (int, optional): The length of the downstream sequence. Defaults to 30.

    Returns:
        str: The downstream sequence.
    """
    chrom = row['chr']
    pos = row['pos']
    ref_len = len(row['ref'])
    downstream_start = pos + ref_len
    downstream_end = downstream_start + n - 1
    return get_nucleotides_in_interval(chrom, downstream_start, downstream_end)
    

def import_pyensembl(grch):
    if grch not in [37, 38]:
        raise ValueError("grch must be either 37 or 38")

    ens_release = 75 if grch == 37 else 109
    
    import os

    from pyensembl import EnsemblRelease
    os.environ['PYENSEMBL_CACHE_DIR'] = PYENSEMBL_CACHE_DIR
    assembly = EnsemblRelease(ens_release)
    assembly.download()
    assembly.index()
    
    return assembly


def generate_transcript_id_and_gene_name_columns(df, assembly):
    
    df['transcript_id'] = df.apply(lambda x: assembly.transcript_ids_at_locus(x['chr'], x['pos']), axis=1)
    # df["transcript_name"] = df.apply(lambda x: assembly.transcript_names_at_locus(x['chr'], x['pos']), axis=1)
    df["gene_name"] = df.apply(lambda x: assembly.gene_names_at_locus(x['chr'], x['pos']), axis=1)
    



def classify_and_save(df, file_name, save_path="results"):
    """
    Classifies mutations, saves to files, and logs the output.

    Parameters:
    - df (DataFrame): The DataFrame containing the mutation data.
    - file_name (str): The name of the file to be used as a prefix for the output files.
    - save_path (str): The directory where the output files will be saved. Default is "results".

    Returns:
    - None

    Side Effects:
    - Creates CSV files for case 1 and case 2 mutations (if any) in the specified save_path.
    - Logs information about the number of case 2 mutations found.
    """


    # Select case 1 and case 2 mutations
    case_1 = df[df.is_mirna == 0][["id", "wt_seq", "mut_seq"]]
    case_2 = df[df.is_mirna == 1][["id", "wt_seq", "mut_seq"]]

    # Create the directory if it does not exist
    os.makedirs(save_path, exist_ok=True)

    # Save case 1 mutations to a CSV file
    case_1_file = os.path.join(save_path, f"{file_name}_case_1.csv")
    case_1.to_csv(case_1_file, index=False)
    logging.info(f"Saved {len(case_1)} case 1 mutations to {case_1_file}")

    # Save case 2 mutations to a CSV file, if any
    if not case_2.empty:
        case_2_file = os.path.join(save_path, f"{file_name}_case_2.csv")
        case_2.to_csv(case_2_file, index=False)
        logging.info(f"Saved {len(case_2)} case 2 mutations to {case_2_file}")
    else:
        logging.info(f"No case 2 mutations were found for {file_name}")


def create_results_dataframe(wt_array, mutated_array):
    """Create a consolidated dataframe from the results."""
    wt_result_df = pd.DataFrame(wt_array)
    mut_result_df = pd.DataFrame(mutated_array)
    df = pd.concat([wt_result_df, mut_result_df])
    
    colnames = ["mrna_start", "mrna_end", "mrna_dot_bracket_5to3", "mirna_start", "mirna_end", "mirna_dot_bracket_5to3", "pred_energy", "mutation_id", "mirna_accession", "mrna_sequence", "mirna_sequence", "is_mutated"]
    df.columns = colnames
    
    df["id"] = df["mutation_id"].astype(str) + "_" + df["mirna_accession"].astype(str)
    df.drop(columns=["mutation_id"], inplace=True)
    return df

def generate_alignment_string_from_dot_bracket(df):
    """
    Generate an alignment string for each row in the DataFrame based on the miRNA dot-bracket structure.

    Args:
        df (pandas.DataFrame): A DataFrame containing columns 'mirna_start', 'mirna_dot_bracket_5to3', 'mirna_sequence', and 'mirna_end'.

    Returns:
        pandas.DataFrame: The input DataFrame with a new column 'alignment_string' added.
    """
    alignment_strings = []

    for _, row in df.iterrows():
        start_string = "0" * row.mirna_start
        mid_string = "".join("1" if char == ")" else "0" for char in row.mirna_dot_bracket_5to3)
        end_string = "0" * (len(row.mirna_sequence) - row.mirna_end - 1)

        alignment_string = start_string + mid_string + end_string
        alignment_strings.append(alignment_string)

    df["alignment_string"] = alignment_strings

    return df



def generate_match_count_columns(df):
    """
    Generate two new columns in the DataFrame: 'pred_num_basepairs' and 'pred_seed_basepairs'.

    Args:
        df (pandas.DataFrame): A DataFrame containing the 'alignment_string' column.

    Returns:
        pandas.DataFrame: The input DataFrame with two new columns added.
    """
    def count_ones(alignment_string, seed=False):
        """
        Count the number of '1' characters in the alignment string.

        Args:
            alignment_string (str): The alignment string.
            seed (bool, optional): If True, count only the '1' characters in the seed region (positions 2-7).

        Returns:
            int: The number of '1' characters in the alignment string or seed region.
        """
        if seed:
            return alignment_string[1:7].count("1")
        else:
            return alignment_string.count("1")

    df["pred_num_basepairs"] = df["alignment_string"].apply(count_ones)
    df["pred_seed_basepairs"] = df["alignment_string"].apply(lambda x: count_ones(x, seed=True))

    return df

def generate_ta_sps_columns(df):
    """
    Add 'ta_log10' and 'sps_mean' columns to the input DataFrame based on the miRNA seed sequence.

    Args:
        df (pandas.DataFrame): A DataFrame containing the 'mirna_sequence' column.

    Returns:
        pandas.DataFrame: The input DataFrame with 'ta_log10' and 'sps_mean' columns added.
    """
    # Generate temporary seed column
    df["seed"] = df["mirna_sequence"].str[1:8].str.replace("T", "U")
    # Read ta sps data
    ta_sps_df = pd.read_csv(TA_SPS_CSV, usecols=["seed_8mer", "ta_log10", "sps_mean"])
    ta_sps_df = ta_sps_df.rename(columns={"seed_8mer": "seed"})
    # Merge dataframes on seed column
    df = df.merge(ta_sps_df, on="seed", how="left")
    # Drop temporary column
    df.drop(columns=["seed"], inplace=True)

    return df


def generate_mre_sequence_for_vcf(vcf_df):
    """
    Generate the miRNA response element (MRE) sequence for each row in the input DataFrame.

    Args:
        vcf_df (pandas.DataFrame): A DataFrame containing columns 'mrna_sequence', 'mrna_end', and 'mirna_start'.

    Returns:
        pandas.DataFrame: The input DataFrame with a new column 'mre_region' added.
    """
    # Calculate miRNA length
    vcf_df["mirna_length"] = vcf_df["mirna_sequence"].str.len()

    # Calculate MRE coordinates
    vcf_df["mre_end"] = vcf_df["mrna_end"] + vcf_df["mirna_start"]
    vcf_df["mre_start"] = vcf_df["mre_end"] - vcf_df["mirna_length"]

    # Ensure MRE start is not negative
    vcf_df["mre_start"] = vcf_df["mre_start"].apply(lambda x: max(x, 0))

    # Extract MRE sequence
    vcf_df["mre_region"] = vcf_df.apply(lambda row: row["mrna_sequence"][row["mre_start"]:row["mre_end"]], axis=1)

    # Drop temporary column
    vcf_df = vcf_df.drop(columns=["mirna_length"])

    return vcf_df



def generate_important_sites(df):
    """
    Generate columns for important miRNA target site features.

    Args:
        df (pandas.DataFrame): A DataFrame containing the 'mre_region' and 'alignment_string' columns.

    Returns:
        pandas.DataFrame: The input DataFrame with additional columns for important site features.
    """
    # Anchor site
    df["anchor_a"] = (df["mre_region"].str[-1] == "A").astype(int)

    # Seed match features
    df["6mer_seed"] = (df["alignment_string"].str[1:7].str.count("0") == 0).astype(int)
    df["match_8"] = (df["alignment_string"].str[7] == "1").astype(int)
    df["6mer_seed_1_mismatch"] = (df["alignment_string"].str[1:7].str.count("0") == 1).astype(int)
    df["empty_seed"] = (df["alignment_string"].str[1:8].str.count("1") == 0).astype(int)

    # Compensatory and supplementary sites
    df["compensatory_site"] = (df["alignment_string"].str[12:17].str.count("0") == 0).astype(int)
    df["supplementary_site"] = (df["alignment_string"].str[12:16].str.count("0") == 0).astype(int)
    df["supplementary_site_2"] = (df["alignment_string"].str[16:21].str.count("0") == 0).astype(int)

    # Consecutive match
    df["9_consecutive_match_anywhere"] = (df["alignment_string"].str.contains("1{" + str(9) + ",}")).astype(int)

    return df

def generate_mirna_conservation_column(df):
    """
    Add a 'mirna_conservation' column to the input DataFrame based on miRNA conservation data.

    Args:
        df (pandas.DataFrame): A DataFrame containing the 'mirna_accession' column.

    Returns:
        pandas.DataFrame: The input DataFrame with a 'mirna_conservation' column added.
    """
    mirna_df = (pd.read_csv(MIRNA_CSV, 
                            usecols=["mirna_accession", "conservation"])
                            .rename(columns={"conservation": "mirna_conservation"})
                            [["mirna_accession", "mirna_conservation"]])

    df = df.merge(mirna_df, on="mirna_accession", how="left")
    return df


def generate_seed_type_columns(df):
    """
    Generate columns for different types of miRNA seed matches.

    Args:
        df (pandas.DataFrame): A DataFrame containing columns for various miRNA target site features.

    Returns:
        pandas.DataFrame: The input DataFrame with additional columns for seed match types.
    """
    # Canonical seed matches
    has_anchor_a = (df['anchor_a'] == 1)
    has_6mer_seed = (df['6mer_seed'] == 1)
    has_match_8 = (df['match_8'] == 1)
    no_supplementary_sites = (df['supplementary_site'] == 0) & (df['supplementary_site_2'] == 0)

    df['seed_8mer'] = (has_anchor_a & has_6mer_seed & has_match_8).astype(int)
    df['seed_7mer_a1'] = (has_anchor_a & has_6mer_seed & ~has_match_8).astype(int)
    df['seed_7mer_m8'] = (~has_anchor_a & has_6mer_seed & has_match_8 & no_supplementary_sites).astype(int)

    # Non-canonical seed matches
    has_compensatory_site = (df['compensatory_site'] == 1)
    has_6mer_seed_1_mismatch = (df['6mer_seed_1_mismatch'] == 1)
    has_supplementary_site = (df['supplementary_site'] == 1)
    has_supplementary_site_2 = (df['supplementary_site_2'] == 1)
    has_empty_seed = (df['empty_seed'] == 1)
    has_9_consecutive_match = (df['9_consecutive_match_anywhere'] == 1)
    has_many_basepairs = (df['pred_num_basepairs'] > 10)

    df['seed_compensatory'] = (has_compensatory_site & has_6mer_seed_1_mismatch & has_match_8).astype(int)
    df['seed_clash_2'] = (has_supplementary_site & has_6mer_seed & has_match_8).astype(int)
    df['seed_clash_3'] = (has_supplementary_site_2 & has_6mer_seed & has_match_8).astype(int)
    df['seed_clash_4'] = (has_empty_seed & has_9_consecutive_match).astype(int)
    df['seed_clash_5'] = (has_many_basepairs & ~has_6mer_seed).astype(int)

    return df

def calculate_au_content(sequence):
    au_count = sequence.count(
        'A') + sequence.count('T') + sequence.count('U')
    return None if len(sequence) == 0 else au_count / len(sequence)


def generate_mre_au_content_column(df):
    df["mre_au_content"] = df['mre_region'].apply(calculate_au_content)

    return df


def generate_au_content_column_for_vcf(vcf_df):

    vcf_df["local_au_content"] = vcf_df['mrna_sequence'].apply(calculate_au_content)
    
    return vcf_df


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

def make_predictions(df_with_features):
    """Make predictions using the XGBoost regressor."""
    df_filtered = filter_columns_for_xgb_prediction(df_with_features)
    return make_predictions_regressor(df_with_features, df_filtered)



def handle_output_dir(output_dir):
    """Create output directory if it doesn't exist."""
    os.makedirs(output_dir, exist_ok=True)



def prepare_jobs_from_df(df):
    # Load the MIRNA_CSV file
    mirna_df = pd.read_csv(MIRNA_CSV)

    # Extract relevant data from the input DataFrame
    wt_sequences = df["wt_seq"].tolist()
    mutated_sequences = df["mut_seq"].tolist()
    identifiers = df["id"].tolist()
    mirna_identifiers = mirna_df["mirna_accession"].tolist()
    mirna_sequences = mirna_df["sequence"].tolist()

    # Create pairs of wild-type sequences and miRNA sequences
    wt_pairs = [(wt, mirna) for wt in wt_sequences for mirna in mirna_sequences]
    mutated_pairs = [(mutated, mirna) for mutated in mutated_sequences for mirna in mirna_sequences]
    id_pairs = [(identifier, mirna_id) for identifier in identifiers for mirna_id in mirna_identifiers]

    # Create sets of job tuples
    wt_jobs_set = {wt_pair + id_pair for wt_pair, id_pair in zip(wt_pairs, id_pairs)}
    mutated_jobs_set = {mutated_pair + id_pair for mutated_pair, id_pair in zip(mutated_pairs, id_pairs)}

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


 
def run_jobs_multithreaded(job_list, binary_value):
    """
    Run jobs in a multithreaded manner using ProcessPoolExecutor.

    Args:
        job_list (list): A list of jobs to be executed. Each job should be a tuple or list of arguments to be passed to the `run_rnaduplex_multithreaded` function.
        binary_value (any): A binary value to be appended to the result of each job. 0 for wild-type, 1 for mutated.

    Returns:
        list: A list of tuples, where each tuple contains the result of a job and the binary value.
    """
    with ProcessPoolExecutor() as executor:
        future_jobs = [executor.submit(
            run_rnaduplex_multithreaded, *job) for job in job_list]

        results = []

        for future in as_completed(future_jobs):
            result = future.result()
            result_with_binary = result + (binary_value,)
            results.append(result_with_binary)

    return results
 
 
def save_results_to_disk(results, result_file):
    """
    Save the results to a file.

    Args:
        results (list): A list of tuples, where each tuple contains the result of a job and the binary value.
        result_file (str): The path to the file where the results should be saved.
    """
    with open(result_file, 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(results)
 

def get_memory_limit():
    mem = psutil.virtual_memory()
    total_memory = mem.total
    available_memory = mem.available
    memory_limit = min(total_memory, available_memory)
    return memory_limit


def print_memory_usage():
    # Get the current process
    process = psutil.Process()

    # Get the memory info for the current process
    mem_info = process.memory_info()

    # Get the total memory and used memory in bytes
    total_memory = psutil.virtual_memory().total
    used_memory = mem_info.rss

    # Calculate the memory usage percentage
    memory_percentage = (used_memory / total_memory) * 100

    # Print the memory usage information
    print(f"Used Memory: {used_memory / 1024 / 1024:.2f} MB")
    print(f"Memory Usage: {memory_percentage:.2f}%")

 
def main():
    # Get the memory limit
    memory_limit = get_memory_limit()
    print(f"Memory limit: {memory_limit / 1024 / 1024:.2f} MB")
    print_memory_usage()

    args = parse_cli_arguments()
    vcf_file_path = args.vcf
    start = args.start
    end = args.end
    output_dir = args.output_dir

    vcf_file_name = os.path.basename(vcf_file_path).split(".")[0]

    handle_output_dir(output_dir)
    rnaduplex_results_file = os.path.join(output_dir, f"{vcf_file_name}_{start}_{end}_rnaduplex.csv")

    df = load_vcf_into_df(vcf_file_path)
    print("Loaded VCF")
    print_memory_usage()

    if end == -1:
        end = len(df)

    df = df[start:end]
    print(f"Number of mutations: {len(df)}")
    print_memory_usage()

    augment_id(df)
    df = validate_ref_nucleotides(df)
    print("Augmented ID and validated reference nucleotides")
    print_memory_usage()

    df = generate_is_mirna_column(df, grch=37)
    print("Generated is_mirna column")
    print_memory_usage()

    df['upstream_seq'] = df.apply(lambda row: get_upstream_sequence(row, NUCLEOTIDE_OFFSET), axis=1)
    df['downstream_seq'] = df.apply(lambda row: get_downstream_sequence(row, NUCLEOTIDE_OFFSET), axis=1)
    df['wt_seq'] = df["upstream_seq"] + df["ref"] + df["downstream_seq"]
    df['mut_seq'] = df["upstream_seq"] + df["alt"] + df["downstream_seq"]
    print("Generated sequence columns")
    print_memory_usage()

    grch37 = import_pyensembl(37)
    print("Imported pyensembl")
    print_memory_usage()

    generate_transcript_id_and_gene_name_columns(df, grch37)
    print("Generated transcript ID and gene name columns")
    print_memory_usage()

    classify_and_save(df, vcf_file_name)
    print("Classified and saved mutations")
    print_memory_usage()

    df = pd.read_csv(os.path.join(output_dir, f"{vcf_file_name}_case_1.csv"))
    print("Loaded case 1 mutations")
    print_memory_usage()

    wt_jobs, mutated_jobs = prepare_jobs_from_df(df)
    del df
    gc.collect()
    print("Prepared jobs from dataframe")
    print_memory_usage()

    wt_result_array = run_jobs_multithreaded(wt_jobs, 0)
    print("Ran wild-type jobs")
    print_memory_usage()

    mut_result_array = run_jobs_multithreaded(mutated_jobs, 1)
    del wt_jobs, mutated_jobs
    gc.collect()
    print("Ran mutated jobs")
    print_memory_usage()

    save_results_to_disk(wt_result_array, rnaduplex_results_file)
    print("Saved wild-type results")
    print_memory_usage()

    save_results_to_disk(mut_result_array, rnaduplex_results_file)
    print("Saved mutated results")
    print_memory_usage()

    df = create_results_dataframe(wt_result_array, mut_result_array)
    print("Created results dataframe")
    del wt_result_array, mut_result_array
    gc.collect()
    print_memory_usage()

    df = generate_alignment_string_from_dot_bracket(df)
    df = generate_match_count_columns(df)
    df = generate_ta_sps_columns(df)
    df = generate_mre_sequence_for_vcf(df)
    df = generate_important_sites(df)
    df = generate_mirna_conservation_column(df)
    df = generate_seed_type_columns(df)
    df = generate_mre_au_content_column(df)
    df = generate_au_content_column_for_vcf(df)
    print("Generated feature columns")
    print_memory_usage()

    pred_df = make_predictions(df)
    del df
    gc.collect()
    print("Made predictions")
    print_memory_usage()

    meaningful_results_file = os.path.join(output_dir, f"{vcf_file_name}_{start}_{end}_meaningful_results.csv")

    pred_df[pred_df.pred_difference_binary != 0].to_csv(meaningful_results_file, index=False)
    del pred_df
    gc.collect()
    print("Saved meaningful results. Exiting. ###############################################")
    print_memory_usage()

if __name__ == "__main__":
    main()
