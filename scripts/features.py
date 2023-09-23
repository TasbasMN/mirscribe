import pandas as pd
import numpy as np
from scripts.globals import *

def calculate_au_content(sequence):
    au_count = sequence.count(
        'A') + sequence.count('T') + sequence.count('U')
    return None if len(sequence) == 0 else au_count / len(sequence)

def slice_column(row):
    return row['full_sequence_of_transcript'][row['temp_start']:row['temp_end']]

def generate_alignment_string_from_dot_bracket(df):
    full_strings = []
    for _, row in df.iterrows():
        start_string = (row.mirna_start) * "0"
        mid_string = row["mirna_dot_bracket_5to3"].replace(".", "0").replace(")", "1")
        end_string = (len(row.mirna_sequence) - row.mirna_end -1) * "0"
        
        full_string = start_string + mid_string + end_string
        full_strings.append(full_string)

    df["alignment_string"] = full_strings

    return df


def generate_match_count_columns(df):

    def count_ones(str, seed=False):
        return str[1:7].count("1") if seed else str.count("1")

    df["pred_num_basepairs"] = df["alignment_string"].apply(count_ones)
    df["pred_seed_basepairs"] = df["alignment_string"].apply(
        count_ones, seed=True)

    return df

def generate_ta_sps_columns(df):
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

def generate_important_sites(df):
    df["anchor_a"] = (df["mre_region"].str[-1] == "A").astype(int)
    df["6mer_seed"] = (
        df["alignment_string"].str[1:7].str.count("0") == 0).astype(int)
    df["match_8"] = (df["alignment_string"].str[7] == "1").astype(int)
    df["6mer_seed_1_mismatch"] = (
        df["alignment_string"].str[1:7].str.count("0") == 1).astype(int)

    df["compensatory_site"] = (
        df["alignment_string"].str[12:17].str.count("0") == 0).astype(int)

    df["supplementary_site"] = (
        df["alignment_string"].str[12:16].str.count("0") == 0).astype(int)
    df["supplementary_site_2"] = (
        df["alignment_string"].str[16:21].str.count("0") == 0).astype(int)
    df["empty_seed"] = (
        df["alignment_string"].str[1:8].str.count("1") == 0).astype(int)

    df["9_consecutive_match_anywhere"] = (df["alignment_string"]
                                          .str
                                          .contains("1{" + str(9) + ",}")
                                          .astype(int))

    return df

def generate_mirna_conservation_column(df):
    mirna_df = (pd.read_csv(MIRNA_CSV, 
                            usecols=["mirna_accession", "conservation"])
                            .rename(columns={"conservation": "mirna_conservation"})
                            [["mirna_accession", "mirna_conservation"]])

    df = df.merge(mirna_df, on="mirna_accession", how="left")
    return df


def generate_seed_type_columns(df):
    df['seed_8mer'] = ((df['anchor_a'] == 1) & (
        df['6mer_seed'] == 1) & (df['match_8'] == 1)).astype(int)
    df['seed_7mer_a1'] = ((df['anchor_a'] == 1) & (
        df['6mer_seed'] == 1) & (df['match_8'] == 0)).astype(int)
    df['seed_7mer_m8'] = ((df['anchor_a'] == 0) & (df['6mer_seed'] == 1) & (df['match_8'] == 1) & (
        df['supplementary_site'] == 0) & (df['supplementary_site_2'] == 0)).astype(int)
    df['seed_compensatory'] = ((df['compensatory_site'] == 1) & (
        df['6mer_seed_1_mismatch'] == 1) & (df['match_8'] == 1)).astype(int)

    df['seed_clash_2'] = ((df['supplementary_site'] == 1) & (
        df['6mer_seed'] == 1) & (df['match_8'] == 1)).astype(int)
    df['seed_clash_3'] = ((df['supplementary_site_2'] == 1) & (
        df['6mer_seed'] == 1) & (df['match_8'] == 1)).astype(int)
    df['seed_clash_4'] = ((df['empty_seed'] == 1) & (
        df['9_consecutive_match_anywhere'] == 1)).astype(int)
    df['seed_clash_5'] = ((df['pred_num_basepairs'] > 10)
                          & (df['6mer_seed'] == 0)).astype(int)

    return df

def generate_mre_au_content_column(df):
    df["mre_au_content"] = df['mre_region'].apply(calculate_au_content)

    return df

def generate_local_au_content_column(df):

    # clip sets negative start indices to 0
    df["temp_start"] = (df.mre_start - 30).clip(0)

    # this handles mre_ends extending off the mRNA transcript
    df["temp_extended_end"] = df.mre_end + 30
    df["temp_transcript_length"] = df.full_sequence_of_transcript.str.len()
    df["temp_end"] = df[["temp_transcript_length",
                         "temp_extended_end"]].min(axis=1)

    # temp col for calculating au content
    df["temp_col_for_calculating_au_content"] = df.apply(slice_column, axis=1)

    df["local_au_content"] = df['temp_col_for_calculating_au_content'].apply(
        calculate_au_content)

    df.drop(["temp_start", "temp_extended_end", "temp_transcript_length",
            "temp_end", "temp_col_for_calculating_au_content"], axis=1, inplace=True)

    return df

def find_mres_in_close_proximity(df):
    # Calculate the middle of MRE
    df["middle_of_mre"] = ((df["mre_start"] + df["mre_end"]) // 2).astype(int)

    # Find the indices where the length of each group is greater than or equal to 2
    group_indices = df.groupby(["enst", "mirna_accession"]).indices
    filtered_indices = [indices for indices in group_indices.values() if len(indices) >= 2]

    # Create a dictionary to store the results
    result = {}

    # Loop through the filtered indices and perform the computation
    for indices in filtered_indices:
        # Get the values corresponding to the indices
        values = df.iloc[indices]["middle_of_mre"].values

        # Create a mask for the condition
        mask = np.logical_and(
            np.abs(values[:, None] - values) > (13 + 4),
            np.abs(values[:, None] - values) < (35 + 4)
        )

        # Get the indices where the condition is satisfied
        i, j = np.where(mask)

        # Add the pairs to the result dictionary
        for k in range(len(i)):
            # Create a key using the 'enst' and 'mirna_accession' values from the DataFrame
            key = (df.iloc[indices[i[k]]]['enst'], df.iloc[indices[i[k]]]['mirna_accession'])

            # If the key is not in the result dictionary, add it with an empty list as the value
            if key not in result:
                result[key] = []

            # Get the pair values
            a, b = values[i[k]], values[j[k]]

            # Check if the pair (a, b) or (b, a) is already in the list associated with the key
            if (a, b) not in result[key] and (b, a) not in result[key]:
                # Append the pair of values to the list associated with the key
                result[key].append((a, b))

    # Create a set of tuples for faster membership checking
    data_set = {(enst, mirna, i) for (enst, mirna), values in result.items() for i in values[0]}

    # Use vectorized operations to check membership and assign the result
    df['another_mre_in_close_proximity'] = np.where(
        df[['enst', 'mirna_accession', 'middle_of_mre']].apply(tuple, axis=1).isin(data_set),
        1,
        0
    )

    return df

def generate_is_mirna_column(df, grch):
    coords = pd.read_csv(f"data/mirna_coordinates/grch{grch}_coordinates.csv")
    df['is_mirna'] = 0
    df["mirna_accession"] = None
    
    # Iterate over each mutation in the mutations dataframe
    for index, row in df.iterrows():
        mutation_chr = row['chr']
        mutation_start = row['pos']

        # Check if the mutation falls into any of the miRNAs
        matching_rnas = coords[(coords['chr'] == mutation_chr) & (coords['start'] <= mutation_start) & (coords['end'] >= mutation_start)]

        if not matching_rnas.empty:
            # Update the 'is_mirna' column to 1 for the current mutation
            df.at[index, 'is_mirna'] = 1
            df.at[index, 'mirna_accession'] = matching_rnas['mirna_accession'].values[0]
    return df

#####
# vcf functions

def generate_positions_from_id(vcf_df):
    vcf_df['chr'] = vcf_df['id'].str.split('_').str[0]

    vcf_df['start_coordinate'] = vcf_df['id'].str.split('_').str[1].astype(int) - 30 + vcf_df["mrna_start"]
    vcf_df['end_coordinate'] = vcf_df['id'].str.split('_').str[1].astype(int) - 30 + vcf_df["mrna_end"]
    
    return vcf_df

def generate_mre_sequence_for_vcf(vcf_df):

    def slice_column(row):
        return row["mrna_sequence"][row["mre_start"]:row["mre_end"]]
    
    # getting mirna length
    vcf_df["mirna_length"] = vcf_df["mirna_sequence"].str.len()

    # using mirna length to figure out mre coordinates
    vcf_df["mre_end"] = vcf_df["mrna_end"] + vcf_df["mirna_start"]
    vcf_df["mre_start"] = vcf_df["mre_end"] - vcf_df["mirna_length"]

    # some start values might be lower than zero, so we need to adjust
    vcf_df["mre_start"] = vcf_df["mre_start"].apply(lambda x: max(x, 0))

    # creating mre sequence column
    vcf_df["mre_region"] = vcf_df.apply(slice_column, axis=1)

    # dropping temp column
    vcf_df.drop(columns=["mirna_length"], inplace=True)
    
    return vcf_df

def generate_au_content_column_for_vcf(vcf_df):

    vcf_df["local_au_content"] = vcf_df['mrna_sequence'].apply(calculate_au_content)
    
    return vcf_df

def apply_pipeline(df):
    df = generate_positions_from_id(df)
    df = generate_alignment_string_from_dot_bracket(df)
    df = generate_match_count_columns(df)
    df = generate_ta_sps_columns(df)
    df = generate_mre_sequence_for_vcf(df)
    df = generate_important_sites(df)
    df = generate_mirna_conservation_column(df)
    df = generate_seed_type_columns(df)
    df = generate_mre_au_content_column(df)
    df = generate_au_content_column_for_vcf(df)
    return df
    # df_filtered = filter_columns_for_xgb_prediction(df)
    # return make_predictions_regressor(df, df_filtered)


