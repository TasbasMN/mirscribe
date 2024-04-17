import argparse
import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from scripts.globals import *
from scripts.utils import *
from scripts.features import *


def handle_output_dir(output_dir):
    """Create output directory if it doesn't exist."""
    os.makedirs(output_dir, exist_ok=True)


def prepare_output_name(input_file_name, start, end):
    """Prepare output file name."""
    file_name = os.path.basename(input_file_name).split("_")[0]
    return f"{file_name}_{start}_{end}.csv"


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

def add_features(rnaduplex_results_df):
    """Apply the feature creation pipeline in parallel."""
    df_chunks = split_df_to_num_thread_chunks(rnaduplex_results_df)
    with ProcessPoolExecutor(max_workers=NUM_CORES) as executor:
        results = list(executor.map(apply_pipeline, df_chunks))
    return pd.concat(results)

def make_predictions(df_with_features):
    """Make predictions using the XGBoost regressor."""
    df_filtered = filter_columns_for_xgb_prediction(df_with_features)
    return make_predictions_regressor(df_with_features, df_filtered)