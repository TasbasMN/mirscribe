import argparse
import csv
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from socket import gethostname

import pandas as pd

from scripts.globals import *
from scripts.utils import *
from scripts.features import *

# parse cli
################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("-s", default=0, type=int, help="start index")
parser.add_argument("-e", default=1, type=int, help="end index")
parser.add_argument("--vcf", help="VCF file to analyze")
parser.add_argument("--output_dir", default="./results",
                    help="Output directory for runtime file, defaults to ./results")
args = parser.parse_args()
start = args.s
end = args.e
vcf_file = args.vcf
output_dir = args.output_dir
result_csv_file = f"results/sana_results_{start}_{end}.csv"


# importing dfs
################################################################################
df = pd.read_csv(vcf_file)
mirna_df = pd.read_csv(MIRNA_CSV)
df = df[start:end]


# prepare jobs
wt_jobs, mutated_jobs = prepare_jobs_from_df(df)
handle_target_file(result_csv_file)


# run jobs multithreaded
wt_result_array = run_jobs_multithreaded(wt_jobs, 0, result_csv_file)
mut_result_array = run_jobs_multithreaded(mutated_jobs, 1, result_csv_file)


# reading results
wt_result_df = pd.DataFrame(wt_result_array)
mut_result_df = pd.DataFrame(mut_result_array)
rnaduplex_results_df = pd.concat([wt_result_df, mut_result_df])


# preprocess for prediction
colnames = ["mrna_start", "mrna_end", "mrna_dot_bracket_5to3", "mirna_start", "mirna_end", "mirna_dot_bracket_5to3", "pred_energy", "mutation_id", "mirna_accession", "mrna_sequence", "mirna_sequence", "is_mutated" ]
rnaduplex_results_df.columns = colnames
rnaduplex_results_df["id"] = rnaduplex_results_df["mutation_id"].astype(str) + "_" + rnaduplex_results_df["mirna_accession"].astype(str)
rnaduplex_results_df.drop(columns=["mutation_id"], inplace=True)
df_chunks = split_df_to_num_thread_chunks(rnaduplex_results_df)


# running prediction pipeline
with ProcessPoolExecutor(max_workers=NUM_CORES) as executor:
    results = list(executor.map(apply_pipeline, df_chunks))
    
pipeline_results_df = pd.concat(results)
filtered_df = filter_columns_for_xgb_prediction(pipeline_results_df)
df = make_predictions_regressor(pipeline_results_df, filtered_df)

df.to_csv(f"results/sana_results_{start}_{end}_with_prediction.csv", index=False)

df2 = df[df.pred_difference_binary != 0]
df2.to_csv(f"results/sana_results_{start}_{end}_with_prediction_only_meaningful_results.csv", index=False)