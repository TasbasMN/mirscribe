import argparse
import gc
from concurrent.futures import ProcessPoolExecutor

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
file_name = vcf_file.split("/")[1].split("_")[0]
output_dir = args.output_dir
result_csv_file = f"results/{file_name}_{start}_{end}.csv"


# importing dfs
################################################################################
df = pd.read_csv(vcf_file)
mirna_df = pd.read_csv(MIRNA_CSV)
df = df[start:end]


# prepare jobs
wt_jobs, mutated_jobs = prepare_jobs_from_df(df)

# free up RAM
del df
gc.collect()
print("gc 1 done")

handle_target_file(result_csv_file)


# run jobs multithreaded
wt_result_array = run_jobs_multithreaded(wt_jobs, 0, result_csv_file)
mut_result_array = run_jobs_multithreaded(mutated_jobs, 1, result_csv_file)

print("results are written to the disk.")

# reading results
wt_result_df = pd.DataFrame(wt_result_array)
mut_result_df = pd.DataFrame(mut_result_array)
rnaduplex_results_df = pd.concat([wt_result_df, mut_result_df])

# free up RAM
del wt_result_df, mut_result_df, wt_result_array, mut_result_array
gc.collect()
print("gc 2 done")


# preprocess for prediction
colnames = ["mrna_start", "mrna_end", "mrna_dot_bracket_5to3", "mirna_start", "mirna_end", "mirna_dot_bracket_5to3", "pred_energy", "mutation_id", "mirna_accession", "mrna_sequence", "mirna_sequence", "is_mutated" ]
rnaduplex_results_df.columns = colnames
rnaduplex_results_df["id"] = rnaduplex_results_df["mutation_id"].astype(str) + "_" + rnaduplex_results_df["mirna_accession"].astype(str)
rnaduplex_results_df.drop(columns=["mutation_id"], inplace=True)
df_chunks = split_df_to_num_thread_chunks(rnaduplex_results_df)

# free up RAM
del rnaduplex_results_df
gc.collect()
print("gc 3 done")


# running prediction pipeline
with ProcessPoolExecutor(max_workers=NUM_CORES) as executor:
    results = list(executor.map(apply_pipeline, df_chunks))

pipeline_results_df = pd.concat(results)

# free up RAM
del results
gc.collect()
print("gc 4 done")


filtered_df = filter_columns_for_xgb_prediction(pipeline_results_df)
df = make_predictions_regressor(pipeline_results_df, filtered_df)

# non-results for debugging purposes
df[df.pred_difference_binary == 0].to_csv(f"results/{file_name}_{start}_{end}_other_results.csv", index=False)

# meaningful results
df[df.pred_difference_binary != 0].to_csv(f"results/{file_name}_{start}_{end}_meaningful_results.csv", index=False)
