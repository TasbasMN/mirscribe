import gc
from scripts.main_helpers import *


def parse_cli_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="VCF file to analyze")
    parser.add_argument("-s", "--start", default=0, type=int, help="start index")
    parser.add_argument("-e", "--end", default=1, type=int, help="end index")
    parser.add_argument("-o", "--output_dir", default="./results", help="Output directory for runtime file, defaults to ./results")
    return parser.parse_args()

def main():
    
    args = parse_cli_arguments()
    input_file_name = args.vcf
    start = args.start
    end = args.end
    output_dir = args.output_dir
    
    handle_output_dir(output_dir)
    rnaduplex_results_file = os.path.join(output_dir, prepare_output_name(input_file_name, start, end))

    df = pd.read_csv(input_file_name)[start:end]
    wt_jobs, mutated_jobs = prepare_jobs_from_df(df)
    del df
    gc.collect()

    wt_result_array = run_jobs_multithreaded(wt_jobs, 0)
    mut_result_array = run_jobs_multithreaded(mutated_jobs, 1)
    del wt_jobs, mutated_jobs
    gc.collect()

    save_results_to_disk(wt_result_array, rnaduplex_results_file)
    save_results_to_disk(mut_result_array, rnaduplex_results_file)

    df = create_results_dataframe(wt_result_array, mut_result_array)
    del wt_result_array, mut_result_array
    gc.collect()

    df_with_features = add_features(df)
    del df
    gc.collect()

    df_with_predictions = make_predictions(df_with_features)
    del df_with_features
    gc.collect()

    meaningful_results_file = f"{output_dir}/{os.path.splitext(os.path.basename(rnaduplex_results_file))[0]}_meaningful_results.csv"
    df_with_predictions[df_with_predictions.pred_difference_binary != 0].to_csv(meaningful_results_file, index=False)
    del df_with_predictions
    gc.collect()

if __name__ == "__main__":
    main()
