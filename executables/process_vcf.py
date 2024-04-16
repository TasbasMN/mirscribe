import os
import argparse
import pandas as pd
from pyensembl import EnsemblRelease
from scripts.features import generate_is_mirna_column
from scripts.utils import get_nucleotide_at_position, get_nucleotides_in_interval
from scripts.globals import PYENSEMBL_CACHE_DIR
from tqdm import tqdm

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
    """
    if not file_path.endswith('.vcf'):
        raise ValueError("Please provide a valid .vcf file.")

    try:
        return pd.read_csv(
            file_path,
            sep="\t",
            header=None,
            names=["chr", "pos", "id", "ref", "alt"],
            comment='#',
        )
    except FileNotFoundError as e:
        raise FileNotFoundError(f"The file at {file_path} was not found.") from e
    except pd.errors.EmptyDataError as e:
        raise ValueError(
            f"The file at {file_path} is empty or not in the expected VCF format."
        ) from e


def augment_id(df):
    """Augment the 'id' column with additional information."""
    df['id'] += '_' + df['chr'] + '_' + df['pos'].astype(str) + '_' + df['ref'] + '_' + df['alt']


def compute_allele_lengths(df):
    """Compute lengths of reference and alternate alleles."""
    df["ref_len"] = df["ref"].apply(len)
    df["alt_len"] = df["alt"].apply(len)


def fetch_nucleotides(df):
    """Fetch nucleotides based on the length of the affected sequences."""
    df.loc[df['ref_len'] > 1, 'fetched_nucleotides'] = df.apply(
        lambda x: get_nucleotides_in_interval(x['chr'], x['pos'], x["pos"]+x["ref_len"]-1), axis=1
    )
    df.loc[df['ref_len'] == 1, 'fetched_nucleotides'] = df.apply(
        lambda x: get_nucleotide_at_position(x['chr'], x['pos']), axis=1
    )
    df.loc[df['ref_len'] == 0, 'fetched_nucleotides'] = ""


def compare_nucleotides(df):
    """Compare the fetched nucleotides with the reference nucleotides."""
    df["is_nucleotides_same"] = df["fetched_nucleotides"] == df["ref"]


def prepare_sequences(df):
    """Prepares reference indices, extended sequences, wild type, and mutated sequences."""
    df["ref_start"] = df["pos"]
    df["ref_end"] = df["pos"] + df["ref_len"]
    df['upstream_sequence'] = df.apply(lambda row: get_nucleotides_in_interval(row['chr'], row['ref_start'] - 30, row["ref_start"] - 1), axis=1)
    df['downstream_sequence'] = df.apply(lambda row: get_nucleotides_in_interval(row['chr'], row['ref_end'], row["ref_end"] + 29), axis=1)
    df['sequence'] = df["upstream_sequence"] + df["ref"] + df["downstream_sequence"]
    df['mutated_sequence'] = df["upstream_sequence"] + df["alt"] + df["downstream_sequence"]


def import_pyensembl(grch):
    """Import PyEnsembl and download required data."""
    if grch not in [37, 38]:
        raise ValueError("grch must be either 37 or 38")

    ens_release = 75 if grch == 37 else 109
    assembly = EnsemblRelease(ens_release)
    os.environ['PYENSEMBL_CACHE_DIR'] = PYENSEMBL_CACHE_DIR
    assembly.download()
    assembly.index()
    return assembly


def generate_transcript_id_and_gene_name_columns(df, assembly):
    """Generate transcript IDs and gene names."""
    df['transcript_id'] = df.apply(lambda x: assembly.transcript_ids_at_locus(x['chr'], x['pos']), axis=1)
    df["gene_name"] = df.apply(lambda x: assembly.gene_names_at_locus(x['chr'], x['pos']), axis=1)


def classify_and_save(df, file_name, save_path="results"):
    """Classify mutations, save to files, and log the output."""
    case_1 = df[df.is_mirna == 0][["id", "sequence", "mutated_sequence"]]
    case_2 = df[df.is_mirna == 1][["id", "sequence", "mutated_sequence"]]

    os.makedirs(save_path, exist_ok=True)

    case_1.to_csv(f"{save_path}/{file_name}_case_1.csv", index=False)
    if not case_2.empty:
        print(f"{len(case_2)} case 2 mutation(s) were found on {file_name}.")
        case_2.to_csv(f"{save_path}/{file_name}_case_2.csv", index=False)
    else:
        print(f"No case 2 mutations were found on {file_name}.")




def main():
    parser = argparse.ArgumentParser(description="Process VCF files in a directory.")
    parser.add_argument("directory_path", type=str, help="Path to the directory containing VCF files")
    parser.add_argument("--save_dir", type=str, default="vcf_results", help="Directory to save the output files")
    args = parser.parse_args()

    directory_path = args.directory_path
    save_path = args.save_dir

    os.makedirs(save_path, exist_ok=True)

    # Import and initialize the PyEnsembl assembly
    grch37 = import_pyensembl(37)

    # Get a list of all VCF files in the directory
    vcf_files = [f for f in os.listdir(directory_path) if f.endswith(".vcf")]

    # Use tqdm to create a progress bar
    for filename in tqdm(vcf_files, desc="Processing VCF files", unit="file"):
        file_path = os.path.join(directory_path, filename)
        df = load_vcf_into_df(file_path)
        augment_id(df)
        compute_allele_lengths(df)
        fetch_nucleotides(df)
        compare_nucleotides(df)
        prepare_sequences(df)
        df = generate_is_mirna_column(df, grch=37)
        generate_transcript_id_and_gene_name_columns(df, grch37)
        classify_and_save(df, filename.split(".")[0], save_path)


if __name__ == "__main__":
    main()