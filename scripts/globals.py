from socket import gethostname
from multiprocessing import cpu_count


def get_RNAduplex_location():
    hostname = gethostname()
    if hostname == "nazo":
        return "/usr/local/bin/RNAduplex"
    else:
        return "/truba/home/mtasbas/miniconda3/envs/venv/bin/RNAduplex"
    
def get_pyensembl_cache_location():
    hostname = gethostname()
    if hostname == "nazo":
        return "/home/nazif/thesis/data"
    else:
        return "/truba/home/mtasbas/data"


RNADUPLEX_LOCATION = get_RNAduplex_location()
PYENSEMBL_CACHE_DIR = get_pyensembl_cache_location()

XGB_PIPELINE_DIR = "scripts/1_xgboost_train/results"
SANA_DIR = "data/sana"
GRCH37_DIR = "data/fasta/grch37"

CLASH_CSV = "data/clash/clash_parsed.csv"
TA_SPS_CSV = "data/ta_sps/ta_sps.csv"
MIRNA_CSV = "data/mirna/mirna.csv"

NUM_CORES = cpu_count()
XGB_MODEL = "models/model_with_no_close_proximity.json"

QUANTILE_RANGE = 0.25

MIRNA_COORDS_DIR = "data/mirna_coordinates"

NUCLEOTIDE_OFFSET = 30