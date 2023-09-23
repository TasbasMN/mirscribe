from socket import gethostname
from multiprocessing import cpu_count


def get_RNAduplex_location():
    hostname = gethostname()
    if hostname == "nazif-pc":
        return "/usr/bin/RNAduplex"
    else:
        return "/truba/home/mtasbas/miniconda3/envs/venv/bin/RNAduplex"
    

RNADUPLEX_LOCATION = get_RNAduplex_location()

PYENSEMBL_CACHE_DIR = "/run/media/nazif/2F946E411BA61D49/data"
XGB_PIPELINE_DIR = "scripts/xgb_pipeline/results"
SANA_DIR = "data/sana"
GRCH37_DIR = "data/fasta/grch37"

CLASH_CSV = "data/clash/clash_parsed.csv"
TA_SPS_CSV = "data/ta_sps/ta_sps.csv"
MIRNA_CSV = "data/mirna/mirna.csv"

NUM_CORES = cpu_count()
XGB_MODEL = "models/model_with_no_close_proximity.json"