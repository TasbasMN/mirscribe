from socket import gethostname

def get_RNAduplex_location():
    hostname = gethostname()
    if hostname == "nazif-pc":
        return "/usr/bin/RNAduplex"
    else:
        return "/truba/home/mtasbas/miniconda3/envs/venv/bin/RNAduplex"
    

RNADUPLEX_LOCATION = get_RNAduplex_location()
PYENSEMBL_CACHE_DIR = "/run/media/nazif/2F946E411BA61D49/data"
XGB_PIPELINE_RESULTS = "scripts/xgb_pipeline/results"
CLASH_CSV = "data/clash/clash_parsed.csv"