from scripts.globals import *

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


def generate_transcript_id_and_gene_name_columns(df, grch):
    
    assembly = import_pyensembl(grch)
    
    df['transcript_id'] = df.apply(lambda x: assembly.transcript_ids_at_locus(x['chr'], x['pos']), axis=1)
    # df["transcript_name"] = df.apply(lambda x: assembly.transcript_names_at_locus(x['chr'], x['pos']), axis=1)
    df["gene_name"] = df.apply(lambda x: assembly.gene_names_at_locus(x['chr'], x['pos']), axis=1)
    
    return df