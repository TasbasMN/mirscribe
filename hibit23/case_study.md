# Case Study for HIBIT23

This is the detailed document of case study presented on the poster.

## Data

The data is "PD3851a" sample taken from COSMIC database, which is a collection of breast cancer tumor mutations, published in Nik-Zainal et al. (2012).

https://cancer.sanger.ac.uk/cosmic/sample/overview?id=1230731

## Methodology

### Importing Data

The vcf dataset is processed by mirscribe pipeline, and results are imported into [processing pipeline](../scripts/truba_results_pipeline/).

Figure 1 below shows the top 10 entries of processed mutations. Each row represents the mutation's effect on a singular miRNA. The "id" column represents the mutation ID, which is encoded in the name_chr_position_ref_alt scheme. The "mirna_accession" column represents the accession of the miRNA that's affected. The "ENST" column contains the transcript ID of the mutation's target. The "gene_name" column contains the name of the aforementioned transcript ID. The "gain" column indicates whether the mutation enabled microRNA binding to that transcript, while the "loss" column indicates the opposite, i.e., whether the mutation resulted in the loss of microRNA binding to that transcript.

![](Pasted%20image%2020231003134616.png)
*Figure 1: Sample rows from processed analysis data*

### mRNA Trends

The primary objective is to group mutations based on their association with mRNA identifiers ('ENST').

	mRNA_trends = df.groupby('ENST_single').agg(
	    total_gains=pd.NamedAgg(column='gain', aggfunc='sum'),
	    total_losses=pd.NamedAgg(column='loss', aggfunc='sum'),
	    total_modifications=pd.NamedAgg(column='id', aggfunc='count')
	).reset_index().sort_values("total_modifications", ascending=False)
