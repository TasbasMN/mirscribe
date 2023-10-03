# Case Study for HIBIT23

This is the detailed document of case study presented on the poster.

## Data

The data is "PD3851a" sample taken from COSMIC database, which is a collection of breast cancer tumor mutations, published in Nik-Zainal et al. (2012).

https://cancer.sanger.ac.uk/cosmic/sample/overview?id=1230731

## Methodology

### Importing Data

The vcf dataset is processed by mirscribe pipeline, and results are imported into [processing pipeline](../scripts/truba_results_pipeline/).

Figure 1 below shows the top 10 entries of processed mutations. Each row represents the mutation's effect on a singular miRNA. The "id" column represents the mutation ID, which is encoded in the name_chr_position_ref_alt scheme. The "mirna_accession" column represents the accession of the miRNA that's affected. The "ENST" column contains the transcript ID of the mutation's target. The "gene_name" column contains the name of the aforementioned transcript ID. The "gain" column indicates whether the mutation enabled miRNA binding to that transcript, while the "loss" column indicates the opposite, i.e., whether the mutation resulted in the loss of miRNA binding to that transcript.

![](Pasted%20image%2020231003134616.png)

*Figure 1: Sample rows from processed analysis data*

### mRNA Trends

The primary objective is to group mutations based on their association with mRNA identifiers ('ENST') and count how many gains and losses that transcript had suffered.

	mRNA_trends = df.groupby('ENST').agg(
	    total_gains=pd.NamedAgg(column='gain', aggfunc='sum'),
	    total_losses=pd.NamedAgg(column='loss', aggfunc='sum'),
	    total_modifications=pd.NamedAgg(column='id', aggfunc='count')
	).reset_index().sort_values("total_modifications", ascending=False)

After aggregating mutation details, gene names and biotypes are added to the dataframe and the dataframe is filtered to only "protein_coding" biotypes, as the vcf has mutations that fall into intergenomic regions such as pseudogenes and antisense regions.

	gene_names = dict(zip(df['ENST'], df['gene_name']))
	mRNA_trends["gene_names"] = mRNA_trends["ENST"].apply(
	    lambda x: gene_names.get(x, None)
	)
	mRNA_trends["biotype"] = mRNA_trends["gene_names"].apply(lambda x: g37.genes_by_name(x)[0].biotype if x is not None else None)
	
	mRNA_trends[mRNA_trends.biotype == 'protein_coding'].head(10)

The top 10 affected transcript ids are shown in Figure 2 below.

![](../Pasted%20image%2020231003143814.png)

*Figure 2: Top 10 most affected protein coding transcript ids*



| gene_name | cancer_linked? | reference             |  Details   |
| --------- | -------------- | --------------------- | --- |
| RBFOX1    | Yes            | Shen et al. (2020)    |   Regulates BBB, downregulated in cancer  |
| PTPRN2    | Yes            | Sorokin et al. (2015) |   Overexpressed in cancer  |
| DMD       |   Yes             |      Wang et al. (2014)                 |  Tumor suppressor   |
| CTNNA2    |       Yes         |           Fanjul-Fernández et al. (2013)            | Tumor suppressor    |
| CDC27     |           Yes     |     Qiu et al. (2016)                  | Regulates cell division    |
| OPCML     |          Yes      |   McKie et al. (2012)                    |    Tumor suppressor |
| GPC5      |        Yes        |    Yuan et al. (2016)                   |   Tumor suppressor  |
| KHDRBS2   |        -        |                       |     |
| TSC22D1   |       Yes         |      Cho et al. (2017)                 |  Tumor suppressor   |
| NRXN3     |       Yes         |           Liu et al. (2021)            |  Regulated by oncogene miRNA miR-431    |