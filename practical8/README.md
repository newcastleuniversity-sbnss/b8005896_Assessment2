###Author

b8005896  // source code for MMB8052 Bioinformatics Lab Report
25/01/2023


###Process

This continuous script contains the entire workflow process of Practical 8. Specific sections include:

-Data installation
- Data preprocessing and quality control (tximport and Salmon and rlog)
- Principal Component analysis and Sample Distance Analysis
- DESEq2 differential analysis
- Volcano plot generation
- Annotation


###Data

Important information about the dataset : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116583/soft/GSE116583_family.soft.gz

Sample table (.csv) data : https://raw.githubusercontent.com/sjcockell/mmb8052/main/practicals/practical_08/data/sample_table.csv

Download count (quant.sf) data: https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/results/counts.z


###Software

FASTQ data were downloaded and imported (October 2022) with sra-tools (fastq-dump v2.8.0).

Salmon (v1.9.0) used for quantification.

Statistical analysis was performed in R studio (version 4.3.1) using R studio add-on packages used from bioconductor (BiocManager 1.30.20):
- tximport(v1.26.1)
- DESeq2(v.1.38.3)
- biomaRt(v2.54.0)
- pheatmap(v.1.0.12)
- tidyverse
- ggrepel(v0.9.2)
- RcolorBrewer (v1.1-3)
- ggplot(v3.4.2)



