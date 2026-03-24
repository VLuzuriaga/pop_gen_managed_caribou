This repository contains R pipelines for the analysis of caribou/reindeer genotype data generated using the _Rangifer_ SNP array. The workflows include data processing, quality control, molecular sexing, population structure analyses, and relatedness estimation. This project was developed by the Toronto Zoo as part of the Caribou Genome Biobanking Project (CCA-23-03-TZ), which aims to support future conservation breeding management of the ex-situ caribou population in Canada. 

Raw genotype data associated with this project have been deposited in the NCBI's Gene Expression Omnibus (GEO) database under accession number: [GSE325110](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325110)

### This is the analyses workflow:
1. [Libraries](https://github.com/VLuzuriaga/pop_gen_managed_caribou/blob/main/1_Libraries_to_load.R)
2. [Quality Control screening](https://github.com/VLuzuriaga/pop_gen_managed_caribou/blob/main/2_QC_analysis.R)
3. [Sex assignment](https://github.com/VLuzuriaga/pop_gen_managed_caribou/blob/main/3_Ho_SNP_in_x-chr.R)
4. [Inference of genetic structure](https://github.com/VLuzuriaga/pop_gen_managed_caribou/blob/main/4_Population_clustering.R)
5. [Genetic relatedness and kinship](https://github.com/VLuzuriaga/pop_gen_managed_caribou/blob/main/5_Genetic_relatedness.R)

<p align="center">
<img width="322" height="323" alt="image" src="https://github.com/user-attachments/assets/43f73f61-8bac-4b0b-81bb-a57a25701f50" />
</p>
For any questions regarding these analyses or the manuscript, please email me: vanessaluzuriagaaveiga(at)gmail.com
