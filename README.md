# Yeast Supplementary

## code

   - **ClusterThreshold.R**:  Network PCC cutoff determination
   - **GO_experiment.R**: GO enrichment analysis for each cluster
   - **SCnorm.R**: RNA-Seq raw counts normalization by SCnorm
   - **SummaryCluster.R**: random test for GO enrichment analysis and summary for clustering performance
   - **downloadMetadata.R**: RNA-Seq metadata download from SRAdb and GEOdb
   - **edgeR.R**: RNA-Seq raw counts normalization by edgeR
   - **kallisto.sh**: RNA-Seq raw reads qunatification by kallisto
   - **phylostratr.R**: determine gene phylostratum level and heatmap
   - **ribo-seq.sh**: Ribo-Seq analysis
   
 ## supplemental file in DataHub
 
   <https://datahub.io/lijing28101/yeast_supplementary>
 
 - **Supplementary Materials.pdf**:
Supplementary figures and tables.
 
 - **divergent_pairs_pair.csv**:
A list all ORFs have bidirectional promoter feature as Vakirlis et al., including genomic location and their corresponding SGD annotated genes with correlation.

 - **geneInfo_cluster112and317.xlsx**:
Gene and ORFs infomation in cluster 112 and cluster 317 for case study.

 - **GOsummary.xlsx**:
Significant GO terms in each cluster. 

 - **phylostratr_heatmap.pdf**:
Gene by gene BLAST result heatmap from phylostratr.

 - **RNA-Seq_rawcounts.txt**:
Raw counts for RNA-Seq analysis from 3,457 samples by kallisto.

 - **Ribo-Seq_rawcounts.csv**:
Raw counts for Ribo-Seq analysis from 302 samples by ribotricer.

 - **Ribo-Seq_metadata.xlsx**:
306 Ribo-Seq metadata download from SRA.

 - **Saccharomyces_cerevisiae.gff**:
Gene model for Saccharomyces cerevisiae, including unannotated ORFs.

 - **faa folder**:
Including all protein sequence for 7 Saccharomyces species.

 - **Saccharomyces_sequence_source.xlsx**:
Saccharomyces sequence source and infomation.

 - **species.xlsx**:
124 species used for phylostratr

 - **UTR_match.csv**:
All ORFs located within UTR and their corresponding annotated genes.

 - **yeast_mog.zip**:
mog file for virsualization. 

  Unzip this file and click mog1.7.7 to open the mog file.

  Click  "Open another project" under "Open an Existing Project" to select S.cerevisiae_RNA-seq_3457_27.mog

  For the first time user, it may ask you to choose Data file and metadata. Data file is JL_2019-10-21_cpminfo.txt, metadata is JL_2018-02-22_metadata3457.csv. 

  27 columns are non-graph0related infomation. And then click "OK" to open the mog project.

  Genes and ORFs infomation are shown in the main windows, mog manual download from here: https://github.com/urmi-21/MetaOmGraph/tree/master/manual

 
