

### Scripts description
1. Downloading.sh - downloading tsv-files for all phenotypes from UK Biobank 
2. Data_prep.sh - preprocessing downloaded data and merging files into one table
3. tsv_parse.py - merging rows in tsv-file (mearging SNP associated with several phenotypes)  
4. Data_Visualization.r - visualisation data (distributions of SNPs eith signifficant associations)
5. Data_plot_preparation.sh, Data_preparation_for_plot.py - data preprocessing for Manhattan plot  
6. Data_preparation_for_heatmap.py - data preprocessing for heatmap  
7. Heatmap_array_build.py -  building heatmap with phenotypes  
8. Heatmap_viz.r - heatmap visualizing and hierarchical clustering
9. Clustered_manhattan_prep.py - script for data preparation for Manhattan plot with associated clusters 
10. Manhattan_plot_viz.r - Manhattan plot visualizing  
11. MAF_preparation.py, MAF_for_all_SNP.py, Specifying_data_for_phen_corr.sh - data preprocessing for further MAF modelling
12. MAF_correlation.r - calculating MAF correlations and modelling MAF with Poissoin distribution
13. pval_manhattan.py - preparing Manhattan plot with Poissoin modelling p-values
14. creating_cluster_diagram.py, Cluster_annotation.py,SNP_diagram.py - creating treemap of phenotypes
15. Common_SNP_jaccard.py, Common_SNP_bed_and_VCF.py - exracting SNPs common for 2 or more clusters and creating BED-files and VCF-files for further analysis
16. eqtl_analysis.py - calculating the number of significant eqtls for data
17. CLUMPING.sh - performing PLINK clumping on phenotypes and clusters of phenotypes
18. LSEA - tool for performing hypergeometric test on independent loci
18. HLA_enrich.py - performing enrichment for HLA-locus and HIST1 gene cluster
19. enrichment_analysis.py - summarizing enrichment data
21. var_stats.py - analysis of protein altering variants enrichments for clusters
22. graph_preparation.py - creating network of shared loci between clusters
23. find_cliques.py - searcing for dense subgraphs from obtained network
24. grapg_and_clust_diags.r - visualising enrichments and network of shared loci between clusters
