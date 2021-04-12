# sup45-ts_readthrough
Analysis pipeline, codes, and processed data for Mangkalaphiban et al.<br/>
DOI: [pending]<br/>
Preprint: https://doi.org/10.1101/2020.12.15.422930<br/>

Raw sequencing data generated in this study are deposited and available at Gene Expression Omnibus (GEO) under accession number GSE162780.<br/>
Numerical data underlying the plots and the R codes used to generate them are in the folder _**Figures**_<br/>

## Analysis pipeline
1. Sequence alignment
	* Input:
		* Raw sequencing data used in this study and their associated SRR numbers are listed in *Raw_data_SRR.csv*
	* Transcriptome used for sequence alignment is available at https://github.com/Jacobson-Lab/yeast_transcriptome_v5
	* Output: 
		* Transcript abundance files generated by RSEM are in the folder _**RSEM_results**_
		* bam files
2. Calculate read P-site using riboWaltz
	* ***scripts**/read_p-site_riboWaltz.Rmd*
	* riboWaltz: https://github.com/LabTranslationalArchitectomics/riboWaltz
	* Input:
		* bam files from step 1 
		* annotation file (***RData**/genedf_riboWaltz_v5_CDS_corrected.txt*)
		* files containing p-site offset for each read length (***RData**/(sample)_psite_offset_adj.txt*)
	* Output:
		* ***RData/(dataset)**/(sample)_reads_psite_list.txt.gz* contains the following reads information: length of read, 5' & 3' ends + position of the read's P-site relative to the annotation provided in the annotation file. Distance from read's P-site to start and stop codons are calculated. The mRNA region (5'-UTR, CDS, or 3'-UTR) the P-site falls in is assigned.
3. Read count by mRNA region and readthrough efficiency calculation
	* ***scripts**/rt_efficiency.R*
	* Input: 
		* reads_psite_list from step 2
		* ***RData**/next_inframe_stop.txt*
	* Output:
		* ***RData**/rte_f0_cds_m3=33.Rdata*
4. Random forest analyses
	* ***scripts**/random_forest.Rmd*
	* Input:
		* mRNA features (***RData**/feature_file.csv*)
		* Readthrough efficiency calculated for each gene from step 3 (***RData**/rte_f0_cds_m3=33.Rdata*)
	* Output: 
		* Model accuracy
		* Feature importance
