Transcriptomics

Differential gene expression of mouse lung lymphocytes [Control vs OVA-Induced asthma mice]
Database: GEO
Link to FPKM file: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE212406&format=file&file=GSE212406%5Fgene%5Ffpkm%5Fmatrix%2Etxt%2Egz 

workflow:
GEO -> Asthma RNAseq data -> working environment as R -> extract only control vs test mice(72 h) -> Log transformation of FPKM -> QC -> PCA -> Identify upregulated and downregulated genes(bar plot) -> p value calc -> Volcano plot -> Heatmap -> GO analysis

Results:

On performing GO analysis, it comes to light that the highly upregulated genes in asthma mice code for proteins that are responsible for chemotaxis, leukocyte migration and regulation of leukocyte mediated immunity. 

The bar plot in the figures folder depicts the top 20 upregulated genes, based on the log2F values. For eg, Retnla has a log2F value of 9, which means that the gene is upregulated 10*9 times higher in the test samples when compared to the control.

From the given dataset, there are 774 genes whose |log2F| > 1 and -logp > 1.3, ie., these are the significant upregulated genes whose functions were identified using GO analysis.
