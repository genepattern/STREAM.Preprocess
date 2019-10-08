# STREAM.Preprocess

Here we will first normalize the raw gene expression values based on library size. Then the gene expression values will be logarithmized. The mitochondrial genes will be removed.

We can filter out cells based on one of the metrics, including minimum number of genes expressed, minimum percentage of genes expresse, and minimum number of read count for one cell.

We can also filter out genes based on one of the metrics, including minimum number of cells expressing one gene, minimum percentage of cells expressing one gene, and minimum number of read count for one gene.
