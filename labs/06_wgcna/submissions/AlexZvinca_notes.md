I used Spearman correlation with an absolute correlation threshold of 0.7 and a variance threshold of 0.5.

The original dataset contained 500 genes and 200 cells. After filtering, we have 151 genes. The network contained 34 genes and 51 connections, and Louvain detected 9 modules.

Classical clustering groups genes or samples based on overall similarity. Co-expression networks connect genes based on correlations between their expression patterns. Louvain then identifies groups of strongly connected genes. However, correlation does not prove direct biological interaction.