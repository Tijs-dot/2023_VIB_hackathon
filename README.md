# brainstorming how we can leverage statistical rigor and biological relevance in single-cell differential gene expression


Theoretically, you would expect that ignoring the correlation between cells of the same individual will lead to more significant results, but those will be mostly false positives. 
Pseudobulking could resolve this problem, but this a crude method. You loose information by averaging the expression. 

We would run:
The default DE tool(s) in Seurat.
Pseudobulking
Milo (that paper that Valeriya shared).

This seemed a good starting point for further investigation. 

Here there might be some public data: https://www.nature.com/articles/s41467-023-37126-3#data-availability

Interesting tool: miloDE
https://www.biorxiv.org/content/10.1101/2023.03.08.531744v1.full

