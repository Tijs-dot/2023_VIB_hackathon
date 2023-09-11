# brainstorming how we can leverage statistical rigor and biological relevance in single-cell differential gene expression


Theoretically, you would expect that ignoring the correlation between cells of the same individual will lead to more significant results, but those will be mostly false positives. 
Pseudobulking could resolve this problem, but this a crude method. You loose information by averaging the expression. 

We would run:
The default DE tool(s) in Seurat.
Pseudobulking
Milo (that paper that Valeriya shared).

This seemed a good starting point for further investigation. 


Options for data:

- Microglia from the eggen dataset (grn cases) minus ALS microglia
- single cell from sala frigerio 2019 (mouse MG) (ask Nico if we can't access)


Interesting tool: miloDE
https://www.biorxiv.org/content/10.1101/2023.03.08.531744v1.full

