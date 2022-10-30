# UPC-BSG-TaxonomicAssignment

The goal of this project is writing an own implementation of a taxonomic assignment task. 
First, a real NCBI taxonomic reference (nodes.dmp) is transformed into a tree. Afterwords, a file is 
read containing reads of genomes (a set of taxes compatible to a new genomic sequence). An algorithm is then 
implemented to find the best taxonomic assignement for the given genomic sequence based on maximum F-measure. 

1) Extract folders and navigate to this folder in CMD.
2) Download nodes.dmp (Taxonomy dataset) from ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
3) Make sure 'nodes.dmp' (taxonomy dataset) and 'sample.inp' (genome reads) are in the same directory!
4) Execute: 'python full_pipeline.py'

--> To change input (location or name), run 'python full_pipeline.py -h' for instructions. 