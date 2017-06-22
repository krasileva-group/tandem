Given a list of gene identifiers and a gene annotation in .gff3 format, tandem.sh uses bedtools closest to find neighbouring, and potentially overlapping, genes for each gene on the list.

Requirements:
 - latest bedtools version

Usage:
bash tandem.sh -o=<OUTPUT_PREFIX (default: tandem)> --min-overlap=<minimum number of overlapping bases, default: 0 - no overlap required> --max-distance=<maximum distance between neighbouring genes to be reported> <list_of_gene-ids> <gff3>

tandem.sh will generate 5 output files:
 - OUTPUT_PREFIX.FF.tsv
   Neighbouring/overlapping genes in forward/forward orientation
 - OUTPUT_PREFIX.RF.tsv
   Neighbouring/overlapping genes in reverse/forward orientation
 - OUTPUT_PREFIX.FR.tsv
   Neighbouring/overlapping genes in forward/reverse orientation
 - OUTPUT_PREFIX.RR.tsv
   Neighbouring/overlapping genes in reverse/reverse orientation
 - OUTPUT_PREFIX.full.tsv
   Overlapping genes with one gene completely containing the other 
