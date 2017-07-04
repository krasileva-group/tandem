Given a list of gene identifiers and a gene annotation in .gff3 format, tandem.sh uses bedtools closest to find neighbouring, and potentially overlapping, genes for each gene on the list.

Requirements:
 - latest bedtools version

Usage:
bash tandem.sh -k=<check k nearest neighbours, default: 2> -o=<OUTPUT_PREFIX (default: tandem)> --min-overlap=<minimum number of overlapping bases, default: 0 - no overlap required> --max-distance=<maximum distance between neighbouring genes to be reported, default=100bp> <list_of_gene-ids> <annotation_in_gff3-format>

tandem.sh will generate 5 output files:
 - OUTPUT_PREFIX.k<value of k>.FF.tsv
   Neighbouring/overlapping genes in forward/forward orientation
 - OUTPUT_PREFIX.k<value of k>.RF.tsv
   Neighbouring/overlapping genes in reverse/forward orientation
 - OUTPUT_PREFIX.k<value of k>.FR.tsv
   Neighbouring/overlapping genes in forward/reverse orientation
 - OUTPUT_PREFIX.k<value of k>.RR.tsv
   Neighbouring/overlapping genes in reverse/reverse orientation
 - OUTPUT_PREFIX.k<value of k>.full.tsv
   Overlapping genes with one gene completely containing the other
