orfanage is a method for identifying the best matching ORF for each transcript in the GTF file based on evidence from reference annotaitons. 
The method is designed to identify cases of known ORFs fitting the query transcript both with and without modifications, introduced by additional exons, 
alternative start and end sites, etc.

orfanage
-c -i -o [-f -p -r ]
Arguments:
	c/--cds	Comma-separated list of GTF filenames with known transcripts and CDS annotations
	f/--filter	Select the best fitting CDS where possible
	i/--input	Input GTF with transcripts to which CDSs are to be ported
	o/--output	Output name
	p/--ppptrack	PhyloCSF track (currently in development)
	r/--reference	Reference fasta
