April 25, 2016
Liheng Che
cheliheng@126.com

-------------------------------------
1.Filter_long_single_copy_seqs.py
-------------------------------------
This script used for identity sing-copy and long sequences using blastn. 

Usage: ./1.Filter_long_single_copy_seqs.py [-Q #] [-D #] [-e #] [-L #] [-C #] [-S #]
  -Q: # = Input file name(query file)
  -D: # = BLASTN database name(DataBase file)
  -e: # = Expectation value(E) threshold for saving hits
  -L: # = Expectation length threshold for saving single-copy sequences
  -C: # = Expectation coverage(the ratio of total length of locally aligned sequences over the length of query sequence) threshold of the second BLAST hit for saving query
  -S: # = Expectation similarity threshold of the second BLAST hit for saving query

The output files including:
1.blastn_results.xml (the blastn result of xml format);
2.single_copy_long_result.fas(the single-copy long sequences leaved);

Example: python 1.Filter_long_single_copy_seqs.py  -Q Tribolium_castaneum_exon_0414_bu.fas -D Tribolium_castaneum.Tcas3.21.dna.toplevel.fa -e 1E-10 -L 600 -C 0.30 -S 0.50

Note:this script calls blastn and makeblastdb,so the BLASTN program is needed.
Installers are available from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/.


-------------------------------------
2.Auto_filter_alignments.py
-------------------------------------
This script used for filtering alignments that suitable for primer design. 

Usage: ./2.Auto_filter_alignments.py [-f #] [-o #] [-e #] [-L #] [-C #] [-S #]
  -f: # = Input Directory Name that contains Multiple Sequence Alignments(MSA);
  -o: # = Output Directory Name that contains filtered alignments suitable for primer design;
  -L: # = Expectation length threshold of the conserved columns.
  -I: # = Expectation Identity threshold of the each conserved blocks
  -A: # = Expectation length threshold apart from two conserved regions.

The output files or directories including:
1.Filter_alignments_results.xls (Summary information of each alignments);
2.Output Directory contains filtered alignments suitable for primer design;

Example: python 2.Auto_filter_alignments.py  -f Input_MSA -o Output_Filter_MSA -L 7 -I 0.80 -A 100


-------------------------------------
3.Auto_exclusion_introns.py
-------------------------------------
This script used to exicise introns roughly in targeted sequences of 95 NPCLs using blastx.

Usage: ./3.Auto_exclusion_introns.py [-q #]
  -q: # = Input file name 
The output files or directories including:
1.blastx_results.xml (the blastx result of xml format);
2.results_exclusion_introns.fas(the sequences exclusion introns);

Example: python 3.Auto_exclusion_introns.py -q CAD_contigs.fasta
This script calls blastx, makeblastdb, score_cal_rev.py and Tribolium_ref_pro_DB_beetles_final.fasta.

Note:this script calls blastx and makeblastdb,so the BLASTX program is needed.
Installers are available from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/.