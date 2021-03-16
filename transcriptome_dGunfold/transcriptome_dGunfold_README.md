* transcriptome_dGunfold program
******************************************************************************************************
*																									 *
*										Dependencies												 *
*																									 *
******************************************************************************************************
Python
Python can be installed using the steps found at https://www.python.org/downloads/ 
Make sure the version of python is 3.8 or higher - this is the current requirement of Vienna software to run. 
So kindly check the most updated version of python required to run Vienna RNAfold command.

Biopython
Biopython can be installed using the steps found at https://biopython.org/wiki/Download

ViennaRNA-software download
The following are the steps to download ViennaRNA software:
1.	Download the source code version using the  Compile from Source Code option from the following site - https://www.tbi.univie.ac.at/RNA/index.html
2.	Install the software based on the OS following the steps on the webpage - https://www.tbi.univie.ac.at/RNA/documentation.html

RNAstructure
Download the RNAstructure folder from the following webpage https://rna.urmc.rochester.edu/RNAstructure.html
Make sure to download in the same folder as the program file.





******************************************************************************************************
*																									 *
*										Input file(s)/parameter(s)									 *
*																									 *
******************************************************************************************************
1. 	genbank file

2. 	TSS file  -	tab delimited with the TSS positions in the first column and strand information in the third column,
				second column can have the peak reads or it can also be left blank
				Make sure there is are column headings in the first line and the data begins at first line itself
				Kindly refer the example file in the repository labeled as tss_example_file.txt
	
3.  operon file - tab delimited with Operon_no. in first column,
				  left_cordinates in the second column,
				  right_cordinates in the third column and
				  strands in the fourth column
				  Make sure there is are column headings in the first line and the data begins at first line itself
				  Kindly refer the example file in the repository labeled as operon_example_file.txt

4.	growth temperature of the organism

5.	mRNA length to be selected for dGmRNA calculation




******************************************************************************************************
*																									 *
*										Program execution											 *
*																									 *
******************************************************************************************************
step 1:
Make the script executable by running the following command:
chmod +x transcriptome_dGunfold_pipeline.sh

step 2:
Run transcriptome_dGunfold_pipeline shell script with following options on linux command prompt as shown in the following example:
./transcriptome_dgunfold_pipeline.sh -g na1000_sequence.gb -o operon_example_file.txt -t tss_example_file.txt -c 28 -r 50

the following option is mandatory:
-g genbank file

the following options are optional:
-o operon file
-t tss file 
-c is the temperature you want the RNAfold program to run. If not provided the default is 37 degreesC.
-r mRNA length to be selected for dGmRNA calculation. If not provided the default is 50 bases. This cannot be an odd number.


******************************************************************************************************
*																									 *
*										Result/output file(s)										 *
*																									 *
******************************************************************************************************
This program will output two files:
1.	results file - this file will be labeled as _output1_final.txt
	this file will have following columns:
	operon_no 				- operon number assigned to the gene
	operon_left_cordinate 	- leftmost coordinate of the operon
	operon_right_cordinate	- rightmost coordinate of the operon
	locus_tag 				- gene name
	feature 				- whether the gene is CDS, ncRNA, etc
	strand 					- +/-
	form 					- if multiple TSS assigned, which form it is
	gene_left_cordinate 	- leftmost coordinate of the gene
	gene_right_cordinate	- rightmost coordinate of the gene
	TSS_pos 				- coordinate of TSS
	5_UTR_length 			- length of 5' untranslated region (UTR)
	sequence 				- sequence retrieved for dGunfold calculation
	dot-bracket_structure	- dot-bracket structure of that sequence
	constrain 				- dot-bracket structure of that sequence such that the ribosome binding site(RBS) is made single stranded. 
							  All the single bases are represented as x
	dGmRNA 					- the minimum free energy (mfe) of the sequence retrieved
	dGinit 					- the minimum free energy (mfe) of the sequence in the constrained form
	dGunfold				- the difference between dGmRNA and dGinit

2.	summary file - this file will be labeled as _output2_summary_table.txt
	this file will have following information about the transcriptome:
	total genes: 	
	total CDSs: 	
	mRNAs with tss: 	
	number of multiforms: 	
	leaderless mRNAs: 	
	leadered mRNAs: 	
	mRNAs having leadered and leaderless forms: 	
	total ncRNAs: 	
	total tRNAs: 	
	total rRNAs: 	
	all mRNAs starting with ATG: 	
	all mRNAs starting with GTG: 	
	all mRNAs starting with CTG: 	
	all mRNAs starting with TTG: 	
	all mRNAs starting with ATA: 	
	all mRNAs starting with ATC: 	
	all mRNAs starting with ATT: 	
	leaderless mRNAs starting with ATG: 	
	leaderless mRNAs starting with GTG: 	
	leaderless mRNAs starting with CTG: 	
	leaderless mRNAs starting with TTG: 	
	leaderless mRNAs starting with ATA: 	
	leaderless mRNAs starting with ATC: 	
	leaderless mRNAs starting with ATT: 	
	leadered mRNAs starting with ATG: 	
	leadered mRNAs starting with GTG: 	
	leadered mRNAs starting with CTG: 	
	leadered mRNAs starting with TTG: 	
	leadered mRNAs starting with ATA: 	
	leadered mRNAs starting with ATC: 	
	leadered mRNAs starting with ATT: 	








******************************************************************************************************
*																									 *
*										Notes														 *
*																									 *
******************************************************************************************************
1. 	RNAfold program is run with all default parameters like:
	a) 	no enforce constraint 
	b) 	isolated base pairing is allowed
	for other options, kindly visit RNAfold Manpage https://www.tbi.univie.ac.at/RNA/RNAfold.1.html
	the options can be incorporated at line 50 and 69 of the program (shell script) for specific analysis by the user 
	
2. 	mRNA length input cannot be an odd number

3.	If both TSS and operon data are used, then the program will execute the following steps:
	a) 	the genes will be assigned to the operons
	b) 	the transcription start sites (TSS) will be assigned to the operons
	c) 	the mRNA sequences will be retrieved based on the mRNA length specified and the TSS position if its the first gene of the operon. For instance,
		if the mRNA length is kept default of 50 bases, there are following possibilities:
		(i)		if length of entire_sequence is < 50:
				then program will retrieve entire sequence
		(ii)	if length of entire_sequence is >= 50 and length of 5'UTR is < 25:
				then program will retrieve 50 bases from start position
		(iii) 	if length of entire_sequence is >= 50 and length of 5'UTR is >= 25 and length of coding region is < 25:
				say if the coding region is only 20 bases and 5' UTR is 70 bases (total 90 bases), then the program will retrieve from 
				41th position to last position
		(iv)	if length of 5'UTR is >= 25 and length of coding region >= 25:
				then program will retrieve 50 bases such that the start codon is at the center, i.e., 26th position
	d) 	if the gene is not the first gene of the operon (i.e., it is internal gene), then based on the mRNA length specified the sequence will be retrieved
		such that the start codon is at the center. For instance, if the default mRNA length of 50 bases is selected the start codon will be at 26th position
		in the sequence
	
4. If only TSS data is used, then the program will execute the following steps:
	a)	the transcription start sites (TSS) will be assigned to the genes directly
	b)	the mRNA sequences will be retrieved based on the mRNA length specified and the TSS position. Same as point 3c)

5. If only operon data is used, then the program will execute the following steps:
	a)	the genes will be assigned to the operons
	b)	however, since the TSS information is not available, it will retrieve the mRNA sequences based on the mRNA length specified for all 
		genes (first or internal) such that the start codon is at the center
	c) 	the only advantage of having operon info is that all the genes will be well annotated.
	
6. If none of the data is used, then the program will execute the following steps:
	a)	the mRNA sequences for all the genes will be retrieved based on the mRNA length specified such that the start codon is at the center
	













******************************************************************************************************
*																									 *
*										End											 				 *
*																									 *
******************************************************************************************************
