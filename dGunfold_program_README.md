######################################################################################################
######################################################################################################
# 										dGunfold_program											 #
######################################################################################################
######################################################################################################


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
The following are the steps to download RNAstructure software:
1.	Go to the web-page https://rna.urmc.rochester.edu/RNAstructure.html and select Download RNAstructure to run software locally option.
2.	It will take you to the registration page. 
3.	After registration, it will take you to a web-page where you can download the Text (Command-Line) Interfaces which are provided as GZipped tarball. 
	Download either Linux 64-bit or 32-bit depending on your system. However, the transcriptome_dGunfold program has been written and tested with 32-bit version.
4.	Extract the files from the tarball as explained on https://rna.urmc.rochester.edu/Overview/Installation_Instructions_Linux.html. this will create a folder RNAstructure.
5.	To build the package, follow steps on https://rna.urmc.rochester.edu/Text/Building.html     
6.	After the package is built, if the downloading and installing site is not transcriptome_dGunfold folder,
	make sure the RNAstructure folder is moved/copied to the transcriptome_dGunfold folder. 
7. 	For any other documentation, kindly refer to https://rna.urmc.rochester.edu/RNAstructureHelp.html







######################################################################################################
######################################################################################################
# 										transcriptome_dGunfold_program								 #
######################################################################################################
######################################################################################################


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
./transcriptome_dGunfold_pipeline.sh -g genbank_example_file.gb -o operon_example_file.txt -t tss_example_file.txt -c 28 -r 50

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
1.	results file - this file will be labeled as transcriptome_dGunfold_output1_final.txt
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

2.	summary file - this file will be labeled as transcriptome_dGunfold_output2_summary_table.txt
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
	leadered mRNAs with 1 nt leader length: 	
	leadered mRNAs with 2 nts leader length: 	
	leadered mRNAs with 3 nts leader length: 	
	leadered mRNAs with 4 nts leader length: 	
	leadered mRNAs with 5 nts leader length: 	
	leadered mRNAs with 6 nts leader length: 	
	leadered mRNAs with 7 nts leader length: 	
	leadered mRNAs with 8 nts leader length: 	
	leadered mRNAs with 9 nts leader length: 	
	leadered mRNAs with 10 nts leader length: 	
	leadered mRNAs with 11 nts leader length: 	
	leadered mRNAs with 12 nts leader length: 	
	leadered mRNAs with 13 nts leader length: 	
	leadered mRNAs with 14 nts leader length: 	
	leadered mRNAs with 15 nts leader length: 	
	leadered mRNAs with 16 nts leader length: 	
	leadered mRNAs with 17 nts leader length: 	
	leadered mRNAs with 18 nts leader length: 	
	leadered mRNAs with 19 nts leader length: 	
	leadered mRNAs with 20 nts leader length: 	








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
	
	
4. If none of the data is used, then the program will execute the following steps:
	a)	the mRNA sequences for all the genes will be retrieved based on the mRNA length specified such that the start codon is at the center
	




	
######################################################################################################
######################################################################################################
# 										sequences_dGunfold_program							     	 #
######################################################################################################
######################################################################################################


******************************************************************************************************
*																									 *
*										Input file(s)/parameter(s)									 *
*																									 *
******************************************************************************************************
1.	sequences file  -	list of all the sequences one below other 
						Kindly refer the example file in the repository labeled as seqs_example_file.txt
						
	
2.	start site positions file -	list of all start codon positions with respect to the input sequences in the same order as that of sequences
								Kindly refer the example file in the repository labeled as starts_positions_example_file.txt
								
3. growth temperature of the organism

4.	mRNA length to be selected for dGmRNA calculation


******************************************************************************************************
*																									 *
*										Program execution											 *
*																									 *
******************************************************************************************************
step 1:
Make the script executable by running the following command:
chmod +x sequences_dGunfold_pipeline.sh

step 2:
Run sequences_dGunfold_pipeline shell script with following options on linux command prompt as shown in the following example:
./sequences_dGunfold_pipeline.sh -s seqs_example_file.txt -t start_positions_example_file.txt -c 28 -r 50

the following option is mandatory:
-s is the file with all the sequences you want the dGunfold to be calculated.
-t is the file with all the start positions for the sequences in the -f file.

the following options are optional:
-c is the temperature you want the RNAfold program to run. If not provided the default is 37 degreesC.
-r mRNA length to be selected for dGmRNA calculation. If not provided the default is 50 bases. This cannot be an odd number.

******************************************************************************************************
*																									 *
*										Result/output file(s)										 *
*																									 *
******************************************************************************************************
This program will output one file:
1.	results file - this file will be labeled as sequences_dGunfold_output1_final.txt
	this file will have following columns:
	sequence 				- sequence used for dGunfold calculation
	dot-bracket_structure	- dot-bracket structure of that sequence
	constrain 				- dot-bracket structure of that sequence such that the ribosome binding site(RBS) is made single stranded. 
							  All the single bases are represented as x
	dGmRNA 					- the minimum free energy (mfe) of the sequence retrieved
	dGinit 					- the minimum free energy (mfe) of the sequence in the constrained form
	dGunfold				- the difference between dGmRNA and dGinit









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

3.  the mRNA sequences will be retrieved from the input sequences based on the mRNA length specified and the start codon positions specified. For instance,
		if the mRNA length is kept default of 50 bases, there are following possibilities:
		(i)		if length of input_sequence is < 50:
				then program will use entire input sequence
		(ii)	if length of input_sequence is >= 50 and length of region upstream of start codon is < 25:
				then program will retrieve 50 bases from start position
		(iii) 	if length of input_sequence is >= 50 and length of region upstream of start codon is >= 25 and length of coding region is < 25:
				say if the coding region is only 20 bases and length of region upstream of start codon is 70 bases (total 90 bases), 
				then the program will retrieve from 41th position to last position
		(iv)	if length of length of region upstream of start codon is >= 25 and length of coding region >= 25:
				then program will retrieve 50 bases such that the start codon is at the center, i.e., 26th position









######################################################################################################
######################################################################################################
# 										TIR_optimization_program									 #
######################################################################################################
######################################################################################################

******************************************************************************************************
*																									 *
*										Input file(s)/parameter(s)									 *
*																									 *
******************************************************************************************************
1.	sequence  -	the leaderless mRNA sequence whose TIR needs to be optimized
				the sequence directly needs to be inputted at command line, kindly refer below in program execution section
								
2.	growth temperature of the organism

3.	Translation initiation region (TIR) length to be selected for creating synonymous mutations

4.	mRNA length to be selected for dGmRNA calculation

******************************************************************************************************
*																									 *
*										Program execution											 *
*																									 *
******************************************************************************************************
step 1:
Make the script executable by running the following command:
chmod +x sequences_dGunfold_pipeline.sh


step2 (optional - RNAstructure program should work after correct installation, but if not): 
Make the dot2ct program of RNAstructure executable by changing directory to RNAstructure/tests/scripts 
and then running the following command:
chmod +x dot2ct.sh


step 3:
Run tir_optimization shell script with following options on linux command prompt as shown in the following example:
./tir_optimization.sh -s  atgatcgaggagatcgtagtagatgatgtagagagagg -l 5 -c 28 -r 50

the following option is mandatory:
-s is the leaderless mRNA sequence whose TIR needs to be optimized
-l is the length of TIR (except for start codon, synonymous mutations will be made in all other codons in all combinations). If not provided, the default is 5.

the following options are optional:
-c is the temperature you want the RNAfold program to run. If not provided the default is 37 degreesC.
-r mRNA length to be selected for dGmRNA calculation. If not provided the default is 50 bases.


******************************************************************************************************
*																									 *
*										Result/output file(s)										 *
*																									 *
******************************************************************************************************
This program will output one file:
1.	results file - this file is labeled as tir_optimization_output1_final.txt
	this file has following columns:
	sequence 				- all the sequences with synonymous mutations in the TIR.
	dot-bracket_structure	- dot-bracket structure of that sequence
	constrain 				- dot-bracket structure of that sequence such that the ribosome binding site(RBS) is made single stranded. 
							  All the single bases are represented as x
	dGmRNA 					- the minimum free energy (mfe) of the sequence
	dGinit 					- the minimum free energy (mfe) of the sequence in the constrained form
	dGunfold				- the difference between dGmRNA and dGinit









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
	

2. TIR length cannot be greater than mRNA length.













******************************************************************************************************
*																									 *
*										End											 				 *
*																									 *
******************************************************************************************************
