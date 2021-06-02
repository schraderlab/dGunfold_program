* sequences_dGunfold program
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
2.	Install the software based on the OS following the steps on the web-page - https://www.tbi.univie.ac.at/RNA/documentation.html

RNAstructure
The following are the steps to download RNAstructure software:
1.	Go to the web-page https://rna.urmc.rochester.edu/RNAstructure.html and select Download RNAstructure to run software locally option.
2.	It will take you to the registration page. 
3.	After registration, it will take you to a web-page where you can download the Text (Command-Line) Interfaces which are provided as GZipped tarball. 
	Download either Linux 64-bit or 32-bit depending on your system. However, the sequences_dGunfold program has been written and tested with 32-bit version.
4.	Extract the files from the tarball as explained on https://rna.urmc.rochester.edu/Overview/Installation_Instructions_Linux.html
5.	To build the package, follow steps on https://rna.urmc.rochester.edu/Text/Building.html     
6.	After the package is built, if the downloading and installing site is not sequences_dGunfold folder,
	make sure the RNAstructure folder is moved/copied to the sequences_dGunfold folder. 
7. 	For any other documentation, kindly refer to https://rna.urmc.rochester.edu/RNAstructureHelp.html


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
1.	results file - this file will be labeled as _output1_final.txt
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







******************************************************************************************************
*																									 *
*										End											 				 *
*																									 *
******************************************************************************************************