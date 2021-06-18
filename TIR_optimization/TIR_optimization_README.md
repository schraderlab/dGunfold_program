# TIR_optimization
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
1.	results file - this file is labeled as final.txt
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