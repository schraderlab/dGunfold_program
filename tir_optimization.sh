#-------------------------------------------------------------------------------------------------------------------
#flags
#-------------------------------------------------------------------------------------------------------------------
#!/bin/bash
#!/usr/bin/env bash
# Flags

while getopts :s:l:c:r: aflag; do
  case $aflag in
    s) sequence=$OPTARG;;
    l) length=$OPTARG;;
    c) temperature=$OPTARG;;
    r) mrna_length=$OPTARG;;
    ?) echo "Invalid option $OPTARG";;

  esac
done

#-------------------------------------------------------------------------------------------------------------------
#variables
#-------------------------------------------------------------------------------------------------------------------
var1=$sequence
var2=${length:-5}
var3=${temperature:-37}
var4=${mrna_length:-50}

#--------------------------------------------------------------------------------------------------------------------
#pipeline
#--------------------------------------------------------------------------------------------------------------------
python3 tir_optimization_1_a_TIR_mutations.py $var1 $var2 $var4

if [ $? == 1 ]; then

echo "ERROR: The TIR region may contain a stop codon"

else

python3 tir_optimization_1_b_generating_seqs_file.py >seqs_file.txt

python3 tir_optimization_1_c_generating_start_sites_file.py >start_sites_file.txt

RNAfold --temp $var3 <seqs_file.txt >2_sequences_vienna_output.txt

python3 tir_optimization_2_dot_bracket_format.py seqs_file.txt >2_dot_bracket_format.txt #this not be used if running fold command

wc -l <seqs_file.txt >1_number_of_sequences.txt
input="1_number_of_sequences.txt"
while IFS= read -r line
do
  echo $(expr $line \* 3)
  for i in $(seq 1 3 $(expr $line \* 3)); do sed -n "$i",$(($i+2))p '2_dot_bracket_format.txt' >seq_$i.txt; done
  export DATAPATH=RNAstructure/data_tables/
  for i in $(seq 1 3 $(expr $line \* 3)); do RNAstructure/exe/./dot2ct seq_$i.txt ct_$i; done
  for i in $(seq 1 3 $(expr $line \* 3)); do cat ct_$i >> 2_ct.txt; done
  for i in $(seq 1 3 $(expr $line \* 3)); do rm seq_$i.txt; done
  for i in $(seq 1 3 $(expr $line \* 3)); do rm ct_$i; done
done < "$input"

python3 tir_optimization_3_constraints.py start_sites_file.txt seqs_file.txt >3_sequences_with_constraint.txt

RNAfold --temp $var3 -C <3_sequences_with_constraint.txt >4_sequences_with_constraint_vienna_output.txt

python3 tir_optimization_4_calculate_ddG.py >tir_optimization_output1_final.txt

rm 1_number_of_sequences.txt
rm 2_ct.txt
rm 2_dot_bracket_format.txt
rm 2_sequences_vienna_output.txt
rm 3_sequences_with_constraint.txt
rm 4_sequences_with_constraint_vienna_output.txt
rm constraints_file.pkl
rm TIRs_file.pkl
rm start_positions_file.pkl
rm rna.ps
rm seqs_file.txt
rm start_sites_file.txt

fi

#-------------------------------------------------------------------------------------------------------------------
#end
#-------------------------------------------------------------------------------------------------------------------
