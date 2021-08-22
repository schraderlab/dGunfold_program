#-------------------------------------------------------------------------------------------------------------------
#flags
#-------------------------------------------------------------------------------------------------------------------
#!/bin/bash
#!/usr/bin/env bash
# Flags

while getopts :s:t:c:r: aflag; do
  case $aflag in
    s) seqs_file=$OPTARG;;
    t) start_sites_file=$OPTARG;;
    c) temperature=$OPTARG;;
    r) mrna_length=$OPTARG;;
    ?) echo "Invalid option $OPTARG";;

  esac
done


#-------------------------------------------------------------------------------------------------------------------
#variables
#-------------------------------------------------------------------------------------------------------------------
var1=$seqs_file
var2=$start_sites_file
var3=${temperature:-37}
var4=${mrna_length:-50}

#--------------------------------------------------------------------------------------------------------------------
#pipeline
#--------------------------------------------------------------------------------------------------------------------
python3 sequences_dGunfold_1_retrieving_sequences.py $var1 $var2 $var4 >1_list_of_sequences.txt
# python3 1_retrieving_50bases.py $var1 $var2 >1_list_of_sequences.txt

RNAfold --temp $var3 <1_list_of_sequences.txt >2_sequences_vienna_output.txt
#
python3 sequences_dGunfold_2_dot_bracket_format.py 1_list_of_sequences.txt >2_dot_bracket_format.txt #this not be used if running fold command

wc -l <$var1 >1_number_of_sequences.txt
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

python3 sequences_dGunfold_3_constraints.py >3_sequences_with_constraint.txt

RNAfold --temp $var3 -C <3_sequences_with_constraint.txt >4_sequences_with_constraint_vienna_output.txt

python3 sequences_dGunfold_4_calculate_ddG.py >sequences_dGunfold_output1_final.txt

rm 1_number_of_sequences.txt
rm 1_list_of_sequences.txt
rm 2_ct.txt
rm 2_dot_bracket_format.txt
rm 2_sequences_vienna_output.txt
rm 3_sequences_with_constraint.txt
rm 4_sequences_with_constraint_vienna_output.txt
rm constraints_file.pkl
rm start_positions_2_file.pkl
rm rna.ps

#-------------------------------------------------------------------------------------------------------------------
#end
#-------------------------------------------------------------------------------------------------------------------
