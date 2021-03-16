#-------------------------------------------------------------------------------------------------------------------
#flags
#-------------------------------------------------------------------------------------------------------------------
#!/bin/bash
#!/usr/bin/env bash
# Flags

while getopts :g:t:c:o:r: aflag; do
  case $aflag in
    g) genbank_file=$OPTARG;;
    t) tss_file=$OPTARG;;
    c) temperature=$OPTARG;;
    o) operon_file=$OPTARG;;
    r) mrna_length=$OPTARG;;
    ?) echo "Invalid option $OPTARG";;

  esac
done

#-------------------------------------------------------------------------------------------------------------------
#variables
#-------------------------------------------------------------------------------------------------------------------
var1=$genbank_file
var2=${tss_file:-notpresent.txt}
var3=${temperature:-37}
var4=${operon_file:-notpresent.txt}
var5=${mrna_length:-50}

# #--------------------------------------------------------------------------------------------------------------------
# #pipeline
# #--------------------------------------------------------------------------------------------------------------------
if [ -e $var2 -a -e $var4 ]; then #both files
  python3 1_c_read_gb_file_with_both.py $var1 $var2 $var4 $var5

elif [ -e $var4 ]; then #only operon file
  python3 1_d_read_gb_file_with_operon.py $var1 $var4 $var5

elif [ -e $var2 ]; then #only tss file
  python3 1_a_read_gb_file.py $var1 $var2 $var5

else #neither file
  python3 1_b_read_gb_file.py $var1 $var5

fi

python3 1_e_generating_seqs_file.py >1_list_of_sequences.txt

python3 1_f_generating_start_sites_file.py >1_list_of_start_positions.txt

RNAfold --temp $var3 <1_list_of_sequences.txt >2_sequences_vienna_output.txt

python3 2_dot_bracket_format.py 1_list_of_sequences.txt >2_dot_bracket_format.txt

wc -l <1_list_of_sequences.txt >1_number_of_sequences.txt
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

python3 3_dGunfold.py 1_list_of_sequences.txt 1_list_of_start_positions.txt >3_sequences_with_constraint.txt

RNAfold --temp $var3 -C <3_sequences_with_constraint.txt >4_sequences_with_constraint_vienna_output.txt

python3 4_calculate_ddG.py >_output1_final.txt

python3 5_summary_table.py >_output2_summary_table.txt


#--------------------------------------------------------------------------------------------------------------------
rm 1_list_of_sequences.txt
rm 1_list_of_start_positions.txt
rm 1_number_of_sequences.txt
rm 2_ct.txt
rm 2_dot_bracket_format.txt
rm 2_sequences_vienna_output.txt
rm 3_sequences_with_constraint.txt
rm 4_sequences_with_constraint_vienna_output.txt
rm constraints_file.pkl
rm distances_file.pkl
rm distances_new_file.pkl
rm features_file.pkl
rm features_new_file.pkl
rm forms_file.pkl
rm forms_new_file.pkl
rm left_cordinates_file.pkl
rm left_cordinates_new_file.pkl
rm locus_tags_file.pkl
rm locus_tags_new_file.pkl
rm right_cordinates_file.pkl
rm right_cordinates_new_file.pkl
rm strands_file.pkl
rm strands_new_file.pkl
rm true_file.pkl
rm true_new_file.pkl
rm fifty_bases_file.pkl
rm start_codons_file.pkl
rm start_codons_new_file.pkl
rm start_positions_file.pkl
rm rna.ps

if [ -e operon_nos_file.pkl ]; then
  rm operon_nos_file.pkl
fi

if [ -e operon_left_cordinates_file.pkl ]; then
  rm operon_left_cordinates_file.pkl
fi

if [ -e operon_right_cordinates_file.pkl ]; then
  rm operon_right_cordinates_file.pkl
fi


#####################################################################################################################################################################################################################################################################

                                         # end

#####################################################################################################################################################################################################################################################################
