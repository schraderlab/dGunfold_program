#####################################################################################################################################
#
#                                   libraries
#
#####################################################################################################################################
import sys
import pickle







#####################################################################################################################################
#
#                                   sys_argv files
#
#####################################################################################################################################
seqs_file = sys.argv[1]







#####################################################################################################################################
#
#                                   step 1 - reading DNA sequences
#
#####################################################################################################################################
seq_list = []
F_1 = open(seqs_file)
f_1 = F_1.readlines()

#print(f_1)
for sequence in f_1:
    seq_list.append(sequence.strip("\n"))





#####################################################################################################################################
#
#                                   step 2 - reading RNAfold output file and formatting it in fasta format
#
#####################################################################################################################################
F = open ("2_sequences_vienna_output.txt")
f= F.readlines()

fields=[]
vienna_sequences = []

for line in f:
    field = line.strip("\n").split(" ")
    fields.append(field)

for i in range(1, len(f), 2):
    vienna_sequences.append(fields[i][0])




#####################################################################################################################################
#
#                                   step 3 - printing in fasta format
#
#####################################################################################################################################
for i in range (len(seq_list)):
    print(">\n" + seq_list[i] + "\n"+ vienna_sequences[i])









#####################################################################################################################################
#
#                                   end
#
#####################################################################################################################################