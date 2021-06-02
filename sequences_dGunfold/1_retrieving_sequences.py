#####################################################################################################################################

#                                   libraries

#####################################################################################################################################
import sys
import pickle








#####################################################################################################################################

#                                   sys_argv files

#####################################################################################################################################
seqs_file = sys.argv[1]
start_positions_file = sys.argv[2]
mrna_length = int(sys.argv[3])







#####################################################################################################################################

#                             step 1 - seqs file

#####################################################################################################################################
seq_list = []
F_1 = open(seqs_file)
f_1 = F_1.readlines()

for sequence in f_1:
    seq_list.append(sequence.strip("\n"))






#####################################################################################################################################

#                             step 2 - indexing AUG

#####################################################################################################################################
F= open(start_positions_file)
f = F.readlines()
fields = []
for line in f:
    # fields.append(line.strip("\n").split("\t")) #split_function_creates_new_list
    fields.append(line.strip("\n"))
# print(fields)
# print(len(fields))

start_positions =[]

for i in range(0, len(fields), 1):
    start_positions.append(fields[i])
# print(start_positions)
# print(len(start_positions))





#####################################################################################################################################
#
#                             step 3 - retrieving sequences and updating start positions
#
#####################################################################################################################################
retrieved_seqs = []
start_positions_2 = []
for i in range(0, len(seq_list)):
    # for j in start_positions:
    if len(seq_list[i]) <= mrna_length:
        retrieved_seqs.append(seq_list[i])
        start_positions_2.append(int(start_positions[i]))
    if len(seq_list[i]) > mrna_length and int(start_positions[i]) <= int(mrna_length/2)+1:
        retrieved_seqs.append(seq_list[i][:mrna_length])
        start_positions_2.append(int(start_positions[i]))
    if len(seq_list[i]) > mrna_length and int(start_positions[i]) > int(mrna_length/2)+1:
        if len(seq_list[i])-int(start_positions[i])+1 < int(mrna_length/2):
            retrieved_seqs.append(seq_list[i][len(seq_list[i])-(len(seq_list[i])-int(start_positions[i])+1)-(mrna_length-len(seq_list[i])-int(start_positions[i])+1):])
            start_positions_2.append(mrna_length-(len(seq_list[i])-int(start_positions[i])+1)+1)
        if len(seq_list[i])-int(start_positions[i])+1 >= int(mrna_length/2):
            retrieved_seqs.append(seq_list[i][int(start_positions[i])-(int(mrna_length/2)+1):int(start_positions[i])+(int(mrna_length/2)-1)])
            start_positions_2.append(int(mrna_length/2)+1)


#####################################################################################################################################
#
#                           step 4 - dumping start positions and retrieved_seqs
#
######################################################################################################################################
with open("start_positions_2_file.pkl", "wb") as start_positions_2_file:
    pickle.dump(start_positions_2, start_positions_2_file)
# this will open in 3_dGunfold script as start_positions_2 data


#####################################################################################################################################
#
#                           step 5 - printing retrieved_seqs
#
######################################################################################################################################
for i in range (len(retrieved_seqs)):
    print(retrieved_seqs[i])
























#####################################################################################################################################
#
#                           end
#
######################################################################################################################################
