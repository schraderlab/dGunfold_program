#####################################################################################################################################
#
#                                   libraries
#
#####################################################################################################################################
import collections
from collections import defaultdict
import itertools
import re
import pickle
import sys

#####################################################################################################################################
#
#                                   sys_argv files
#
#####################################################################################################################################








#####################################################################################################################################
#
#                           step 1 - retrieving start_positions data
#
######################################################################################################################################
with open("start_positions_2_file.pkl", "rb") as start_positions_2_file:
    start_positions_2 = pickle.load(start_positions_2_file)






#####################################################################################################################################
#
#                          step 2 - extracting sequences and vienna formats
#
######################################################################################################################################
F1 = open ("2_sequences_vienna_output.txt")
f1= F1.readlines()

fields1 = []
seq_list = []
vienna_sequences = []

for line1 in f1:
    field1 = line1.strip("\n").split(" ")
    fields1.append(field1)


for i in range(0, len(f1), 2):
    seq_list.append(fields1[i][0])


for i in range(1, len(f1), 2):
    vienna_sequences.append(fields1[i][0])





#####################################################################################################################################
#
#                          step 3 - pairing indexes
#
######################################################################################################################################
num_of_sequences = len(seq_list)
# print(num_of_sequences)

F3 = open ("2_ct.txt")
f3 = F3.readlines()
# print(f)
fields3=[]
indexes = []
sequences = [[] for i in range(num_of_sequences)] # number of list = no of sequences # can also input sequences file and use its length
sequences_unpacked = []
sequences_2 = [[] for i in range(num_of_sequences)]
index_pairs = [[] for i in range(num_of_sequences)]

for line3 in f3:
    field3 = line3.strip("\n").replace("       ", ";").replace("      ", ";").replace("     ", ";").replace("    ", ";").replace("   ", ";").replace("  ", ";").replace(" ", ";").split(";")
    fields3.append(field3)
    a = field3[0]

#--------------------------putting fields in sequences-----------------------------------------
for i in range(0, len(fields3), 1):
    if fields3[i][2] == "":
        indexes.append(i)

indexes.append(len(fields3)+1)
# print(indexes)

index = 0
while index < len(indexes)-1:
    sequences[index].append(fields3[indexes[index]:indexes[index+1]])
    index += 1

for i in range(len(sequences)):
    sequences_unpacked.append(sequences[i][0])
# print(sequences_unpacked)
# print(len(sequences_unpacked))
# print(sequences_unpacked[1])

#-------------------------------removing ENERGY line-----------------------------------
for i in range(len(sequences_unpacked)):
    for j in range(len(sequences_unpacked[i])):
        if sequences_unpacked[i][j][2] != "":
            sequences_2[i].append(sequences_unpacked[i][j])
# print(sequences_2)

# #-----------------------------pairing indexes-----------------------------------------
for i in range(len(sequences_2)):
    for j in range(len(sequences_2[i])):
        index_pairs[i].append((sequences_2[i][j][1], sequences_2[i][j][5]))














#####################################################################################################################################
#
#                          step 4 - constraints
#
######################################################################################################################################
RBS_window = 12
# --------below version---------this opens 25 bases if RBS_window is 12 bases------------------------------------
fields2 = []
for vienna_sequence in vienna_sequences:
    fields2.append(list(vienna_sequence.strip("\n")))

for i, j in enumerate(index_pairs):
    for k, l in enumerate(j):
        if int(l[1]) == 0:
            fields2[i][k] = "x"
        if int(l[1]) != 0 and int(l[0]) < int(l[1]) and (int(l[0]) < int(start_positions_2[i])-RBS_window)                                                    and     (int(l[1]) < int(start_positions_2[i])-RBS_window or int(l[1]) > int(start_positions_2[i])+RBS_window):
            pass
        if int(l[1]) != 0 and int(l[0]) < int(l[1]) and (int(l[0]) < int(start_positions_2[i])-RBS_window)                                                    and     (int(l[1]) >= int(start_positions_2[i])-RBS_window and int(l[1]) <= int(start_positions_2[i])+RBS_window):
            fields2[i][k] = "x"
        if int(l[1]) != 0 and int(l[0]) < int(l[1]) and (int(l[0]) >= int(start_positions_2[i])-RBS_window and int(l[0]) <= int(start_positions_2[i])+RBS_window)   and     (int(l[1]) >= int(start_positions_2[i])-RBS_window or int(l[1]) <= int(start_positions_2[i])+RBS_window): # no need of second condition
            fields2[i][k] = "x"
        if int(l[1]) != 0 and int(l[0]) < int(l[1]) and (int(l[0]) > int(start_positions_2[i])+RBS_window)                                                    and     (int(l[1]) > int(start_positions_2[i])+RBS_window):                                            # no need of second condition
            pass

        if int(l[1]) != 0 and int(l[0]) > int(l[1]) and (int(l[1]) < int(start_positions_2[i])-RBS_window)                                                    and     (int(l[0]) < int(start_positions_2[i])-RBS_window or int(l[0]) > int(start_positions_2[i])+RBS_window):
            pass
        if int(l[1]) != 0 and int(l[0]) > int(l[1]) and (int(l[1]) < int(start_positions_2[i])-RBS_window)                                                    and     (int(l[0]) >= int(start_positions_2[i])-RBS_window and int(l[0]) <= int(start_positions_2[i])+RBS_window):
            fields2[i][k] = "x"
        if int(l[1]) != 0 and int(l[0]) > int(l[1]) and (int(l[1]) >= int(start_positions_2[i])-RBS_window and int(l[1]) <= int(start_positions_2[i])+RBS_window)   and     (int(l[0]) >= int(start_positions_2[i])-RBS_window or int(l[0]) <= int(start_positions_2[i])+RBS_window): # no need of second condition
            fields2[i][k] = "x"
        if int(l[1]) != 0 and int(l[0]) > int(l[1]) and (int(l[1]) > int(start_positions_2[i])+RBS_window)                                                    and     (int(l[0]) > int(start_positions_2[i])+RBS_window):                                            # no need of second condition
            pass


# #---------------------------------------------------------------------------------------
constraints = [[] for i in range(len(fields2))]

for i in range(len(fields2)):
    # print(len(fields2))
    constraints[i] = "".join(fields2[i])



with open("constraints_file.pkl", "wb") as constraints_file:
    pickle.dump(constraints, constraints_file)





#####################################################################################################################################
#
#                          step 5 - sequences with constraints
#
######################################################################################################################################
for i in range(0, len(seq_list), 1):
    print(seq_list[i],constraints[i], sep="\n")






#####################################################################################################################################
#
#                          end
#
######################################################################################################################################
