#####################################################################################################################################################################################################################################################################

                                       # libraries

#####################################################################################################################################################################################################################################################################
import collections
from collections import defaultdict
import itertools
import re
import pickle
import sys

#####################################################################################################################################################################################################################################################################

                                       # sys_argv

#####################################################################################################################################################################################################################################################################
seqs_file = sys.argv[1]
start_sites_file = sys.argv[2]






# #####################################################################################################################################
# #
# #                           step 1 - retrieving start_positions data from pkl file
# #
# ######################################################################################################################################
# with open("start_positions_file.pkl", "rb") as start_positions_file:
#     start_positions = pickle.load(start_positions_file)



#####################################################################################################################################

#                             step 1 - retrieving start_positions data from txt file

#####################################################################################################################################
F= open(start_sites_file)
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








#####################################################################################################################################
#
#                          step 2 - extracting vienna format
#
######################################################################################################################################
F1 = open ("2_sequences_vienna_output.txt")
f1= F1.readlines()
# print(f)
fields1=[]
vienna_sequences = []

for line1 in f1:
    field1 = line1.strip("\n").split(" ")
    fields1.append(field1)

for i in range(1, len(f1), 2):
    vienna_sequences.append(fields1[i][0])

# for vienna_sequence in vienna_sequences:
#     print(vienna_sequence)






#####################################################################################################################################
#
#                          step 3 - pairing indexes
#
######################################################################################################################################
# F2 = open ("1_sequences.txt")
F2 = open (seqs_file)
f2= F2.readlines()
num_of_sequences = len(f2)
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
        # if sequences_2[i][j][5] != str(0):
        index_pairs[i].append((sequences_2[i][j][1], sequences_2[i][j][5]))
            # print(sequences_2[i][j][1])

# print(index_pairs)
# with open("3_index_pairs.pkl", "wb") as f2:
#     pickle.dump(index_pairs, f2)
#need to save this index_pairs list data and use it for next file












#####################################################################################################################################
#
#                          step 4 - constraints
#
######################################################################################################################################
RBS_window = 12
# --------below correct version---------this opens 25 bases if constrain is 12 bases------------------------------------
fields2 = []
for vienna_sequence in vienna_sequences:
    fields2.append(list(vienna_sequence.strip("\n")))

for i, j in enumerate(index_pairs):
    # print(type(i))
    for k, l in enumerate(j):
        # print(type(l[0]))
        # print(l)
        if int(l[1]) == 0:
            fields2[i][k] = "x"
        if int(l[1]) != 0 and int(l[0]) < int(l[1]) and (int(l[0]) < int(start_positions[i])-RBS_window)                                                    and     (int(l[1]) < int(start_positions[i])-RBS_window or int(l[1]) > int(start_positions[i])+RBS_window):
            pass
        if int(l[1]) != 0 and int(l[0]) < int(l[1]) and (int(l[0]) < int(start_positions[i])-RBS_window)                                                    and     (int(l[1]) >= int(start_positions[i])-RBS_window and int(l[1]) <= int(start_positions[i])+RBS_window):
            fields2[i][k] = "x"
        if int(l[1]) != 0 and int(l[0]) < int(l[1]) and (int(l[0]) >= int(start_positions[i])-RBS_window and int(l[0]) <= int(start_positions[i])+RBS_window)   and     (int(l[1]) >= int(start_positions[i])-RBS_window or int(l[1]) <= int(start_positions[i])+RBS_window): # no need of second condition
            fields2[i][k] = "x"
        if int(l[1]) != 0 and int(l[0]) < int(l[1]) and (int(l[0]) > int(start_positions[i])+RBS_window)                                                    and     (int(l[1]) > int(start_positions[i])+RBS_window):                                            # no need of second condition
            pass

        if int(l[1]) != 0 and int(l[0]) > int(l[1]) and (int(l[1]) < int(start_positions[i])-RBS_window)                                                    and     (int(l[0]) < int(start_positions[i])-RBS_window or int(l[0]) > int(start_positions[i])+RBS_window):
            pass
        if int(l[1]) != 0 and int(l[0]) > int(l[1]) and (int(l[1]) < int(start_positions[i])-RBS_window)                                                    and     (int(l[0]) >= int(start_positions[i])-RBS_window and int(l[0]) <= int(start_positions[i])+RBS_window):
            fields2[i][k] = "x"
        if int(l[1]) != 0 and int(l[0]) > int(l[1]) and (int(l[1]) >= int(start_positions[i])-RBS_window and int(l[1]) <= int(start_positions[i])+RBS_window)   and     (int(l[0]) >= int(start_positions[i])-RBS_window or int(l[0]) <= int(start_positions[i])+RBS_window): # no need of second condition
            fields2[i][k] = "x"
        if int(l[1]) != 0 and int(l[0]) > int(l[1]) and (int(l[1]) > int(start_positions[i])+RBS_window)                                                    and     (int(l[0]) > int(start_positions[i])+RBS_window):                                            # no need of second condition
            pass

# --------below correct version---------this opens 24 bases if constrain is 12 bases-------------
# for i, j in enumerate(index_pairs):
#     # print(type(i))
#     for k, l in enumerate(j):
#         # print(type(l[0]))
#         # print(l)
#         if int(l[1]) == 0:
#             fields2[i][k] = "x"
#         if int(l[1]) != 0 and int(l[0]) < int(l[1]) and (int(l[0]) < int(start_positions[i])-RBS_window)                                                    and     (int(l[1]) < int(start_positions[i])-RBS_window or int(l[1]) >= int(start_positions[i])+RBS_window):
#             pass
#         if int(l[1]) != 0 and int(l[0]) < int(l[1]) and (int(l[0]) < int(start_positions[i])-RBS_window)                                                    and     (int(l[1]) >= int(start_positions[i])-RBS_window and int(l[1]) < int(start_positions[i])+RBS_window):
#             fields2[i][k] = "x"
#         if int(l[1]) != 0 and int(l[0]) < int(l[1]) and (int(l[0]) >= int(start_positions[i])-RBS_window and int(l[0]) < int(start_positions[i])+RBS_window):   #and     (int(l[1]) >= int(start_positions[i])-RBS_window or int(l[1]) < int(start_positions[i])+RBS_window): # no need of second condition
#             fields2[i][k] = "x"
#         if int(l[1]) != 0 and int(l[0]) < int(l[1]) and (int(l[0]) >= int(start_positions[i])+RBS_window):                                                  #and     (int(l[1]) >= int(start_positions[i])+RBS_window):                                            # no need of second condition
#             pass
#
#         if int(l[1]) != 0 and int(l[0]) > int(l[1]) and (int(l[1]) < int(start_positions[i])-RBS_window)                                                    and     (int(l[0]) < int(start_positions[i])-RBS_window or int(l[0]) >= int(start_positions[i])+RBS_window):
#             pass
#         if int(l[1]) != 0 and int(l[0]) > int(l[1]) and (int(l[1]) < int(start_positions[i])-RBS_window)                                                    and     (int(l[0]) >= int(start_positions[i])-RBS_window and int(l[0]) < int(start_positions[i])+RBS_window):
#             fields2[i][k] = "x"
#         if int(l[1]) != 0 and int(l[0]) > int(l[1]) and (int(l[1]) >= int(start_positions[i])-RBS_window and int(l[1]) < int(start_positions[i])+RBS_window):    #and     (int(l[0]) >= int(start_positions[i])-RBS_window or int(l[0]) < int(start_positions[i])+RBS_window): # no need of second condition
#             fields2[i][k] = "x"
#         if int(l[1]) != 0 and int(l[0]) > int(l[1]) and (int(l[1]) >= int(start_positions[i])+RBS_window):                                                   #and     (int(l[0]) >= int(start_positions[i])+RBS_window):                                            # no need of second condition
#             pass
#----ends here
# #---------------------------------------------------------------------------------------
constraints = [[] for i in range(len(fields2))]

for i in range(len(fields2)):
    # print(len(fields2))
    constraints[i] = "".join(fields2[i]) # .join function will remove the internal list brackets

# print(constraints)
# print(len(constraints))
# for constraint in constraints:
#     print(constraint)
# print(constraints[20])

with open("constraints_file.pkl", "wb") as constraints_file:
    pickle.dump(constraints, constraints_file)











#####################################################################################################################################
#
#                          step 5 - sequences with constraints
#
######################################################################################################################################
seq_list = []

for sequence in f2:
    seq_list.append(sequence.strip("\n") + "\n")
# print(seq_list)
# print(len(seq_list))

for i in range(0, len(seq_list), 1):
    print(seq_list[i],constraints[i], sep="")










#####################################################################################################################################
#
#                          end
#
######################################################################################################################################
