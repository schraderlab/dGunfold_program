#####################################################################################################################################################################################################################################################################

                                       #libraries

#####################################################################################################################################################################################################################################################################
from itertools import repeat
import sys
import pickle
from operator import itemgetter, attrgetter
from Bio.Data import CodonTable
from Bio.Seq import Seq
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]







#####################################################################################################################################################################################################################################################################

                                       #sys_argv

#####################################################################################################################################################################################################################################################################
genbank = sys.argv[1]
tss = sys.argv[2]
operon = sys.argv[3]
mrna_length = int(sys.argv[4])










#####################################################################################################################################################################################################################################################################

                                       #file 1-extracting features from genbank file

#####################################################################################################################################################################################################################################################################
F = open(genbank)
f = F.readlines()
fields=[]
# print(f[:100])

for line in f:
    field = line.strip('\n')
    fields.append(field)
# print(fields[80121:80125])
# print(fields[85][:9])

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                       # genome sequence
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#----------------------------ORIGIN index-----------------------------------------------------------
for i in range(len(fields)):
    if fields[i][:6] == "ORIGIN":
        origin_index = i
# print(origin_index)

#----------------------------sequence-----------------------------------------------------------
seq = []
for i in range(origin_index, len(fields)-2, 1):
    seq.append(fields[i][10:].replace(" ", ""))

seq_2 = "".join(seq)
# print(seq_2)
# print(seq_2[3])
# print(len(seq_2))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                       # features
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------gene_count------------------------------------------------------------------------
gene_count = 0
for i in range(len(fields)):
    if fields[i][:9] == "     gene":
        gene_count += 1
# print(gene_count)

#-----------------------------gene_indexes------------------------------------------------------------------------
gene_indexes = []
for i in range(len(fields)):
    if fields[i][:9] == "     gene":
        gene_indexes.append(i)
gene_indexes.append(origin_index)
# print(gene_indexes)
# print(len(gene_indexes))

#-----------------------------all_features------------------------------------------------------------------------
all_features = []

for i, j in enumerate(gene_indexes):
    if i != len(gene_indexes)-1:
        all_features.append(fields[gene_indexes[i]:gene_indexes[i+1]])
        # coordinates.append()
# print(all_features[0])

for i, j in enumerate(all_features):
    for k, l in enumerate(all_features[i]):
        all_features[i][k] = l.replace(' ', '').replace('"', '')
# print(len(all_features))
# print(all_features)
# print(all_features[0])
# print(all_features[2299])
# print(all_features[0][4][-5:-1])

#-----------------------------coordinates and other info------------------------------------------------------------------------
left_cordinates = [[] for i in range(gene_count)]
right_cordinates = [[] for i in range(gene_count)]
features = [[] for i in range(gene_count)]
locus_tags = [[] for i in range(gene_count)]
strands = [[] for i in range(gene_count)]


dot_indexes = []
for i, j in enumerate(all_features):
    dot_indexes.append(j[0].index('..'))
# print(dot_indexes[0])

for i,j in enumerate(dot_indexes):
# for i, j in enumerate(all_features):
    if all_features[i][0][4:14] == "complement":
        strands[i].append("-")
        left_cordinates[i].append(all_features[i][0][15:j])
        right_cordinates[i].append(all_features[i][0][j+2:-1])
    else:
        strands[i].append("+")
        left_cordinates[i].append(all_features[i][0][4:j])
        right_cordinates[i].append(all_features[i][0][j+2:])
# print(left_cordinates[2299])
# print(right_cordinates[2299])

for i, j in enumerate(all_features):
    for k, l in enumerate(j):
        if l[:11] == "/locus_tag=" and len(locus_tags[i]) == 0:
            locus_tags[i].append(l[11:])
        if j[0][4:] in l and k!=0 and len(features[i]) == 0:
            features[i].append(l[:len(l)-len(j[0][4:])])
# for i in range(len(features)):
#     if len(features[i]) <1:
#         print(i)
for i in range(len(features)):
    if len(features[i]) < 1:
        features[i].append("N/A")
# print(features[2299])


for i in range(len(locus_tags)):
    if len(locus_tags[i]) < 1:
        locus_tags[i].append("N/A")

# for i in range(len(strands)):
#     if len(strands[i]) < 1:
#         strands[i].append("N/A")
#
# for i in range(len(left_cordinates)):
#     if len(left_cordinates[i]) < 1:
#         left_cordinates[i].append("N/A")
#
# for i in range(len(right_cordinates)):
#     if len(right_cordinates[i]) < 1:
#         right_cordinates[i].append("N/A")

with open("left_cordinates_file.pkl", "wb") as left_cordinates_file:
    pickle.dump(left_cordinates, left_cordinates_file)

with open("right_cordinates_file.pkl", "wb") as right_cordinates_file:
    pickle.dump(right_cordinates, right_cordinates_file)

with open("features_file.pkl", "wb") as features_file:
    pickle.dump(features, features_file)

with open("locus_tags_file.pkl", "wb") as locus_tags_file:
    pickle.dump(locus_tags, locus_tags_file)

with open("strands_file.pkl", "wb") as strands_file:
    pickle.dump(strands, strands_file)


# print(strands)
# print(locus_tags)
# print(features)
# print(left_cordinates[0:5])
# print(int(left_cordinates[1][0]))
# print(right_cordinates[0:5])
# print(int(right_cordinates[1][0]))

# for i in  range (0, gene_count, 1):
#     if len(features[i]) != 1:
#         print(all_features[i])
#         print(all_features[i+1])

# for i in range(len(features)):
#     if len(features[i]) < 1:
#         print(features[i])


######################################################################################################################################################################################################################################################################

                                       #file 2 - operon file

#####################################################################################################################################################################################################################################################################
F_operon = open(operon)
f_operon = F_operon.readlines()
fields_operon=[]
for line in f_operon:
    field_operon = line.strip("\n").split()
    fields_operon.append(field_operon)
# print(fields_operon)
# print(fields_operon[0])


operon_nos = [[] for i in range(gene_count)]
operon_left_cordinates = [[] for i in range(gene_count)]
operon_right_cordinates = [[] for i in range(gene_count)]

for i in range(gene_count):
    for j in range(0, len(fields_operon), 1):
        if (strands[i][0] == '+' and
            fields_operon[j][3]  == '+' and
            int(left_cordinates[i][0]) >= int(fields_operon[j][1]) and
            int(right_cordinates[i][0]) <= int(fields_operon[j][2])):
            operon_nos[i].append(fields_operon[j][0])
            operon_left_cordinates[i].append(fields_operon[j][1])
            operon_right_cordinates[i].append(fields_operon[j][2])

        if (strands[i][0] == '-' and
            fields_operon[j][3]  == '-' and
            int(left_cordinates[i][0]) >= int(fields_operon[j][1]) and
            int(right_cordinates[i][0]) <= int(fields_operon[j][2])):
            operon_nos[i].append(fields_operon[j][0])
            operon_left_cordinates[i].append(fields_operon[j][1])
            operon_right_cordinates[i].append(fields_operon[j][2])


for i in range(gene_count):
    if len(operon_nos[i]) < 1:
        operon_nos[i].append("N/A")
        operon_left_cordinates[i].append("N/A")
        operon_right_cordinates[i].append("N/A")


# print(locus_tags[9][0])
# print(left_cordinates[9][0])
# print(right_cordinates[9][0])
# print(strands[9][0])

# print(fields_operon[2257][0])
# print(fields_operon[2257][1])
# print(fields_operon[2257][2])
# print(fields_operon[2257][3])

# print(operon_nos[9][0])
# print(operon_left_cordinates[9][0])
# print(operon_right_cordinates[9][0])

# print(operon_nos)
# print(operon_left_cordinates)
# print(operon_right_cordinates)
# print(len(operon_nos))
# print(len(operon_left_cordinates))
# print(len(operon_right_cordinates))


with open("operon_nos_file.pkl", "wb") as operon_nos_file:
    pickle.dump(operon_nos, operon_nos_file)

with open("operon_left_cordinates_file.pkl", "wb") as operon_left_cordinates_file:
    pickle.dump(operon_left_cordinates, operon_left_cordinates_file)

with open("operon_right_cordinates_file.pkl", "wb") as operon_right_cordinates_file:
    pickle.dump(operon_right_cordinates, operon_right_cordinates_file)











#####################################################################################################################################################################################################################################################################

                                       #file 3 - TSS file

#####################################################################################################################################################################################################################################################################
F2 = open(tss)
f2 = F2.readlines()
fields_tss=[]
for line in f2:
    field_tss = line.strip("\n").split()
    fields_tss.append(field_tss)
# print(fields_tss)
# print(fields_tss[0])

fields2 = []
fields3 = []
for i in range(len(fields_tss)):
    if fields_tss[i][2] == "+":
        fields2.append(fields_tss[i])
    if fields_tss[i][2] == "-":
        fields3.append(fields_tss[i])

# --------------------------sorted function-----------------------
#https://docs.python.org/3/howto/sorting.html
# final_2 = sorted(final, key=itemgetter(5))

#--------------------------sort function-------------------------
#https://www.geeksforgeeks.org/python-list-sort/
#https://www.programiz.com/python-programming/methods/list/sort
def sortSecond(val):
    return int(val[0])
fields2.sort(key = sortSecond)
fields3.sort(key = sortSecond, reverse = True)

# print(fields2)
# print(fields3)














#####################################################################################################################################################################################################################################################################

                                       # gene sequences based on TSS wrt operon info when you need both to incorporate operon and gene

#####################################################################################################################################################################################################################################################################

# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                        # for operon_true
# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# operon_true = [[] for i in range(gene_count)]
# for j in  range (0, len(fields2), 1):
#     for i in range (0, gene_count, 1):
#         if strands[i][0] == "+" and operon_left_cordinates[i][0] != 'N/A' and float(fields2[j][0]) <= float(operon_left_cordinates[i][0]) and float(fields2[j][0]) >= float(operon_left_cordinates[i][0]) - 300: #number changed from 300 to 50 for ncRNAs
#             operon_true[i].append(fields2[j][0])
# for j in  range (0, len(fields3), 1):
#     for i in range (0, gene_count, 1):
#         if strands[i][0] == "-" and operon_right_cordinates[i][0] != 'N/A' and float(fields3[j][0]) >= float(operon_right_cordinates[i][0]) and float(fields3[j][0]) <= float(operon_right_cordinates[i][0]) + 300: # for negative start coordinate is higher and stop is lower so going in direction from left to right-the way wanted in the final version
#             operon_true[i].append(fields3[j][0])
# # print(true)
#
# operon_true_lengths = []
# for i in range (0, len(operon_true), 1):
#     operon_true_lengths.append(len(operon_true[i]))
# # print(operon_true_lengths)
# # print(max(operon_true_lengths))
#
# for i in range (0, len(operon_true), 1):
#     if len(operon_true[i]) < max(operon_true_lengths):
#         for j in range(max(operon_true_lengths)-len(operon_true[i])):
#             operon_true[i].append("x")
# # print(operon_true)
# # print(operon_true[2299])
#
# with open("operon_true_file.pkl", "wb") as operon_true_file:
#     pickle.dump(operon_true, operon_true_file)
#
# # # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# #                                        # operon_distances
# # # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# operon_distances = [[] for i in range(gene_count)]
# # for i in range(max(true_lengths)):
#     # distances[0].append("d")
#
# for i in  range (0, len(operon_true), 1):
#     for j in range (0, len(operon_true[i]), 1):
#         if strands[i][0] == "+" and operon_true[i][j] != "x":
#             operon_distances[i].append(int(operon_left_cordinates[i][0])-int(operon_true[i][j]))
#         if strands[i][0] == "+" and operon_true[i][j] == "x":
#             operon_distances[i].append("x")
#
#         if strands[i][0] == "-" and operon_true[i][j] != "x":
#             operon_distances[i].append(int(operon_true[i][j])-int(operon_right_cordinates[i][0]))
#         if strands[i][0] == "-" and operon_true[i][j] == "x":
#             operon_distances[i].append("x")
# # print(operon_distances)
# # print(operon_distances[2299])
# # print(type(operon_distances[44][1]))
# # print(fields[44])
# # print(len(operon_distances))
# # print(len(operon_true))
# # print(len(fields))
# #
# with open("operon_distances_file.pkl", "wb") as operon_distances_file:
#     pickle.dump(operon_distances, operon_distances_file)
# #
# #
# # # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# #                                        # sequences and start positions for first form based on TSS
# # # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# operon_entire_sequences = [[] for i in range(gene_count)]
# # entire_sequences[0].append("entire_seq")
#
# operon_len_entire_sequences = [[] for i in range(gene_count)]
# # len_entire_sequences[0].append("len_entire_seq")
#
# operon_start_codons = [[] for i in range(gene_count)]
# # start_codons[0].append("start_codon")
#
# operon_fifty_bases = [[] for i in range(gene_count)]
# # fifty_bases[0].append("fifty_bases")
#
# operon_start_positions = [[] for i in range(gene_count)]
# # start_positions[0].append("start_position")
#
# operon_forms = [[] for i in range(gene_count)]
#
# for i in range(0, gene_count, 1):
#     # for j in range (0, max(true_lengths), 1):
#     if strands[i][0] == "+" and operon_true[i][0] != "x":# and int(forms[i]) == j+1:
#         operon_entire_sequences[i].append(seq_2[int(operon_true[i][0])-1:int(operon_right_cordinates[i][0])].upper())
#         operon_len_entire_sequences[i].append(len(operon_entire_sequences[i][0]))
#         operon_start_codons[i].append(seq_2[int(operon_left_cordinates[i][0])-1:int(operon_left_cordinates[i][0])+2].upper())
#         operon_forms[i].append("form_1")
#

        # one more condition refer 1a or below in next section

#         if operon_distances[i][0] >= 25 and len(operon_entire_sequences[i][0]) >= 50:
#             operon_fifty_bases[i].append(seq_2[int(operon_left_cordinates[i][0])-26:int(operon_left_cordinates[i][0])+24].upper())
#             operon_start_positions[i].append(26)
#
#         if operon_distances[i][0] < 25 and len(operon_entire_sequences[i][0]) >= 50:
#             operon_fifty_bases[i].append(seq_2[int(operon_true[i][0])-1:int(operon_true[i][0])+49].upper())
#             operon_start_positions[i].append(operon_distances[i][0]+1)
#
#         if len(operon_entire_sequences[i][0]) < 50:
#             operon_fifty_bases[i].append(operon_entire_sequences[i][0])
#             operon_start_positions[i].append(operon_distances[i][0]+1)
#
#     elif strands[i][0] == "-" and operon_true[i][0] != "x":# and int(forms[i]) == j+1:
#         operon_entire_sequences[i].append(str(Seq(seq_2[int(operon_left_cordinates[i][0])-1:int(operon_true[i][0])]).reverse_complement().upper()))
#         operon_len_entire_sequences[i].append(len(operon_entire_sequences[i][0]))
#         operon_start_codons[i].append(str(Seq(seq_2[int(operon_right_cordinates[i][0])-3:int(operon_right_cordinates[i][0])]).reverse_complement().upper()))
#         operon_forms[i].append("form_1")
#

        # one more condition refer 1a or below in next section



#         if operon_distances[i][0] >= 25 and len(operon_entire_sequences[i][0]) >= 50:
#             operon_fifty_bases[i].append(str(Seq(seq_2[int(operon_right_cordinates[i][0])-25:int(operon_right_cordinates[i][0])+25]).reverse_complement().upper()))
#             operon_start_positions[i].append(26)
#
#         if operon_distances[i][0] < 25 and len(operon_entire_sequences[i][0]) >= 50:
#             operon_fifty_bases[i].append(str(Seq(seq_2[int(operon_true[i][0])-50:int(operon_true[i][0])]).reverse_complement().upper()))
#             operon_start_positions[i].append(operon_distances[i][0]+1)
#
#         if len(operon_entire_sequences[i][0]) < 50:
#             operon_fifty_bases[i].append(operon_entire_sequences[i][0])
#             operon_start_positions[i].append(operon_distances[i][0]+1)
#
#     else:
#         operon_entire_sequences[i].append("x")
#         operon_len_entire_sequences[i].append("x")
#         operon_forms[i].append("x")
#         operon_fifty_bases[i].append("x")
#         operon_start_positions[i].append("x")
#
# # print(operon_fifty_bases)
#
# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                        # sequences and start positions for diff_forms
# #------------------------------------------------------------------------------------------------------------------------------------------------------------------
# for i in  range (0, len(operon_distances), 1):
#     for j in range (1, len(operon_distances[i]), 1): #starting from 1, since the first form already considered before and starting from form2
#         # if strands[i][0] == "+" and len(features[i])>0 and features[i][0] == "CDS" and distances[i][j] != "x" and distances[i][j] < 25: #note: for all genes one form for sure which is form_1 (the most upstream TSS). ALso if (features[i][0] == "CDS") this condition is mentioned, there seems to be no need of (len(features[i])>0) this condition.
#         if strands[i][0] == "+" and operon_distances[i][j] != "x" and operon_distances[i][j] < 25:
#             operon_entire_sequences[i].append(seq_2[int(operon_true[i][j])-1:int(operon_right_cordinates[i][0])].upper())
#             operon_len_entire_sequences[i].append(len(operon_entire_sequences[i][j]))
#             operon_forms[i].append("form_"+str(j+1))
#
#             if operon_distances[i][j] < 25 and len(operon_entire_sequences[i][j]) >= 50:
#                 operon_fifty_bases[i].append(seq_2[int(operon_true[i][j])-1:int(operon_true[i][j])+49].upper())
#                 operon_start_positions[i].append(operon_distances[i][j]+1)
#
#             if len(operon_entire_sequences[i][j]) < 50:
#                 operon_fifty_bases[i].append(operon_entire_sequences[i][j])
#                 operon_start_positions[i].append(operon_distances[i][j]+1)
#
#         # elif strands[i][0] == "-" and len(features[i])>0 and features[i][0] == "CDS" and distances[i][j] != "x" and distances[i][j] < 25: # and int(forms[i]) == j+1:
#         elif strands[i][0] == "-" and operon_distances[i][j] != "x" and operon_distances[i][j] < 25:
#             operon_entire_sequences[i].append(str(Seq(seq_2[int(operon_left_cordinates[i][0])-1:int(operon_true[i][j])]).reverse_complement().upper()))
#             operon_len_entire_sequences[i].append(len(operon_entire_sequences[i][j]))
#             operon_forms[i].append("form_"+str(j+1))
#
#             if operon_distances[i][j] < 25 and len(operon_entire_sequences[i][j]) >= 50:
#                 operon_fifty_bases[i].append(str(Seq(seq_2[int(operon_true[i][j])-50:int(operon_true[i][j])]).reverse_complement().upper()))
#                 operon_start_positions[i].append(operon_distances[i][j]+1)
#
#             if len(operon_entire_sequences[i][j]) < 50:
#                 operon_fifty_bases[i].append(operon_entire_sequences[i][j])
#                 operon_start_positions[i].append(operon_distances[i][j]+1)
#
#         else:
#             operon_entire_sequences[i].append("x")
#             operon_len_entire_sequences[i].append("x")
#             operon_forms[i].append("x")
#             operon_fifty_bases[i].append("x")
#             operon_start_positions[i].append("x")
#
# # print(operon_fifty_bases)
#
# with open("operon_fifty_bases_file.pkl", "wb") as operon_fifty_bases_file:
#     pickle.dump(operon_fifty_bases, operon_fifty_bases_file)
#
# with open("operon_start_positions_file.pkl", "wb") as operon_start_positions_file:
#     pickle.dump(operon_start_positions, operon_start_positions_file)
#
# with open("operon_start_codons_file.pkl", "wb") as operon_start_codons_file:
#     pickle.dump(operon_start_codons, operon_start_codons_file)
#
# with open("operon_forms_file.pkl", "wb") as operon_forms_file:
#     pickle.dump(operon_forms, operon_forms_file)


#####################################################################################################################################################################################################################################################################

                                       # gene sequences based on TSS wrt operon info - but not combining operon and gene info

#####################################################################################################################################################################################################################################################################

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
                                       # for operon_true
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
true = [[] for i in range(gene_count)]
for j in  range (0, len(fields2), 1):
    for i in range (0, gene_count, 1):
        if strands[i][0] == "+" and operon_left_cordinates[i][0] != 'N/A' and float(fields2[j][0]) <= float(operon_left_cordinates[i][0]) and float(fields2[j][0]) >= float(operon_left_cordinates[i][0]) - 300: #number changed from 300 to 50 for ncRNAs
            true[i].append(fields2[j][0])
for j in  range (0, len(fields3), 1):
    for i in range (0, gene_count, 1):
        if strands[i][0] == "-" and operon_right_cordinates[i][0] != 'N/A' and float(fields3[j][0]) >= float(operon_right_cordinates[i][0]) and float(fields3[j][0]) <= float(operon_right_cordinates[i][0]) + 300: # for negative start coordinate is higher and stop is lower so going in direction from left to right-the way wanted in the final version
            true[i].append(fields3[j][0])
# print(true)

true_lengths = []
for i in range (0, len(true), 1):
    true_lengths.append(len(true[i]))
# print(true_lengths)
# print(max(true_lengths))

for i in range (0, len(true), 1):
    if len(true[i]) < max(true_lengths):
        for j in range(max(true_lengths)-len(true[i])):
            true[i].append("x")
# print(true)
# print(true[2299])

with open("true_file.pkl", "wb") as true_file:
    pickle.dump(true, true_file)

# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                        # distances
# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
distances = [[] for i in range(gene_count)]
# for i in range(max(true_lengths)):
    # distances[0].append("d")

for i in  range (0, len(true), 1):
    for j in range (0, len(true[i]), 1):
        if strands[i][0] == "+" and true[i][j] != "x":
            distances[i].append(int(left_cordinates[i][0])-int(true[i][j]))
        if strands[i][0] == "+" and true[i][j] == "x":
            distances[i].append("x")

        if strands[i][0] == "-" and true[i][j] != "x":
            distances[i].append(int(true[i][j])-int(right_cordinates[i][0]))
        if strands[i][0] == "-" and true[i][j] == "x":
            distances[i].append("x")
# print(distances)
# print(distances[2299])
# print(type(distances[44][1]))
# print(fields[44])
# print(len(distances))
# print(len(true))
# print(len(fields))
#
with open("distances_file.pkl", "wb") as distances_file:
    pickle.dump(distances, distances_file)
#
#
# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                        # sequences and start positions for first form based on TSS
# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
entire_sequences = [[] for i in range(gene_count)]
# entire_sequences[0].append("entire_seq")

len_entire_sequences = [[] for i in range(gene_count)]
# len_entire_sequences[0].append("len_entire_seq")

start_codons = [[] for i in range(gene_count)]
# start_codons[0].append("start_codon")

fifty_bases = [[] for i in range(gene_count)]
# fifty_bases[0].append("fifty_bases")

start_positions = [[] for i in range(gene_count)]
# start_positions[0].append("start_position")

forms = [[] for i in range(gene_count)]



#----------------------for 50 bases-------------------------------------------------------------------------------------------------------------------------------------------
# for i in range(0, gene_count, 1):
#     # for j in range (0, max(true_lengths), 1):
#     if strands[i][0] == "+" and true[i][0] != "x":# and int(forms[i]) == j+1:
#         entire_sequences[i].append(seq_2[int(true[i][0])-1:int(right_cordinates[i][0])].upper())
#         len_entire_sequences[i].append(len(entire_sequences[i][0]))
#         start_codons[i].append(seq_2[int(left_cordinates[i][0])-1:int(left_cordinates[i][0])+2].upper())
#         forms[i].append("form_1")
#
#         if len(entire_sequences[i][0]) < 50:
#             fifty_bases[i].append(entire_sequences[i][0])
#             start_positions[i].append(distances[i][0]+1)
#
#         if len(entire_sequences[i][0]) >= 50 and distances[i][0] < 25:
#             fifty_bases[i].append(seq_2[int(true[i][0])-1:int(true[i][0])+49].upper())
#             start_positions[i].append(distances[i][0]+1)
#
#         if len(entire_sequences[i][0]) >= 50 and distances[i][0] >= 25 and len(entire_sequences[i][0])-distances[i][0] < 25:
#             fifty_bases[i].append(seq_2[int(right_cordinates[i][0])-50:int(right_cordinates[i][0])].upper())
#             start_positions[i].append(50-(len(entire_sequences[i][0])-distances[i][0])+1)
#
#         if distances[i][0] >= 25 and len(entire_sequences[i][0])-distances[i][0] >= 25: #len(entire_sequences[i][0]) >= 50:
#             fifty_bases[i].append(seq_2[int(left_cordinates[i][0])-26:int(left_cordinates[i][0])+24].upper())
#             start_positions[i].append(26)
#
#
#     elif strands[i][0] == "-" and true[i][0] != "x":# and int(forms[i]) == j+1:
#         entire_sequences[i].append(str(Seq(seq_2[int(left_cordinates[i][0])-1:int(true[i][0])]).reverse_complement().upper()))
#         len_entire_sequences[i].append(len(entire_sequences[i][0]))
#         start_codons[i].append(str(Seq(seq_2[int(right_cordinates[i][0])-3:int(right_cordinates[i][0])]).reverse_complement().upper()))
#         forms[i].append("form_1")
#
#         if len(entire_sequences[i][0]) < 50:
#             fifty_bases[i].append(entire_sequences[i][0])
#             start_positions[i].append(distances[i][0]+1)
#
#         if len(entire_sequences[i][0]) >= 50 and distances[i][0] < 25:
#             fifty_bases[i].append(str(Seq(seq_2[int(true[i][0])-50:int(true[i][0])]).reverse_complement().upper()))
#             start_positions[i].append(distances[i][0]+1)
#
#         if len(entire_sequences[i][0]) >= 50 and distances[i][0] >= 25 and len(entire_sequences[i][0])-distances[i][0] < 25:
#             fifty_bases[i].append(str(Seq(seq_2[int(left_cordinates[i][0])-1:int(left_cordinates[i][0])+49]).reverse_complement().upper()))
#             start_positions[i].append(50-(len(entire_sequences[i][0])-distances[i][0])+1)
#
#         if distances[i][0] >= 25 and len(entire_sequences[i][0])-distances[i][0] >= 25: #len(entire_sequences[i][0]) >= 50:
#             fifty_bases[i].append(str(Seq(seq_2[int(right_cordinates[i][0])-25:int(right_cordinates[i][0])+25]).reverse_complement().upper()))
#             start_positions[i].append(26)


#----------------------for user input length-------------------------------------------------------------------------------------------------------------------------------------------
for i in range(0, gene_count, 1):
    # for j in range (0, max(true_lengths), 1):
    if strands[i][0] == "+" and true[i][0] != "x":# and int(forms[i]) == j+1:
        entire_sequences[i].append(seq_2[int(true[i][0])-1:int(right_cordinates[i][0])].upper())
        len_entire_sequences[i].append(len(entire_sequences[i][0]))
        start_codons[i].append(seq_2[int(left_cordinates[i][0])-1:int(left_cordinates[i][0])+2].upper())
        forms[i].append("form_1")

        if len(entire_sequences[i][0]) < mrna_length:
            fifty_bases[i].append(entire_sequences[i][0])
            start_positions[i].append(distances[i][0]+1)

        if len(entire_sequences[i][0]) >= mrna_length and distances[i][0] < int(mrna_length/2):
            fifty_bases[i].append(seq_2[int(true[i][0])-1:int(true[i][0])+mrna_length-1].upper())
            start_positions[i].append(distances[i][0]+1)

        if len(entire_sequences[i][0]) >= mrna_length and distances[i][0] >= int(mrna_length/2) and len(entire_sequences[i][0])-distances[i][0] < int(mrna_length/2):
            fifty_bases[i].append(seq_2[int(right_cordinates[i][0])-mrna_length:int(right_cordinates[i][0])].upper())
            start_positions[i].append(mrna_length-(len(entire_sequences[i][0])-distances[i][0])+1)

        if distances[i][0] >= int(mrna_length/2) and len(entire_sequences[i][0])-distances[i][0] >= int(mrna_length/2): #len(entire_sequences[i][0]) >= 50:
            fifty_bases[i].append(seq_2[int(left_cordinates[i][0])-int(mrna_length/2)-1:int(left_cordinates[i][0])+int(mrna_length/2)-1].upper())
            start_positions[i].append(int(mrna_length/2)+1)


    elif strands[i][0] == "-" and true[i][0] != "x":# and int(forms[i]) == j+1:
        entire_sequences[i].append(str(Seq(seq_2[int(left_cordinates[i][0])-1:int(true[i][0])]).reverse_complement().upper()))
        len_entire_sequences[i].append(len(entire_sequences[i][0]))
        start_codons[i].append(str(Seq(seq_2[int(right_cordinates[i][0])-3:int(right_cordinates[i][0])]).reverse_complement().upper()))
        forms[i].append("form_1")

        if len(entire_sequences[i][0]) < mrna_length:
            fifty_bases[i].append(entire_sequences[i][0])
            start_positions[i].append(distances[i][0]+1)

        if len(entire_sequences[i][0]) >= mrna_length and distances[i][0] < int(mrna_length/2):
            fifty_bases[i].append(str(Seq(seq_2[int(true[i][0])-mrna_length:int(true[i][0])]).reverse_complement().upper()))
            start_positions[i].append(distances[i][0]+1)

        if len(entire_sequences[i][0]) >= mrna_length and distances[i][0] >= int(mrna_length/2) and len(entire_sequences[i][0])-distances[i][0] < int(mrna_length/2):
            fifty_bases[i].append(str(Seq(seq_2[int(left_cordinates[i][0])-1:int(left_cordinates[i][0])+mrna_length-1]).reverse_complement().upper()))
            start_positions[i].append(mrna_length-(len(entire_sequences[i][0])-distances[i][0])+1)

        if distances[i][0] >= int(mrna_length/2) and len(entire_sequences[i][0])-distances[i][0] >= int(mrna_length/2): #len(entire_sequences[i][0]) >= 50:
            fifty_bases[i].append(str(Seq(seq_2[int(right_cordinates[i][0])-int(mrna_length/2):int(right_cordinates[i][0])+int(mrna_length/2)]).reverse_complement().upper()))
            start_positions[i].append(int(mrna_length/2)+1)


    else:
        entire_sequences[i].append("x")
        len_entire_sequences[i].append("x")
        forms[i].append("x")
        fifty_bases[i].append("x")
        start_positions[i].append("x")

# print(fifty_bases)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
                                       # sequences and start positions for diff_forms
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------for 50 bases-------------------------------------------------------------------------------------------------------------------------------------------
# for i in  range (0, len(distances), 1):
#     for j in range (1, len(distances[i]), 1): #starting from 1, since the first form already considered before and starting from form2
#         # if strands[i][0] == "+" and len(features[i])>0 and features[i][0] == "CDS" and distances[i][j] != "x" and distances[i][j] < 25: #note: for all genes one form for sure which is form_1 (the most upstream TSS). ALso if (features[i][0] == "CDS") this condition is mentioned, there seems to be no need of (len(features[i])>0) this condition.
#         if strands[i][0] == "+" and distances[i][j] != "x" and distances[i][j] < 25:
#             entire_sequences[i].append(seq_2[int(true[i][j])-1:int(right_cordinates[i][0])].upper())
#             len_entire_sequences[i].append(len(entire_sequences[i][j]))
#             forms[i].append("form_"+str(j+1))
#
#             if distances[i][j] < 25 and len(entire_sequences[i][j]) >= 50:
#                 fifty_bases[i].append(seq_2[int(true[i][j])-1:int(true[i][j])+49].upper())
#                 start_positions[i].append(distances[i][j]+1)
#
#             if len(entire_sequences[i][j]) < 50:
#                 fifty_bases[i].append(entire_sequences[i][j])
#                 start_positions[i].append(distances[i][j]+1)
#
#         # elif strands[i][0] == "-" and len(features[i])>0 and features[i][0] == "CDS" and distances[i][j] != "x" and distances[i][j] < 25: # and int(forms[i]) == j+1:
#         elif strands[i][0] == "-" and distances[i][j] != "x" and distances[i][j] < 25:
#             entire_sequences[i].append(str(Seq(seq_2[int(left_cordinates[i][0])-1:int(true[i][j])]).reverse_complement().upper()))
#             len_entire_sequences[i].append(len(entire_sequences[i][j]))
#             forms[i].append("form_"+str(j+1))
#
#             if distances[i][j] < 25 and len(entire_sequences[i][j]) >= 50:
#                 fifty_bases[i].append(str(Seq(seq_2[int(true[i][j])-50:int(true[i][j])]).reverse_complement().upper()))
#                 start_positions[i].append(distances[i][j]+1)
#
#             if len(entire_sequences[i][j]) < 50:
#                 fifty_bases[i].append(entire_sequences[i][j])
#                 start_positions[i].append(distances[i][j]+1)
#
#         else:
#             entire_sequences[i].append("x")
#             len_entire_sequences[i].append("x")
#             forms[i].append("x")
#             fifty_bases[i].append("x")
#             start_positions[i].append("x")
#
# # print(fifty_bases)


#----------------------for user input length-------------------------------------------------------------------------------------------------------------------------------------------
for i in  range (0, len(distances), 1):
    for j in range (1, len(distances[i]), 1): #starting from 1, since the first form already considered before and starting from form2
        # if strands[i][0] == "+" and len(features[i])>0 and features[i][0] == "CDS" and distances[i][j] != "x" and distances[i][j] < 25: #note: for all genes one form for sure which is form_1 (the most upstream TSS). ALso if (features[i][0] == "CDS") this condition is mentioned, there seems to be no need of (len(features[i])>0) this condition.
        if strands[i][0] == "+" and distances[i][j] != "x" and distances[i][j] < int(mrna_length/2):
            entire_sequences[i].append(seq_2[int(true[i][j])-1:int(right_cordinates[i][0])].upper())
            len_entire_sequences[i].append(len(entire_sequences[i][j]))
            forms[i].append("form_"+str(j+1))

            if distances[i][j] < int(mrna_length/2) and len(entire_sequences[i][j]) >= mrna_length:
                fifty_bases[i].append(seq_2[int(true[i][j])-1:int(true[i][j])+mrna_length-1].upper())
                start_positions[i].append(distances[i][j]+1)

            if len(entire_sequences[i][j]) < mrna_length:
                fifty_bases[i].append(entire_sequences[i][j])
                start_positions[i].append(distances[i][j]+1)

        # elif strands[i][0] == "-" and len(features[i])>0 and features[i][0] == "CDS" and distances[i][j] != "x" and distances[i][j] < 25: # and int(forms[i]) == j+1:
        elif strands[i][0] == "-" and distances[i][j] != "x" and distances[i][j] < int(mrna_length/2):
            entire_sequences[i].append(str(Seq(seq_2[int(left_cordinates[i][0])-1:int(true[i][j])]).reverse_complement().upper()))
            len_entire_sequences[i].append(len(entire_sequences[i][j]))
            forms[i].append("form_"+str(j+1))

            if distances[i][j] < int(mrna_length/2) and len(entire_sequences[i][j]) >= mrna_length:
                fifty_bases[i].append(str(Seq(seq_2[int(true[i][j])-mrna_length:int(true[i][j])]).reverse_complement().upper()))
                start_positions[i].append(distances[i][j]+1)

            if len(entire_sequences[i][j]) < mrna_length:
                fifty_bases[i].append(entire_sequences[i][j])
                start_positions[i].append(distances[i][j]+1)

        else:
            entire_sequences[i].append("x")
            len_entire_sequences[i].append("x")
            forms[i].append("x")
            fifty_bases[i].append("x")
            start_positions[i].append("x")

# print(fifty_bases)





with open("fifty_bases_file.pkl", "wb") as fifty_bases_file:
    pickle.dump(fifty_bases, fifty_bases_file)

with open("start_positions_file.pkl", "wb") as start_positions_file:
    pickle.dump(start_positions, start_positions_file)

with open("start_codons_file.pkl", "wb") as start_codons_file:
    pickle.dump(start_codons, start_codons_file)

with open("forms_file.pkl", "wb") as forms_file:
    pickle.dump(forms, forms_file)








# #####################################################################################################################################################################################################################################################################
#
#                                        # gene sequences based on TSS wrt gene info
#
# #####################################################################################################################################################################################################################################################################
#
# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                        # for true
# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# true = [[] for i in range(gene_count)]
# for j in  range (0, len(fields2), 1):
#     for i in range (0, gene_count, 1):
#         if strands[i][0] == "+" and float(fields2[j][0]) <= float(left_cordinates[i][0]) and float(fields2[j][0]) >= float(left_cordinates[i][0]) - 300: #number changed from 300 to 50 for ncRNAs
#             true[i].append(fields2[j][0])
# for j in  range (0, len(fields3), 1):
#     for i in range (0, gene_count, 1):
#         if strands[i][0] == "-" and float(fields3[j][0]) >= float(right_cordinates[i][0]) and float(fields3[j][0]) <= float(right_cordinates[i][0]) + 300: # for negative start coordinate is higher and stop is lower so going in direction from left to right-the way wanted in the final version
#             true[i].append(fields3[j][0])
# # print(true)
#
# true_lengths = []
# for i in range (0, len(true), 1):
#     true_lengths.append(len(true[i]))
# # print(lengths3)
# # print(max(lengths3))
#
# for i in range (0, len(true), 1):
#     if len(true[i]) < max(true_lengths):
#         for j in range(max(true_lengths)-len(true[i])):
#             true[i].append("x")
# # print(true[2299])
#
# with open("true_file.pkl", "wb") as true_file:
#     pickle.dump(true, true_file)
#
# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                        # distances
# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# distances = [[] for i in range(gene_count)]
# # for i in range(max(true_lengths)):
#     # distances[0].append("d")
#
# for i in  range (0, len(true), 1):
#     for j in range (0, len(true[i]), 1):
#         if strands[i][0] == "+" and true[i][j] != "x":
#             distances[i].append(int(left_cordinates[i][0])-int(true[i][j]))
#         if strands[i][0] == "+" and true[i][j] == "x":
#             distances[i].append("x")
#
#         if strands[i][0] == "-" and true[i][j] != "x":
#             distances[i].append(int(true[i][j])-int(right_cordinates[i][0]))
#         if strands[i][0] == "-" and true[i][j] == "x":
#             distances[i].append("x")
#
# with open("distances_file.pkl", "wb") as distances_file:
#     pickle.dump(distances, distances_file)
# # print(distances[2299])
# # print(type(distances[44][1]))
# # print(fields[44])
# # print(len(distances))
# # print(len(true))
# # print(len(fields))
#
#
#
# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                        # sequences and start positions for first form based on TSS
# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# entire_sequences = [[] for i in range(gene_count)]
# # entire_sequences[0].append("entire_seq")
#
# len_entire_sequences = [[] for i in range(gene_count)]
# # len_entire_sequences[0].append("len_entire_seq")
#
# start_codons = [[] for i in range(gene_count)]
# # start_codons[0].append("start_codon")
#
# fifty_bases = [[] for i in range(gene_count)]
# # fifty_bases[0].append("fifty_bases")
#
# start_positions = [[] for i in range(gene_count)]
# # start_positions[0].append("start_position")
#
# forms = [[] for i in range(gene_count)]
#
# for i in range(0, gene_count, 1):
#     # for j in range (0, max(true_lengths), 1):
#     if strands[i][0] == "+" and true[i][0] != "x":# and int(forms[i]) == j+1:
#         entire_sequences[i].append(seq_2[int(true[i][0])-1:int(right_cordinates[i][0])].upper())
#         len_entire_sequences[i].append(len(entire_sequences[i][0]))
#         start_codons[i].append(seq_2[int(left_cordinates[i][0])-1:int(left_cordinates[i][0])+2].upper())
#         forms[i].append("form_1")
#
        # if distances[i][0] >= 25 and len(entire_sequences[i][0])-distances[i][0] < 25:
            # fifty_bases[i].append(seq_2[int(right_cordinates[i][0])-50:int(right_cordinates[i][0])].upper())
            # start_positions[i].append(50-(len(entire_sequences[i][0])-distances[i][0])+1)

#         if distances[i][0] >= 25 and len(entire_sequences[i][0]) >= 50:
#             fifty_bases[i].append(seq_2[int(left_cordinates[i][0])-26:int(left_cordinates[i][0])+24].upper())
#             start_positions[i].append(26)
#
#         if distances[i][0] < 25 and len(entire_sequences[i][0]) >= 50:
#             fifty_bases[i].append(seq_2[int(true[i][0])-1:int(true[i][0])+49].upper())
#             start_positions[i].append(distances[i][0]+1)
#
#         if len(entire_sequences[i][0]) < 50:
#             fifty_bases[i].append(entire_sequences[i][0])
#             start_positions[i].append(distances[i][0]+1)
#
#     elif strands[i][0] == "-" and true[i][0] != "x":# and int(forms[i]) == j+1:
#         entire_sequences[i].append(str(Seq(seq_2[int(left_cordinates[i][0])-1:int(true[i][0])]).reverse_complement().upper()))
#         len_entire_sequences[i].append(len(entire_sequences[i][0]))
#         start_codons[i].append(str(Seq(seq_2[int(right_cordinates[i][0])-3:int(right_cordinates[i][0])]).reverse_complement().upper()))
#         forms[i].append("form_1")
#
        # if distances[i][0] >= 25 and len(entire_sequences[i][0])-distances[i][0] < 25:
        #     fifty_bases[i].append(str(Seq(seq_2[int(left_cordinates[i][0])-1:int(left_cordinates[i][0])+49]).reverse_complement().upper()))
        #     start_positions[i].append(50-(len(entire_sequences[i][0])-distances[i][0])+1)

#         if distances[i][0] >= 25 and len(entire_sequences[i][0]) >= 50:
#             fifty_bases[i].append(str(Seq(seq_2[int(right_cordinates[i][0])-25:int(right_cordinates[i][0])+25]).reverse_complement().upper()))
#             start_positions[i].append(26)
#
#         if distances[i][0] < 25 and len(entire_sequences[i][0]) >= 50:
#             fifty_bases[i].append(str(Seq(seq_2[int(true[i][0])-50:int(true[i][0])]).reverse_complement().upper()))
#             start_positions[i].append(distances[i][0]+1)
#
#         if len(entire_sequences[i][0]) < 50:
#             fifty_bases[i].append(entire_sequences[i][0])
#             start_positions[i].append(distances[i][0]+1)
#
#     else:
#         entire_sequences[i].append("x")
#         len_entire_sequences[i].append("x")
#         forms[i].append("x")
#         fifty_bases[i].append("x")
#         start_positions[i].append("x")
#
#
# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                        # sequences and start positions for diff_forms
# #--------------------------------------------------------------------------------------------------------------------------------------------------------
# for i in  range (0, len(distances), 1):
#     for j in range (1, len(distances[i]), 1): #starting from 1, since the first form already considered before and starting from form2
#         # if strands[i][0] == "+" and len(features[i])>0 and features[i][0] == "CDS" and distances[i][j] != "x" and distances[i][j] < 25: #note: for all genes one form for sure which is form_1 (the most upstream TSS). ALso if (features[i][0] == "CDS") this condition is mentioned, there seems to be no need of (len(features[i])>0) this condition.
#         if strands[i][0] == "+" and distances[i][j] != "x" and distances[i][j] < 25:
#             entire_sequences[i].append(seq_2[int(true[i][j])-1:int(right_cordinates[i][0])].upper())
#             len_entire_sequences[i].append(len(entire_sequences[i][j]))
#             forms[i].append("form_"+str(j+1))
#
#             if distances[i][j] < 25 and len(entire_sequences[i][j]) >= 50:
#                 fifty_bases[i].append(seq_2[int(true[i][j])-1:int(true[i][j])+49].upper())
#                 start_positions[i].append(distances[i][j]+1)
#
#             if len(entire_sequences[i][j]) < 50:
#                 fifty_bases[i].append(entire_sequences[i][j])
#                 start_positions[i].append(distances[i][j]+1)
#
#         # elif strands[i][0] == "-" and len(features[i])>0 and features[i][0] == "CDS" and distances[i][j] != "x" and distances[i][j] < 25: # and int(forms[i]) == j+1:
#         elif strands[i][0] == "-" and distances[i][j] != "x" and distances[i][j] < 25:
#             entire_sequences[i].append(str(Seq(seq_2[int(left_cordinates[i][0])-1:int(true[i][j])]).reverse_complement().upper()))
#             len_entire_sequences[i].append(len(entire_sequences[i][j]))
#             forms[i].append("form_"+str(j+1))
#
#             if distances[i][j] < 25 and len(entire_sequences[i][j]) >= 50:
#                 fifty_bases[i].append(str(Seq(seq_2[int(true[i][j])-50:int(true[i][j])]).reverse_complement().upper()))
#                 start_positions[i].append(distances[i][j]+1)
#
#             if len(entire_sequences[i][j]) < 50:
#                 fifty_bases[i].append(entire_sequences[i][j])
#                 start_positions[i].append(distances[i][j]+1)
#
#         else:
#             entire_sequences[i].append("x")
#             len_entire_sequences[i].append("x")
#             forms[i].append("x")
#             fifty_bases[i].append("x")
#             start_positions[i].append("x")
#
# # print(fifty_bases)
#
# with open("fifty_bases_file.pkl", "wb") as fifty_bases_file:
#     pickle.dump(fifty_bases, fifty_bases_file)
#
# with open("start_positions_file.pkl", "wb") as start_positions_file:
#     pickle.dump(start_positions, start_positions_file)
#
# with open("start_codons_file.pkl", "wb") as start_codons_file:
#     pickle.dump(start_codons, start_codons_file)
#
# with open("forms_file.pkl", "wb") as forms_file:
#     pickle.dump(forms, forms_file)
#
#











#####################################################################################################################################################################################################################################################################

                                       # gene sequences as leaderless

#####################################################################################################################################################################################################################################################################







#####################################################################################################################################################################################################################################################################

                                       # gene sequences as leadered

#####################################################################################################################################################################################################################################################################












# # # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# # # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# # #------------------------------other_info----------------------------------------------------------------------------------------------------------------------
# # # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# # # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# # TSS_types = [[] for i in range(len(fields))]
# # TSS_types[0].append("TSS_type")
# # numbers_of_TSS = [[] for i in range(len(fields))]
# # numbers_of_TSS[0].append("number_of_TSS")
# # internals = [[] for i in range(len(fields))]
# # internals[0].append("internal")
# # flags = [[] for i in range(len(fields))]
# # flags[0].append("flag")
# # non_coding_RNAs = [[] for i in range(len(fields))]
# # non_coding_RNAs[0].append("non_coding_RNA")
# # RNAs = [[] for i in range(len(fields))]
# # RNAs[0].append("RNA")
# #
# # for i in range(1, len(fields), 1):
# #     TSS_types[i].append("NA")
# #     numbers_of_TSS[i].append("NA")
# #     internals[i].append("NA")
# #     flags[i].append("NA")
# #     non_coding_RNAs[i].append("NA")
# #     RNAs[i].append("NA")



#####################################################################################################################################################################################################################################################################

                                       # end

#####################################################################################################################################################################################################################################################################
