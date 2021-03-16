#####################################################################################################################################################################################################################################################################

                                       # libraries

#####################################################################################################################################################################################################################################################################
from itertools import repeat
import sys
import pickle
from Bio.Data import CodonTable
from Bio.Seq import Seq
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]










#####################################################################################################################################################################################################################################################################

                                       # sys_argv

#####################################################################################################################################################################################################################################################################
genbank = sys.argv[1]
operon = sys.argv[2]
mrna_length = int(sys.argv[3])











#####################################################################################################################################################################################################################################################################

                                       # file 1-extracting features from genbank file

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
# print(all_features)
# print(len(all_features))
# print(all_features[0])
# print(all_features[1])

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


for i, j in enumerate(all_features):
    for k, l in enumerate(j):
        if l[:11] == "/locus_tag=" and len(locus_tags[i]) == 0:
            locus_tags[i].append(l[11:])
        if j[0][4:] in l and k!=0 and len(features[i]) == 0:
            features[i].append(l[:len(l)-len(j[0][4:])])

for i in range(len(features)):
    if len(features[i]) < 1:
        features[i].append("N/A")
# print(features[2299])
# print(features)
# print(len(features))

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
# for i in range(len(features)):
#     if len(features[i]) < 1:
#         print(features[i])
# print(features)
# print(len(features))
# print(left_cordinates[0:5])
# print(int(left_cordinates[1][0]))

# print(right_cordinates[0:5])
# print(int(right_cordinates[1][0]))

# for i in  range (0, gene_count, 1):
#     if len(features[i]) != 1:
#         print(all_features[i])
#         print(all_features[i+1])






#####################################################################################################################################################################################################################################################################

                                       # file 2 - operon file

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

                                       # file 3 - TSS file

#####################################################################################################################################################################################################################################################################










#####################################################################################################################################################################################################################################################################

                                       # gene sequences based on TSS wrt operon info

#####################################################################################################################################################################################################################################################################







#####################################################################################################################################################################################################################################################################

                                       # gene sequences based on TSS wrt gene info

#####################################################################################################################################################################################################################################################################








#####################################################################################################################################################################################################################################################################

                                       # gene sequences as leaderless

#####################################################################################################################################################################################################################################################################
#
# #-----------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                       # for true and distances
# #-----------------------------------------------------------------------------------------------------------------------------------------------------------------
#
# true = [["N/A"] for i in range(gene_count)]
# distances = [[0] for i in range(gene_count)]
#
# # print(true)
# # print(len(true))
#
# with open("true_file.pkl", "wb") as true_file:
#     pickle.dump(true, true_file)
#
# with open("distances_file.pkl", "wb") as distances_file:
#     pickle.dump(distances, distances_file)
#
#
#
#
#
#
#
#
# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                        # sequences for all genes from 1st position (all as leaderless)
# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# sequences = [[] for i in range(gene_count)]
# for i in range(gene_count):
#     if strands[i][0] == "+":
#         sequences[i].append(seq_2[int(left_cordinates[i][0])-1:int(right_cordinates[i][0])].upper())
#     if strands[i][0] == "-":
#         sequences[i].append(str(Seq(seq_2[int(left_cordinates[i][0])-1:int(right_cordinates[i][0])]).reverse_complement().upper()))
# # print(strands[1][0])
# # print(seq_2[int(left_cordinates[1][0])-1:int(right_cordinates[1][0])])
# # print(sequences)
# # print(len(sequences))
#
# fifty_bases_as_leaderless = [[] for i in range(gene_count)]
# for i in range(gene_count):
#     if strands[i][0] == "+" and len(sequences[i][0]) >= 50:
#         fifty_bases_as_leaderless[i].append(seq_2[int(left_cordinates[i][0])-1:int(left_cordinates[i][0])+49].upper())
#     if strands[i][0] == "+" and len(sequences[i][0]) < 50:
#         fifty_bases_as_leaderless[i].append(sequences[i][0])
#
#     if strands[i][0] == "-" and len(sequences[i][0]) >= 50:
#         fifty_bases_as_leaderless[i].append(str(Seq(seq_2[int(right_cordinates[i][0])-50:int(right_cordinates[i][0])]).reverse_complement().upper()))
#     if strands[i][0] == "-" and len(sequences[i][0]) < 50:
#         fifty_bases_as_leaderless[i].append(sequences[i][0])
#
# start_positions_as_leaderless = [[1] for i in range(gene_count)]
# # print(start_positions_as_leaderless)
# # for i in range(gene_count):
# #     # print(len(fifty_bases_as_leaderless[i][0]))
# #     print(fifty_bases_as_leaderless[i][0])
#
# start_codons_as_leaderless = [[] for i in range(gene_count)]
# for i in range(gene_count):
#     if strands[i][0] == "+":
#         start_codons_as_leaderless[i].append(seq_2[int(left_cordinates[i][0])-1:int(left_cordinates[i][0])+2].upper())
#     if strands[i][0] == "-":
#         start_codons_as_leaderless[i].append(str(Seq(seq_2[int(right_cordinates[i][0])-3:int(right_cordinates[i][0])]).reverse_complement().upper()))
#
# forms = [["form_1"] for i in range(gene_count)]
#
#
# with open("fifty_bases_file.pkl", "wb") as fifty_bases_file:
#     pickle.dump(fifty_bases_as_leaderless, fifty_bases_file)
#
# with open("start_positions_file.pkl", "wb") as start_positions_file:
#     pickle.dump(start_positions_as_leaderless, start_positions_file)
#
# with open("start_codons_file.pkl", "wb") as start_codons_file:
#     pickle.dump(start_codons_as_leaderless, start_codons_file)
#
# with open("forms_file.pkl", "wb") as forms_file:
#     pickle.dump(forms, forms_file)











#####################################################################################################################################################################################################################################################################

                                       # gene sequences as leadered

#####################################################################################################################################################################################################################################################################

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
                                      # for true and distances
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------

true = [["N/A"] for i in range(gene_count)]
distances = [[0] for i in range(gene_count)]

# print(true)
# print(len(true))

with open("true_file.pkl", "wb") as true_file:
    pickle.dump(true, true_file)

with open("distances_file.pkl", "wb") as distances_file:
    pickle.dump(distances, distances_file)








# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
                                       # sequences for all genes (all as leadered)
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
sequences = [[] for i in range(gene_count)]
for i in range(gene_count):
    if strands[i][0] == "+":
        sequences[i].append(seq_2[int(left_cordinates[i][0])-1:int(right_cordinates[i][0])].upper())
    if strands[i][0] == "-":
        sequences[i].append(str(Seq(seq_2[int(left_cordinates[i][0])-1:int(right_cordinates[i][0])]).reverse_complement().upper()))
# print(strands[1][0])
# print(seq_2[int(left_cordinates[1][0])-1:int(right_cordinates[1][0])])
# print(sequences)
# print(len(sequences))

fifty_bases = [[] for i in range(gene_count)]

start_positions = [[] for i in range(gene_count)]
# start_positions[0].append("start_position")

#----------------------for 50 bases-------------------------------------------------------------------------------------------------------------------------------------------
# for i in range(gene_count):
#     if strands[i][0] == "+" and len(sequences[i][0]) >= 25:
#         fifty_bases[i].append(seq_2[int(left_cordinates[i][0])-26:int(left_cordinates[i][0])+24].upper())
#         start_positions[i].append(26)
#     if strands[i][0] == "+" and len(sequences[i][0]) < 25:
#         fifty_bases[i].append(seq_2[int(right_cordinates[i][0])-50:int(right_cordinates[i][0])].upper())
#         start_positions[i].append(50-len(sequences[i][0])+1)
#
#     if strands[i][0] == "-" and len(sequences[i][0]) >= 25:
#         fifty_bases[i].append(str(Seq(seq_2[int(right_cordinates[i][0])-25:int(right_cordinates[i][0])+25]).reverse_complement().upper()))
#         start_positions[i].append(26)
#     if strands[i][0] == "-" and len(sequences[i][0]) < 25:
#         fifty_bases[i].append(str(Seq(seq_2[int(left_cordinates[i][0])-1:int(left_cordinates[i][0])+49]).reverse_complement().upper()))
#         start_positions[i].append(50-len(sequences[i][0])+1)


#----------------------for user input length-------------------------------------------------------------------------------------------------------------------------------------------
for i in range(gene_count):
    if strands[i][0] == "+" and len(sequences[i][0]) >= int(mrna_length/2):
        fifty_bases[i].append(seq_2[int(left_cordinates[i][0])-int(mrna_length/2)-1:int(left_cordinates[i][0])+int(mrna_length/2)-1].upper())
        start_positions[i].append(int(mrna_length/2)+1)
    if strands[i][0] == "+" and len(sequences[i][0]) < int(mrna_length/2):
        fifty_bases[i].append(seq_2[int(right_cordinates[i][0])-mrna_length:int(right_cordinates[i][0])].upper())
        start_positions[i].append(mrna_length-len(sequences[i][0])+1)

    if strands[i][0] == "-" and len(sequences[i][0]) >= int(mrna_length/2):
        fifty_bases[i].append(str(Seq(seq_2[int(right_cordinates[i][0])-int(mrna_length/2):int(right_cordinates[i][0])+int(mrna_length/2)]).reverse_complement().upper()))
        start_positions[i].append(int(mrna_length/2)+1)
    if strands[i][0] == "-" and len(sequences[i][0]) < int(mrna_length/2):
        fifty_bases[i].append(str(Seq(seq_2[int(left_cordinates[i][0])-1:int(left_cordinates[i][0])+mrna_length-1]).reverse_complement().upper()))
        start_positions[i].append(mrna_length-len(sequences[i][0])+1)



start_codons = [[] for i in range(gene_count)]
for i in range(gene_count):
    if strands[i][0] == "+":
        start_codons[i].append(seq_2[int(left_cordinates[i][0])-1:int(left_cordinates[i][0])+2].upper())
    if strands[i][0] == "-":
        start_codons[i].append(str(Seq(seq_2[int(right_cordinates[i][0])-3:int(right_cordinates[i][0])]).reverse_complement().upper()))

forms = [["form_1"] for i in range(gene_count)]


with open("fifty_bases_file.pkl", "wb") as fifty_bases_file:
    pickle.dump(fifty_bases, fifty_bases_file)

with open("start_positions_file.pkl", "wb") as start_positions_file:
    pickle.dump(start_positions, start_positions_file)

with open("start_codons_file.pkl", "wb") as start_codons_file:
    pickle.dump(start_codons, start_codons_file)

with open("forms_file.pkl", "wb") as forms_file:
    pickle.dump(forms, forms_file)



#####################################################################################################################################################################################################################################################################

                                       # end

#####################################################################################################################################################################################################################################################################
