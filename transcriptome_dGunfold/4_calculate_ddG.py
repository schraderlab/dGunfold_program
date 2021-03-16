#####################################################################################################################################################################################################################################################################

                                       # libraries

#####################################################################################################################################################################################################################################################################
import pickle








#####################################################################################################################################################################################################################################################################

                                       # sys_argv

#####################################################################################################################################################################################################################################################################















#####################################################################################################################################################################################################################################################################

                                       # file 1 - sequences, structure, dGrefs

#####################################################################################################################################################################################################################################################################
F = open ("2_sequences_vienna_output.txt")
f= F.readlines()
# print(f)
fields = []
dG_refs = []
sequences = []
dot_bracket_structures = []

# for line in f:
#     fields.append(field)
#     a = field[0]
#     b = field[1]
#     c = field[2]
#     d = field[3]
#     e = field[4]
#     f = field[5]
#     g = field[6]
#     h = field[7]
# print(f)
# print(fields)


for line in f:
    field = line.strip("\n").replace(" ( ", "space").replace(" (", "space").split("space") #-------------------since if the value is single digit, there is a space which may ruin the indexing. therefore first figure out a way to avoid that-------------
    fields.append(field)

# print(fields)

for i in range(1, len(f), 2):
    dG_refs.append(fields[i][1].replace(")", ""))
    dot_bracket_structures.append(fields[i][0])

# print(dG_refs)
# print(len(dG_refs))
# for dG_ref in dG_refs:
#     print(dG_ref)

# print(dot_bracket_structures)
# print(len(dot_bracket_structures))

for i in range(0, len(f), 2):
    sequences.append(fields[i][0])
# print(sequences)
# print(len(sequences))



#####################################################################################################################################################################################################################################################################

                                       # file 2 - dGinits, dGunfolds

#####################################################################################################################################################################################################################################################################
F2 = open ("4_sequences_with_constraint_vienna_output.txt")
f2= F2.readlines()
# print(f2)
fields2 = []
dG_consts = []

for line in f2:
    field2 = line.strip("\n").replace(" ( ", "space").replace(" (", "space").split("space")
    fields2.append(field2)

# print(fields2)

for i in range(1, len(f2), 2):
    dG_consts.append(fields2[i][1].replace(")", ""))

# print(dG_consts)
# print(len(dG_consts))
# for dG_const in dG_consts:
#     print(dG_const)


#---------------------------dGunfold----------------------------------------------------
dG_unfolds = []
for i in range(0, len(dG_refs), 1):
    dG_unfolds.append(float(dG_consts[i])-float(dG_refs[i]))

# for dG_unfold in dG_unfolds:
#     print(round(dG_unfold, 2))



#####################################################################################################################################################################################################################################################################

                                       # contraints pkl file

#####################################################################################################################################################################################################################################################################
with open("constraints_file.pkl", "rb") as constraints_file:
    constraints = pickle.load(constraints_file)
# print(constraints)


# F3 = open ("4_constraints.txt")
# f3= F3.readlines()
# # print(f3)
# fields3 = []
# constraints = []
# for line in f3:
#     field3 = line.strip("\n")
#     constraints.append(field3)

# print(constraints)
# print(len(constraints))


#####################################################################################################################################################################################################################################################################

                                       # genbank features pkl file

#####################################################################################################################################################################################################################################################################
with open("left_cordinates_file.pkl", "rb") as left_cordinates_file:
    left_cordinates = pickle.load(left_cordinates_file)

with open("right_cordinates_file.pkl", "rb") as right_cordinates_file:
    right_cordinates = pickle.load(right_cordinates_file)

with open("features_file.pkl", "rb") as features_file:
    features = pickle.load(features_file)

with open("locus_tags_file.pkl", "rb") as locus_tags_file:
    locus_tags = pickle.load(locus_tags_file)

with open("strands_file.pkl", "rb") as strands_file:
    strands = pickle.load(strands_file)



#####################################################################################################################################################################################################################################################################

                                       # tss features pkl file

#####################################################################################################################################################################################################################################################################

with open("true_file.pkl", "rb") as true_file:
    true = pickle.load(true_file)
# print(true)
# print(len(true))

with open("distances_file.pkl", "rb") as distances_file:
    distances = pickle.load(distances_file)
# print(distances)
# print(len(distances))

with open("start_codons_file.pkl", "rb") as start_codons_file:
    start_codons = pickle.load(start_codons_file)

with open("forms_file.pkl", "rb") as forms_file:
    forms = pickle.load(forms_file)

# for i in range(len(features)):
#     if len(features[i]) < 1:
#         print(i, features[i])
# print(true[2299])

# print(left_cordinates)
# print(len(left_cordinates))

#####################################################################################################################################################################################################################################################################

                                       # operon features pkl file

#####################################################################################################################################################################################################################################################################

import os.path
from os import path
# print(path.exists('operon_nos_file.pkl'))
# print(os.path.isfile('operon_nos_file.pkl'))
# print(path.isfile('operon_nos_file.pkl'))

if path.exists('operon_nos_file.pkl'):
    with open("operon_nos_file.pkl", "rb") as operon_nos_file:
        operon_nos = pickle.load(operon_nos_file)
# print(operon_nos)
# print(len(operon_nos))

if path.exists('operon_left_cordinates_file.pkl'):
    with open("operon_left_cordinates_file.pkl", "rb") as operon_left_cordinates_file:
        operon_left_cordinates = pickle.load(operon_left_cordinates_file)

if path.exists('operon_right_cordinates_file.pkl'):
    with open("operon_right_cordinates_file.pkl", "rb") as operon_right_cordinates_file:
        operon_right_cordinates = pickle.load(operon_right_cordinates_file)
# print(operon_right_cordinates)

if path.exists('operon_true_file.pkl'):
    with open("operon_true_file.pkl", "rb") as operon_true_file:
        operon_true = pickle.load(operon_true_file)

if path.exists('operon_distances_file.pkl'):
    with open("operon_distances_file.pkl", "rb") as operon_distances_file:
        operon_distances = pickle.load(operon_distances_file)

if path.exists('operon_fifty_bases_file.pkl'):
    with open("operon_fifty_bases_file.pkl", "rb") as operon_fifty_bases_file:
        operon_fifty_bases = pickle.load(operon_fifty_bases_file)

if path.exists('operon_start_positions_file.pkl'):
    with open("operon_start_positions_file.pkl", "rb") as operon_start_positions_file:
        operon_start_positions = pickle.load(operon_start_positions_file)

if path.exists('operon_start_codons_file.pkl'):
    with open("operon_start_codons_file.pkl", "rb") as operon_start_codons_file:
        operon_start_codons = pickle.load(operon_start_codons_file)

if path.exists('operon_forms_file.pkl'):
    with open("operon_forms_file.pkl", "rb") as operon_forms_file:
        operon_forms = pickle.load(operon_forms_file)


# print(locus_tags[9][0])
# print(left_cordinates[9][0])
# print(right_cordinates[9][0])
# print(strands[9][0])

# print(operon_nos[9][0])
# print(operon_left_cordinates[9][0])
# print(operon_right_cordinates[9][0])









#------------------------------final printing------------------------------------------------------------
# dGunfolds_final = [[] for i in range(len(dG_unfolds))]
# for i in range(len(dG_unfolds)):
#     # print(sequences[i], dot_bracket_structures[i], constraints[i], round(float(dG_refs[i]), 2), round(float(dG_consts[i]), 2), round(float(dG_unfolds[i]), 2), sep="\t")
#     dGunfolds_final[i].append(sequences[i])
#     dGunfolds_final[i].append(dot_bracket_structures[i])
#     dGunfolds_final[i].append(constraints[i])
#     dGunfolds_final[i].append(round(float(dG_refs[i]), 2))
#     dGunfolds_final[i].append(round(float(dG_consts[i]), 2))
#     dGunfolds_final[i].append(round(float(dG_unfolds[i]), 2))



# properties_count = 0
# for i in range(len(true)):
#     for j in range(len(true[i])):
#         if j == 0 and true[i][j] != "x":
#             properties_count += 1
#         elif j > 0 and true[i][j] != "x" and distances[i][j] < 25 and features[i][0] == "CDS":
#             properties_count += 1
# print(properties_count)
# print(len(dG_unfolds))

# properties = [[] for i in range(properties_count)]








#####################################################################################################################################################################################################################################################################

                                       # genbank and operon features new dataset - this is needed since for multiforms all these data needs to be repeated

#####################################################################################################################################################################################################################################################################

left_cordinates_new = []
right_cordinates_new = []
features_new = []
locus_tags_new = []
strands_new = []
true_new = []
distances_new = []
start_codons_new = []
forms_new = []
operon_nos_new = []
operon_left_cordinates_new = []
operon_right_cordinates_new = []



for i in range(len(true)):
    for j in range(len(true[i])):
        # if (j == 0 and true[i][j] != "x") or (j > 0 and true[i][j] != "x" and distances[i][j] < 25 and features[i][0] == "CDS"): # these conditions should exactly meet the conditions for selecting sequences
        if (j == 0 and true[i][j] != "x") or (j > 0 and true[i][j] != "x" and distances[i][j] < 25):
            locus_tags_new.append(locus_tags[i][0])
            features_new.append(features[i][0])
            strands_new.append(strands[i][0])
            left_cordinates_new.append(left_cordinates[i][0])
            right_cordinates_new.append(right_cordinates[i][0])
            true_new.append(true[i][j])
            distances_new.append(distances[i][j])
            start_codons_new.append(start_codons[i][0])
            forms_new.append(forms[i][j])
            if path.exists('operon_nos_file.pkl'):
                operon_nos_new.append(operon_nos[i][0])
                operon_left_cordinates_new.append(operon_left_cordinates[i][0])
                operon_right_cordinates_new.append(operon_right_cordinates[i][0])
# print(features_new)
# print(len(features_new))
# print(distances_new)
# print(len(distances_new))

with open("left_cordinates_new_file.pkl", "wb") as left_cordinates_new_file:
    pickle.dump(left_cordinates_new, left_cordinates_new_file)

with open("right_cordinates_new_file.pkl", "wb") as right_cordinates_new_file:
    pickle.dump(right_cordinates_new, right_cordinates_new_file)

with open("features_new_file.pkl", "wb") as features_new_file:
    pickle.dump(features_new, features_new_file)

with open("locus_tags_new_file.pkl", "wb") as locus_tags_new_file:
    pickle.dump(locus_tags_new, locus_tags_new_file)

with open("strands_new_file.pkl", "wb") as strands_new_file:
    pickle.dump(strands_new, strands_new_file)

with open("true_new_file.pkl", "wb") as true_new_file:
    pickle.dump(true_new, true_new_file)

with open("distances_new_file.pkl", "wb") as distances_new_file:
    pickle.dump(distances_new, distances_new_file)

with open("start_codons_new_file.pkl", "wb") as start_codons_new_file:
    pickle.dump(start_codons_new, start_codons_new_file)

with open("forms_new_file.pkl", "wb") as forms_new_file:
    pickle.dump(forms_new, forms_new_file)


# # need to put start codons in print function
if path.exists('operon_nos_file.pkl'):
    print('operon_no', 'operon_left_cordinate', 'operon_right_cordinate', 'locus_tag', 'feature', 'strand', 'form', 'gene_left_cordinate', 'gene_right_cordinate', 'TSS_pos', '5_UTR_length', 'sequence', 'dot-bracket_structure', 'constrain', 'dGmRNA', 'dGinit', 'dGunfold')
    for i in range(len(dG_unfolds)):
        # if features_new[i] == "CDS":
        print(operon_nos_new[i], operon_left_cordinates_new[i], operon_right_cordinates_new[i], locus_tags_new[i], features_new[i], strands_new[i], forms_new[i], left_cordinates_new[i], right_cordinates_new[i], true_new[i], distances_new[i],
              sequences[i], dot_bracket_structures[i], constraints[i], round(float(dG_refs[i]), 2), round(float(dG_consts[i]), 2), round(float(dG_unfolds[i]), 2), sep = "\t")

else:
    print('locus_tag', 'feature', 'strand', 'form', 'gene_left_cordinate', 'gene_right_cordinate', 'TSS_pos', '5_UTR_length', 'sequence', 'dot-bracket_structure', 'constrain', 'dGmRNA', 'dGinit', 'dGunfold')
    for i in range(len(dG_unfolds)):
        # if features_new[i] == "CDS":
        print(locus_tags_new[i], features_new[i], strands_new[i], forms_new[i], left_cordinates_new[i], right_cordinates_new[i], true_new[i], distances_new[i],
              sequences[i], dot_bracket_structures[i], constraints[i], round(float(dG_refs[i]), 2), round(float(dG_consts[i]), 2), round(float(dG_unfolds[i]), 2), sep = "\t")
