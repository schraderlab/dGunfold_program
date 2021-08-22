#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#----------------------------libraries----------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
import pickle
from operator import itemgetter, attrgetter
#------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------sys_argv files-----------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
F = open ("2_sequences_vienna_output.txt")
f= F.readlines()
# print(f)
fields = []
dG_refs = []
sequences = []
dot_bracket_structures = []

for line in f:
    field = line.strip("\n").replace(" ( ", "space").replace(" (", "space").split("space")
    fields.append(field)
# print(fields)

for i in range(1, len(f), 2):
    dG_refs.append(fields[i][1].replace(")", ""))
    dot_bracket_structures.append(fields[i][0])
# print(dG_refs)
# print(dot_bracket_structures)


for i in range(0, len(f), 2):
    sequences.append(fields[i][0])
# print(sequences)
# print(len(sequences))
#---------------------------dGconst----------------------------------------------------
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

#---------------------------dGunfold----------------------------------------------------
dG_unfolds = []
for i in range(0, len(dG_refs), 1):
    dG_unfolds.append(float(dG_consts[i])-float(dG_refs[i]))
# for dG_unfold in dG_unfolds:
#     print(round(dG_unfold, 2))


#----------------------contraints----------------------------
with open("constraints_file.pkl", "rb") as constraints_file:
    constraints = pickle.load(constraints_file)
# print(constraints)

#------------------------------final printing------------------------------------------------------------
final = [[] for i in range(len(dG_unfolds))]
for i in range(len(dG_unfolds)):
    final[i].append(sequences[i])
    final[i].append(dot_bracket_structures[i])
    final[i].append(constraints[i])
    final[i].append(round(float(dG_refs[i]), 2))
    final[i].append(round(float(dG_consts[i]), 2))
    final[i].append(round(float(dG_unfolds[i]), 2))

final_2 = sorted(final, key=itemgetter(5))
# print(final_2)

for i in range(len(dG_unfolds)):
    # print(sequences[i], dot_bracket_structures[i], constraints[i], round(float(dG_refs[i]), 2), round(float(dG_consts[i]), 2), round(float(dG_unfolds[i]), 2), sep="\t")
    print(*final_2[i], sep="\t")






















#--------------------------------end------------------------------------------------------------------------
