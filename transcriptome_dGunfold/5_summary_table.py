#####################################################################################################################################################################################################################################################################

                                       # libraries

#####################################################################################################################################################################################################################################################################
import pickle







#####################################################################################################################################################################################################################################################################

                                       # sys_argv

#####################################################################################################################################################################################################################################################################








#####################################################################################################################################################################################################################################################################

                                       # step 1 - opening all data

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

with open("true_file.pkl", "rb") as true_file:
    true = pickle.load(true_file)

with open("distances_file.pkl", "rb") as distances_file:
    distances = pickle.load(distances_file)

with open("start_codons_file.pkl", "rb") as start_codons_file:
    start_codons = pickle.load(start_codons_file)

with open("forms_file.pkl", "rb") as forms_file:
    forms = pickle.load(forms_file)


#####################################################################################################################################################################################################################################################################

                                       # step 2 - opening all new data

#####################################################################################################################################################################################################################################################################
with open("left_cordinates_new_file.pkl", "rb") as left_cordinates_new_file:
    left_cordinates_new = pickle.load(left_cordinates_new_file)

with open("right_cordinates_new_file.pkl", "rb") as right_cordinates_new_file:
    right_cordinates_new = pickle.load(right_cordinates_new_file)

with open("features_new_file.pkl", "rb") as features_new_file:
    features_new = pickle.load(features_new_file)

with open("locus_tags_new_file.pkl", "rb") as locus_tags_new_file:
    locus_tags_new = pickle.load(locus_tags_new_file)

with open("strands_new_file.pkl", "rb") as strands_new_file:
    strands_new = pickle.load(strands_new_file)

with open("true_new_file.pkl", "rb") as true_new_file:
    true_new = pickle.load(true_new_file)

with open("distances_new_file.pkl", "rb") as distances_new_file:
    distances_new = pickle.load(distances_new_file)

with open("start_codons_new_file.pkl", "rb") as start_codons_new_file:
    start_codons_new = pickle.load(start_codons_new_file)

with open("forms_new_file.pkl", "rb") as forms_new_file:
    forms_new = pickle.load(forms_new_file)








#####################################################################################################################################################################################################################################################################

                                       # step 3 - gene counts

#####################################################################################################################################################################################################################################################################
print("total genes: ", len(left_cordinates), sep ="\t")

#-----------------------------------------------------------------------------------------------------
total_CDSs_count = 0
for i in range(len(features)):
    if features[i][0] == "CDS":
        total_CDSs_count += 1
print("total CDSs: ", total_CDSs_count, sep ="\t")

#-----------------------------------------------------------------------------------------------------
mRNAs_with_tss = []
for i in range(len(distances_new)):
    if features_new[i] == "CDS" and locus_tags_new[i] not in mRNAs_with_tss:
        mRNAs_with_tss.append(locus_tags_new[i])
print("mRNAs with tss: ", len(mRNAs_with_tss), sep ="\t")

#-----------------------------------------------------------------------------------------------------
# genes_with_multiforms_count = 0
number_of_multiforms = []
for i in range(len(distances_new)):
    if forms_new[i] != "form_1" and features_new[i] == "CDS":# and locus_tags_new[i] not in number_of_multiforms:
        number_of_multiforms.append(locus_tags_new[i])
print("number of multiforms: ", len(number_of_multiforms), sep ="\t")

#-----------------------------------------------------------------------------------------------------
leaderless_count = 0
for i in range(len(distances_new)):
    if distances_new[i] == 0 and features_new[i] == "CDS":
        leaderless_count += 1
print("leaderless mRNAs: ", leaderless_count, sep ="\t")

#-----------------------------------------------------------------------------------------------------
leadered_count = 0
for i in range(len(distances_new)):
    if distances_new[i] > 0 and features_new[i] == "CDS":
        leadered_count += 1
print("leadered mRNAs: ", leadered_count, sep ="\t")

#-----------------------------------------------------------------------------------------------------
leadered = []

for i in range(len(distances_new)):
    if distances_new[i] > 0 and features_new[i] == "CDS" and locus_tags_new[i] not in leadered:
        leadered.append(locus_tags_new[i])

mRNAs_with_leadered_leaderless_forms = []
for i in range(len(distances_new)):
    if distances_new[i] == 0 and features_new[i] == "CDS" and locus_tags_new[i] in leadered:
        mRNAs_with_leadered_leaderless_forms.append(locus_tags_new[i])


print("mRNAs having leadered and leaderless forms: ", len(mRNAs_with_leadered_leaderless_forms), sep ="\t")
if len(mRNAs_with_leadered_leaderless_forms) == 0:
    pass
else:
    print("The following are genes that are transcribed as leadered as well as leaderless:")
    for i in range(len(mRNAs_with_leadered_leaderless_forms)):
        print(mRNAs_with_leadered_leaderless_forms[i])
#-----------------------------------------------------------------------------------------------------
ncRNAs_count = 0
for i in range(len(features)):
    if features[i][0] == "ncRNA":
        ncRNAs_count += 1
print("total ncRNAs: ", ncRNAs_count, sep ="\t")

#-----------------------------------------------------------------------------------------------------
tRNAs_count = 0
for i in range(len(features)):
    if features[i][0] == "tRNA":
        tRNAs_count += 1
print("total tRNAs: ", tRNAs_count, sep ="\t")

#-----------------------------------------------------------------------------------------------------
rRNAs_count = 0
for i in range(len(features)):
    if features[i][0] == "rRNA":
        rRNAs_count += 1
print("total rRNAs: ", rRNAs_count, sep ="\t")













#####################################################################################################################################################################################################################################################################

                                       # step 4 - start codon counts

#####################################################################################################################################################################################################################################################################
#----------------------------all-------------------------------------------------------------------------

ATG_count = 0
for i in range(len(distances_new)):
    if forms_new[i] == "form_1" and start_codons_new[i] == "ATG" and features_new[i] == "CDS":
        ATG_count += 1
print("all mRNAs starting with ATG: ", ATG_count, sep ="\t")

GTG_count = 0
for i in range(len(distances_new)):
    if forms_new[i] == "form_1" and start_codons_new[i] == "GTG" and features_new[i] == "CDS":
        GTG_count += 1
print("all mRNAs starting with GTG: ", GTG_count, sep ="\t")

CTG_count = 0
for i in range(len(distances_new)):
    if forms_new[i] == "form_1" and start_codons_new[i] == "CTG" and features_new[i] == "CDS":
        CTG_count += 1
print("all mRNAs starting with CTG: ", CTG_count, sep ="\t")

TTG_count = 0
for i in range(len(distances_new)):
    if forms_new[i] == "form_1" and start_codons_new[i] == "TTG" and features_new[i] == "CDS":
        TTG_count += 1
print("all mRNAs starting with TTG: ", TTG_count, sep ="\t")

ATA_count = 0
for i in range(len(distances_new)):
    if forms_new[i] == "form_1" and start_codons_new[i] == "ATA" and features_new[i] == "CDS":
        ATA_count += 1
print("all mRNAs starting with ATA: ", ATA_count, sep ="\t")

ATC_count = 0
for i in range(len(distances_new)):
    if forms_new[i] == "form_1" and start_codons_new[i] == "ATC" and features_new[i] == "CDS":
        ATC_count += 1
print("all mRNAs starting with ATC: ", ATC_count, sep ="\t")

ATT_count = 0
for i in range(len(distances_new)):
    if forms_new[i] == "form_1" and start_codons_new[i] == "ATT" and features_new[i] == "CDS":
        ATT_count += 1
print("all mRNAs starting with ATT: ", ATT_count, sep ="\t")

#------------------------------leaderless-----------------------------------------------------------------------
ATG_count_leaderless = 0
for i in range(len(distances_new)):
    if distances_new[i] == 0 and forms_new[i] == "form_1" and start_codons_new[i] == "ATG" and features_new[i] == "CDS":
        ATG_count_leaderless += 1
print("leaderless mRNAs starting with ATG: ", ATG_count_leaderless, sep ="\t")

GTG_count_leaderless = 0
for i in range(len(distances_new)):
    if distances_new[i] == 0 and forms_new[i] == "form_1" and start_codons_new[i] == "GTG" and features_new[i] == "CDS":
        GTG_count_leaderless += 1
print("leaderless mRNAs starting with GTG: ", GTG_count_leaderless, sep ="\t")

CTG_count_leaderless = 0
for i in range(len(distances_new)):
    if distances_new[i] == 0 and forms_new[i] == "form_1" and start_codons_new[i] == "CTG" and features_new[i] == "CDS":
        CTG_count_leaderless += 1
print("leaderless mRNAs starting with CTG: ", CTG_count_leaderless, sep ="\t")

TTG_count_leaderless = 0
for i in range(len(distances_new)):
    if distances_new[i] == 0 and forms_new[i] == "form_1" and start_codons_new[i] == "TTG" and features_new[i] == "CDS":
        TTG_count_leaderless += 1
print("leaderless mRNAs starting with TTG: ", TTG_count_leaderless, sep ="\t")

ATA_count_leaderless = 0
for i in range(len(distances_new)):
    if distances_new[i] == 0 and forms_new[i] == "form_1" and start_codons_new[i] == "ATA" and features_new[i] == "CDS":
        ATA_count_leaderless += 1
print("leaderless mRNAs starting with ATA: ", ATA_count_leaderless, sep ="\t")

ATC_count_leaderless = 0
for i in range(len(distances_new)):
    if distances_new[i] == 0 and forms_new[i] == "form_1" and start_codons_new[i] == "ATC" and features_new[i] == "CDS":
        ATC_count_leaderless += 1
print("leaderless mRNAs starting with ATC: ", ATC_count_leaderless, sep ="\t")

ATT_count_leaderless = 0
for i in range(len(distances_new)):
    if distances_new[i] == 0 and forms_new[i] == "form_1" and start_codons_new[i] == "ATT" and features_new[i] == "CDS":
        ATT_count_leaderless += 1
print("leaderless mRNAs starting with ATT: ", ATT_count_leaderless, sep ="\t")

#---------------------------------leadered--------------------------------------------------------------------
ATG_count_leadered = 0
for i in range(len(distances_new)):
    if distances_new[i] > 0 and forms_new[i] == "form_1" and start_codons_new[i] == "ATG" and features_new[i] == "CDS":
        ATG_count_leadered += 1
print("leadered mRNAs starting with ATG: ", ATG_count_leadered, sep ="\t")

GTG_count_leadered = 0
for i in range(len(distances_new)):
    if distances_new[i] > 0 and forms_new[i] == "form_1" and start_codons_new[i] == "GTG" and features_new[i] == "CDS":
        GTG_count_leadered += 1
print("leadered mRNAs starting with GTG: ", GTG_count_leadered, sep ="\t")

CTG_count_leadered = 0
for i in range(len(distances_new)):
    if distances_new[i] > 0 and forms_new[i] == "form_1" and start_codons_new[i] == "CTG" and features_new[i] == "CDS":
        CTG_count_leadered += 1
print("leadered mRNAs starting with CTG: ", CTG_count_leadered, sep ="\t")

TTG_count_leadered = 0
for i in range(len(distances_new)):
    if distances_new[i] > 0 and forms_new[i] == "form_1" and start_codons_new[i] == "TTG" and features_new[i] == "CDS":
        TTG_count_leadered += 1
print("leadered mRNAs starting with TTG: ", TTG_count_leadered, sep ="\t")

ATA_count_leadered = 0
for i in range(len(distances_new)):
    if distances_new[i] > 0 and forms_new[i] == "form_1" and start_codons_new[i] == "ATA" and features_new[i] == "CDS":
        ATA_count_leadered += 1
print("leadered mRNAs starting with ATA: ", ATA_count_leadered, sep ="\t")

ATC_count_leadered = 0
for i in range(len(distances_new)):
    if distances_new[i] > 0 and forms_new[i] == "form_1" and start_codons_new[i] == "ATC" and features_new[i] == "CDS":
        ATC_count_leadered += 1
print("leadered mRNAs starting with ATC: ", ATC_count_leadered, sep ="\t")

ATT_count_leadered = 0
for i in range(len(distances_new)):
    if distances_new[i] > 0 and forms_new[i] == "form_1" and start_codons_new[i] == "ATT" and features_new[i] == "CDS":
        ATT_count_leadered += 1
print("leadered mRNAs starting with ATT: ", ATT_count_leadered, sep ="\t")


#####################################################################################################################################################################################################################################################################

                                       # step 5 - diff 5'UTR length counts

#####################################################################################################################################################################################################################################################################








#####################################################################################################################################################################################################################################################################

                                       # end

#####################################################################################################################################################################################################################################################################
