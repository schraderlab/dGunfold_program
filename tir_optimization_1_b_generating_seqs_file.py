#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#----------------------------libraries----------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

import pickle
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

with open("TIRs_file.pkl", "rb") as TIRs_file:
    TIRs = pickle.load(TIRs_file)

# print(TIRs)
for i in range(len(TIRs)):
    print(TIRs[i])