#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#----------------------------libraries----------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

import pickle
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

with open("start_positions_file.pkl", "rb") as start_positions_file:
    start_positions = pickle.load(start_positions_file)

# print(start_positions)
for i in range(len(start_positions)):
    print(start_positions[i])