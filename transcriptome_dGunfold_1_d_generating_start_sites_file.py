#####################################################################################################################################################################################################################################################################

                                       # libraries

#####################################################################################################################################################################################################################################################################
import pickle









#####################################################################################################################################################################################################################################################################

                                       # data

#####################################################################################################################################################################################################################################################################
with open("start_positions_file.pkl", "rb") as start_positions_file:
    start_positions = pickle.load(start_positions_file)

# print(start_positions)
for i in range(len(start_positions)):
    for j in range(len(start_positions[i])):
        if(start_positions[i][j]) != "x":# and len(fifty_bases[i][j]) == 0:
            # print(i, j)
            print(start_positions[i][j])




#####################################################################################################################################################################################################################################################################

                                       # end

#####################################################################################################################################################################################################################################################################