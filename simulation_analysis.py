########################libraries #############################################
import numpy as np
import math
import pandas as pd
import random
import matplotlib.pyplot as plt
#import seaborn as sns
import random
import sys
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import os
from my_pipeline import *



####################### params ####################################################
params["main_output_folder"] = "/mnt/c/Users/Elinor/Desktop/תואר שני/simultation outputs/"  #outoutfolder

#####################functions###################################################
def create_folder(main_folder_path,new_folder_name):
    """creates new folder if it doesn't already exists"""
    dirName =new_folder_name
    if not os.path.exists(main_folder_path+dirName):
        os.chdir(main_folder_path) # enable to add sub folder to the main folder
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ")
        return  main_folder_path + new_folder_name+"/"
    else:
        print("Directory " , dirName ,  " already exists")

def graph_for_intresting_col(df,intresting_col):
    """gets the relevant columns and rturns their graphes and the difference between MCMC and all the data"""
    fig, (ax1, ax2) = plt.subplots(2,figsize=(20, 10))

    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=1)

    # Plot  each graph, and manually set the y tick values
##    ax1.scatter([i for i in range(len(df.loc[df['probabilty_res_MCMC'] == "True"]))], df[str(intresting_col)][d#f##['probabilty_res_MCMC'] == "True"])
    plt.title("MCMC probabilities")
    ax1.scatter(df.loc[df['probabilty_res_MCMC'] == "True"].index, df[str(intresting_col)][df['probabilty_res_MCMC'] == "True"])

    ax1.title.set_text('{} number of binders MCMC'.format(str(intresting_col)))
    ax1.set(xlabel='Peptide ', ylabel='{} number of binders MCMC'.format(str(intresting_col)))

    ax2.scatter(df.loc[df['probabilty_res_MCMC'] == "False" ].index, df[str(intresting_col)][df['probabilty_res_MCMC'] =="False"])#ask Tomer !!!!
    plt.title('{} number of all data'.format(str(intresting_col)))
    ax2.set(xlabel='peptide ', ylabel="{} of all data".format(str(intresting_col)))


#################################################### Analysis ###########################

#
# df=simulation_process(my_peptide,"av_of_total_binders")
# df1=simulation_process(my_peptide,"min_rank")
# df2=simulation_process(my_peptide,"median")
#
# df3=simulation_process(my_peptide,"total_super_binders")
# df4=simulation_process(my_peptide,"sb_supertypes")
# df5=simulation_process(my_peptide,"wb_supertypes")
# df6=simulation_process(my_peptide,"nb_supertypes")
#df7=simulation_process(my_peptide,"average")
#df8=simulation_process(my_peptide,"sum_of_binders")
df9=simulation_process(my_peptide,"average_calculated_binders")
df_dict={}
#dict key is the expirment name and the value is the df
#df_dict={"av_of_total_binders":df,"min_rank":df1,"median":df2,"total_super_binders":df3,"sb_super_types":df4,"wb_supertypes":df5,"nb_supertypes":df6,"average":df7}

df_dict={"average":df9}

for expirement in df_dict.keys():

    # # first expirement calculating how the av_of_total_binders changing during MCMC
        expirment_name=expirement+" seed 42"
        new_path=create_folder(params["main_output_folder"],expirment_name)
        print(new_path)
        graph_for_intresting_col(df_dict[expirement],"delta")
        plt.savefig(str(new_path)+" "+expirment_name+" delta.png",dpi=100)
        plt.close()

    # c=delta_graph(i,"probabilty_res_MCMC","all_data_prob")
        # total_bindes_graph(i,"total_binders")
        graph_for_intresting_col(df_dict[expirement],"wb_supertypes")
        plt.savefig(str(new_path)+" "+expirment_name+" wb_supertypes.png",dpi=100)
        plt.close()

#graph_for_intresting_col(i,"wb_supertypes")
        graph_for_intresting_col(df_dict[expirement],"sb_supertypes")
        plt.savefig(str(new_path)+" "+expirment_name+" sb_supertypes.png",dpi=100)
        plt.close()

        graph_for_intresting_col(df_dict[expirement],"nb_supertypes")
        plt.savefig(str(new_path)+" "+expirment_name+" nb_supertypes.png",dpi=100)
        plt.close()

# i["total_supertypes"]=i["sb_supertypes"]+i["wb_supertypes"]
    # graph_for_intresting_col(i,"total_supertypes")


        graph_for_intresting_col(df_dict[expirement],"total_binders")
        plt.savefig(str(new_path)+" "+expirment_name+" total_super_binders.png",dpi=100)
        plt.close()









