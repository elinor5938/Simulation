########################libraries ##############################################
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
params={}
from my_pipeline import *

####################### params ####################################################
params["main_output_folder"]

#####################functions###################################################
def create_folder(main_folder_path,new_folder_name):
    """creates new folder if it doesn't already exists"""
    dirName =new_folder_name
    if not os.path.exists(main_folder_path+dirName):
        os.chdir(main_folder_path) # enable to add sub folder to the main folder
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ")
        return  main_folder_path + new_folder_name+"/"
    # else:
    #     print("Directory " , dirName ,  " already exists")

def graph_for_intresting_col(df,intresting_col):
    """gets the relevant columns and rturns their graphes and the difference between MCMC and all the data"""
    fig, (ax1, ax2) = plt.subplots(2,figsize=(20, 10))

    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=1)

    # Plot  each graph, and manually set the y tick values
##    ax1.scatter([i for i in range(len(df.loc[df['probabilty_res_MCMC'] == "True"]))], df[str(intresting_col)][d#f##['probabilty_res_MCMC'] == "True"])
    plt.title("MCMC probabilities")
    ax1.plot(df.loc[df['probabilty_res_MCMC'] == "True"].index, df[str(intresting_col)][df['probabilty_res_MCMC'] == "True"])

    ax1.title.set_text('{} change in MCMC peptides that got True flag'.format(str(intresting_col)))
    ax1.set(xlabel='Peptide ', ylabel=str(intresting_col))

    ax2.plot(df.loc[df['probabilty_res_MCMC'] == "False" ].index, df[str(intresting_col)][df['probabilty_res_MCMC'] =="False"])
    plt.title('{} change in MCMC peptides that got False flag'.format(str(intresting_col)))
    ax2.set(xlabel='peptide ', ylabel=(intresting_col))


def create_graph(df, x, y):
    """gets the relevant columns and rturns their graphes and the difference between MCMC and all the data"""
    fig, (ax1, ax2) = plt.subplots(2, figsize=(20, 10))

    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=1)

    # Plot  each graph, and manually set the y tick values
    ##    ax1.scatter([i for i in range(len(df.loc[df['probabilty_res_MCMC'] == "True"]))], df[str(intresting_col)][d#f##['probabilty_res_MCMC'] == "True"])
    plt.title("MCMC probabilities")

    ax1.plot(df[x][df['probabilty_res_MCMC'] == "True"].astype(float), df[y][df['probabilty_res_MCMC'] == "True"].astype(float))

    plt.title("association between {} to  {} in MCMC peptides that got True flag".format(x, y))
    ax1.set(xlabel=x, ylabel=y)

    ax2.plot(df[x][df['probabilty_res_MCMC'] == "False"].astype(float), df[y][df['probabilty_res_MCMC'] == "False"].astype(float))
    plt.title("association between {} to  {} in MCMC peptides that got False flag".format(x, y))
    ax2.set(xlabel=x, ylabel=y)


#################################################### Analysis ###########################

#
df8=simulation_process(my_peptide,"sum_of_all_hla")
df9=simulation_process(my_peptide,"average")
df10=simulation_process(my_peptide,"median")

#########################checking why it didn't work on linux#####################################

#dict key is the expirment name and the value is the df
#df_dict={"av_of_total_binders":df,"min_rank":df1,"median":df2,"total_super_binders":df3,"sb_super_types":df4,"wb_supertypes":df5,"nb_supertypes":df6,"average":df7}

df_dict={"sum_of_all_hla":df8,"average":df9,"median":df10}
y_parameters=["delta","WB","SB","NB","total_binders","all_data_prob"]

def simulation_analysis (df_dict,y_parameters):
    """ gets dictionary that it's keys are te parameters of experiment (the value that by it the delta function
    is calculated) and the value it's their dfs.and sanity checks params that will be the y axis
     and return graphs with sanity checks to examine my simulation"""


    for experiment in  df_dict.keys():
        expirment_name = " the col check is " + experiment + " the seed  is 86" # creationg a new folder which will contain the relevant graphs for every experiment
        new_path = create_folder(params["main_output_folder"], expirment_name)
        create_graph(df_dict[experiment],"delta","all_data_prob")
        plt.savefig(str(new_path) + " " + "delta association to probability" , dpi=100)
        create_graph(df_dict[experiment],"delta","all_data_prob")

        plt.close()


        for parameter  in y_parameters:

        # # first expirement calculating how the av_of_total_binders changing during MCMC
           graph_for_intresting_col(df_dict[experiment],parameter)

           plt.savefig(str(new_path)+" "+expirment_name+" "+parameter,dpi=100)
           plt.close()

simulation_analysis(df_dict,y_parameters)
for name in df_dict.keys(): #saving csvs
    df_dict[name].to_csv(name+" seed is 86.csv")
#