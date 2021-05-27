#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import subprocess
import os
import random
import pandas as pd
import math
import numpy as np
import runpy








##################################################################
# pathes
####################################################################
# path_to_tool = "/home/elinorpe/netMHCpan-4.1/"
# path_to_save ='/mnt/c/Users/Elinor/PycharmProjects/pythonProject1/output/'
# output_file = path_to_save +"pred_name.txt"
#
# #output_dir = work_dir
#
# if  os.path.exists(path_to_save) and os.path.exists(output_file) and os.path.exists(path_to_tool) :
#     print('found all paths  :) ')
# else:
#     print('cant find a path  :( ')







#################################################### Params ###########################
params={}
params["probability_function"]=lambda x :1.0/(1+math.exp(x) +0.1) #arbitrary probability function
params["seed"]=random.seed(42)


################################################## functions #########################

def mutation_creator(peptide):
    amino_acid_list=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    index=range(9)
    index_number=random.choice(index)
    old_base=peptide[index_number] #the base was
    random_amin_acid=random.choice(amino_acid_list) #new base
    peptide ="".join((peptide[:index_number],random_amin_acid,peptide[index_number+1:]))
    return peptide




def send_pep_to_pred(peptide):
    all_peps_sent_to_prediction = []
    mutated_list=[]
    open('/mnt/c/Users/Elinor/PycharmProjects/pythonProject1/input/peptides_for_pred.txt', 'w').close() # clearing pred input file
    input_file='/mnt/c//Users/Elinor/PycharmProjects/pythonProject1/input/peptides_for_pred.txt'
    while len(mutated_list)<2:
        with open('/mnt/c//Users/Elinor/PycharmProjects/pythonProject1/input/peptides_for_pred.txt', 'a') as file_object:
            mutated_list.append(peptide)
            file_object.write(peptide+'\n')
            all_peps_sent_to_prediction.append(peptide)
            peptide=mutation_creator(peptide)
            file_object.close()
    path_to_tool = "/home/elinorpe/netMHCpan-4.1/"
    path_to_save = "/mnt/c/Users/Elinor/PycharmProjects/pythonProject1/output/"
    output_file = path_to_save + "pred_name.txt"
    HLA_str = 'HLA-A0101,HLA-A0201,HLA-A0203,HLA-A0206,HLA-A0301,HLA-A1101,HLA-A2301,HLA-A2402,HLA-A2601,HLA-A3001,HLA-A3002,HLA-A3101,HLA-A3201,HLA-A3301,HLA-A6801,HLA-A6802,HLA-B0702,HLA-B0801,HLA-B1501,HLA-B3501,HLA-B4001,HLA-B4403,HLA-B4402,HLA-B5101,HLA-B5301,HLA-B5701,HLA-B5801'
    HLA_first_str = HLA_str[0:999]
    HLA_second_str = HLA_str[1000:]
    comand1 = path_to_tool + "netMHCpan " + "-p " + input_file + " -l " + "9" + " -a " + HLA_first_str + " >" + output_file
    comand2 = path_to_tool + "netMHCpan " + "-p " + input_file + " -l " + "9" + " -a " + HLA_second_str + " >>" + output_file
    raw_output = subprocess.check_output('%s' % comand1, shell=True)
    print("mut list" +str(mutated_list))
    print("sent to prediction")
    print ('Script1 ended')
    print ('Starting script2 ,analyzing ..')
    import pred_analysis
    return pred_analysis.df_creator("/mnt/c/Users/Elinor/PycharmProjects/pythonProject1/output/pred_name.txt")





def simulation(df,function,col_contains_data):
    probabilty_res_MCMC = ["First"]
    all_data_prob = ["First"]
    Delta=["First"]
    """gets df,function, and the column that contains the data with the score im interested to check ,
    insert the data of each row into the column and return df and flag considering the simulation result"""
    delta=np.diff(df[col_contains_data]) #Calculating the delta between two values in the columns
    Delta.append(delta)
    prob_res=params["probability_function"](delta)
    params["seed"]
    random_toss=random.random()
    all_data_prob.append(prob_res)


    if random_toss<=prob_res:
        flag="True"
        probabilty_res_MCMC.append(flag) #true if the results are under the probability calculated
    else: #random_toss>prob_res :
        flag= "False"
        probabilty_res_MCMC.append(flag)

    df["probabilty_res_MCMC"]=probabilty_res_MCMC
    df["all_data_prob"] =all_data_prob
    df["delta"]=Delta
    return df, flag

#
my_peptide= "CDTINCERY"
def simulation_process(peptide,column):
    """gets a peptide and column to calculate the simulation function and return df with all the result """
    appended_data=[]
    while len (appended_data)<20 :# there i decide the stop condition for this process
            df1=pd.DataFrame()
            df1=send_pep_to_pred(peptide)
            if simulation(df1,params["probability_function"],column)[1]=="True":
                old_pep = df1.tail(1).index[0]
                appended_data.append(simulation(df1,params["probability_function"],column)[0])
                peptide=old_pep
            if simulation(df1, params["probability_function"], column) [1]== "False":
                new_pep=df1.head(1).index[0] #genearating new peptide
                appended_data.append(simulation(df1,params["probability_function"],column)[0])
                peptide=new_pep

    appended_data = pd.concat(appended_data)
    return appended_data


simulation_result=simulation_process(my_peptide,"average")
#simulation_result.to_csv("simulation_result_example.csv")

