#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
########################libraries ################
import numpy as np
import math
import pandas as pd
# import random
# import matplotlib.pyplot as plt
# import random
# import sys
# import matplotlib.pyplot as plt
# from matplotlib.collections import LineCollection
# from matplotlib.colors import ListedColormap, BoundaryNorm
import os



def df_creator(path):
    """get netMHCpan output and create df """
    with open(path, "r") as f:
        lines = f.readlines()
    with open("/mnt/c/Users/Elinor/PycharmProjects/pythonProject1/output/pred_without_unwanted_lines.txt", "w") as l:
        for line in lines:
                if line.startswith("#"):
                    continue
                if line.startswith("-"):
                    continue
                if line.startswith("P"):
                    continue
                if line.startswith("HLA-"):
                     continue
                if line.endswith("<= WB\n"):
                    line=line.replace('<= WB\n', '\n')
                if line.endswith("<= SB\n"):
                    line=line.replace('<= SB\n', '\n')
                l.write(line)

    #reading whitespace delimiter file
    mutant_df=pd.read_csv("/mnt/c/Users/Elinor/PycharmProjects/pythonProject1/output/pred_without_unwanted_lines.txt",delim_whitespace=True,skip_blank_lines=True,error_bad_lines=False,warn_bad_lines=False)

    mutant_df=mutant_df[~mutant_df['Peptide'].str.startswith('Peptide')]#removing redundat values
    mutant_df.set_index(["Peptide"],inplace=True)


    ##create df with hla information as columns
    list_of_hla=[]
    for hla in mutant_df["MHC"]:
        list_of_hla.append(hla)
    list_of_hla=list(set(list_of_hla))


    data_hla_as_col=pd.DataFrame()
    for hla in list_of_hla:
        data_hla_as_col[hla]=mutant_df["%Rank_EL"][mutant_df["MHC"]==hla]

    for hla in list_of_hla: #making sre that all the data is float type
        data_hla_as_col[hla].apply(lambda x: float(x))

    ## adding column of interst
    data_hla_as_col["average"]= data_hla_as_col.loc[:,list_of_hla].mean(axis=1).round()
    data_hla_as_col["min_rank"]= data_hla_as_col.loc[:,list_of_hla].min(axis=1)
    data_hla_as_col["median"]=data_hla_as_col.loc[:,list_of_hla].median(axis=1)
    data_hla_as_col = data_hla_as_col.astype(float)
    data_hla_as_col["av_of_total_binders"]=data_hla_as_col[data_hla_as_col.loc[:,list_of_hla]<float(2)].mean(axis=1)


    #
    # data_hla_as_col.head()
    #
    #

    data_hla_as_col=data_hla_as_col.round(decimals=2)
    #
    data_hla_as_col.fillna(0,inplace=True)
    #
    #
    # # %%
    def binding_feedback_func(x):
        """get values and returns the feedback regarding the binding"""
        if x <= 0.5:
            return 'SB'
        elif 0.5<x<=2:
            return 'WB'
        else:
            return 'NB'


    #
    #
    # # %%
    # #adding the cumber of each peptide is sb,nb,wb as columns
    nb=[]
    wb=[]
    sb=[]
    for ind, row in data_hla_as_col.iterrows():
        tmp=row.loc[list_of_hla].apply(binding_feedback_func).str
        nb.append(tmp.contains('NB', regex=False).sum())
        sb.append(tmp.contains('SB', regex=False).sum())
        wb.append(tmp.contains('WB', regex=False).sum())

    data_hla_as_col["NB"]=nb
    data_hla_as_col["SB"]=sb
    data_hla_as_col["WB"]=wb


    supertype_classification={"A01":["HLA-A*0101","HLA-A*2601","HLA-A*3002","HLA-A*3201"],"A02":["HLA-A*0201","HLA-A*0203" ,"HLA-A*0206"],"A03":["HLA-A*1101" ,"HLA-A*3101","HLA-A*0301","HLA-A*3301"],"B07":["HLA-B*0702","HLA-B*3501","HLA-B*5101","HLA-B*5301"],"B08": ["HLA-B*0801"],"B44":["HLA-B*4403", "HLA-B*4402" ,"HLA-B*4001"],"B58":["HLA-B*5801","HLA-B*5701"],"B62":["HLA-B*1501"],"A24":["HLA-A*2301" ,"HLA-A*2402"],"A01A03":[ "HLA-A*6801","HLA-A*6802"],"A03A02":["HLA-A*3001"]}

    # #creating a nested dictionary to look like {all_peptid:{supertype}:binder_classification:all_allels}
    #
    all_peps={}
    all_clasifications={}

    for value in data_hla_as_col.columns:
        data_hla_as_col[value].apply(lambda x: float(x))
    for ind,peptide in  data_hla_as_col.iterrows(): # for each row, each peptide , create dictionary and create count how many times each peptide bind to the supertypes
        pep_bindrs_clasification={}
        for superbinder in supertype_classification.keys(): #iterating over supertypes key name
            id_sb=[]
            id_wb=[]
            id_nb=[]
            binder_classifier={}
            for hla in supertype_classification[superbinder]:
                value=peptide.loc[hla]
                if value<=0.5:
                    id_sb.append(hla)
                    binder_classifier["sb"]=id_sb
    #                pep_dict[sb]=id_sb
                elif 0.5<value<=2:
                    id_wb.append(hla)
                    binder_classifier["wb"]=id_wb

                else :
                    id_nb.append(hla)
                    binder_classifier["nb"]=id_nb

            pep_bindrs_clasification[superbinder]=binder_classifier
        all_peps[ind]=pep_bindrs_clasification


    # # %%
    # #adding counter columns of each binder to the df
    sb_supertypes=[]
    wb_supertypes=[]
    nb_supertypes=[]
    sb=[]
    wb=[]
    for pep in all_peps.keys():
        counter_sb=0
        counter_wb=0
        nb_counter=11
        sb_identity=[]
        wb_identity=[]
        for super_type in supertype_classification.keys():
            for classification in all_peps[pep][super_type]:
                if classification=="sb" and all_peps[pep][super_type]["sb"]!=0: #if it is a strong binder
                    counter_sb+=1
                    sb_identity.append(super_type)

                if classification=="wb" and all_peps[pep][super_type]["wb"]!=0:
                    counter_wb+=1
                    wb_identity.append(super_type)


        sb_supertypes.append(counter_sb)
        wb_supertypes.append(counter_wb)
        nb_supertypes.append(nb_counter-counter_sb-counter_wb)
        sb.append(sb_identity)
        wb.append(wb_identity)


    if data_hla_as_col.iloc[0].name==data_hla_as_col.iloc[1].name: # if both peptides are the same
        data_hla_as_col["total_super_binders"] = (sb_supertypes[0] +wb_supertypes[0])*len(data_hla_as_col)
        data_hla_as_col["sb_supertypes"]=sb_supertypes*len(data_hla_as_col)
        print(sb_supertypes)
        print(data_hla_as_col)
        data_hla_as_col["wb_supertypes"]=wb_supertypes*len(data_hla_as_col)
        data_hla_as_col["nb_supertypes"]=nb_supertypes*len(data_hla_as_col)
        data_hla_as_col["sb_super_type_id"] = sb*len(data_hla_as_col)
        data_hla_as_col["wb_super_type_id"] = wb*len(data_hla_as_col)
        data_hla_as_col["total_binders"]= (data_hla_as_col["SB"] + data_hla_as_col["WB"])*len(data_hla_as_col)

    else:

        data_hla_as_col["nb_supertypes"]=nb_supertypes
        data_hla_as_col["sb_supertypes"]=sb_supertypes
        data_hla_as_col["wb_supertypes"]=wb_supertypes
        data_hla_as_col["sb_super_type_id"]=sb
        data_hla_as_col["wb_super_type_id"]=wb
        data_hla_as_col["total_binders"]=data_hla_as_col["SB"]+data_hla_as_col["WB"]
        data_hla_as_col["total_super_binders"]=data_hla_as_col["sb_supertypes"]+data_hla_as_col["wb_supertypes"]


    return data_hla_as_col


