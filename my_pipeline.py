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
#params["probability_function"]=lambda x :(x+0.15)*10/12
params["seed"]=random.seed(86)
params["seed_mutation"]=random.seed(86) #seed for mutation creator
params["number of peptide for pred"]=5000


################################################fixing changes pc and lynux#############
path_flag=0#0=pc, 1=linux
if path_flag==0:
    main_path = '/mnt/c//Users/Elinor/PycharmProjects/project_elinor'
    path_to_tool = "/home/elinorpe/netMHCpan-4.1/"
    params["main_output_folder"] = os.path.join("/mnt/c/Users/Elinor/Desktop/תואר שני/", "simultation outputs_no_sim_26","")

elif path_flag==1:
    main_path='/home/perr/Desktop/sim/project_elinor/'
    path_to_tool ="/home/perr/netMHCpan-4.1/"
    params["main_output_folder"] = os.path.join(main_path, "5000 pep sum sim july seventh", "")




################################################## functions #########################

def mutation_creator(peptide):
    amino_acid_list=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    index=range(9)
    params["seed_mutation"]
    index_number=random.choice(index)
    old_base=peptide[index_number] #the base was
    random_amin_acid=random.choice(amino_acid_list) #new base
    while old_base==random_amin_acid :#preventing zero movement
        random_amin_acid = random.choice(amino_acid_list)  # new base
    # #else:
    # continue
    peptide ="".join((peptide[:index_number],random_amin_acid,peptide[index_number+1:]))
    position=index_number+1 #becuase index start from zero in phyton
    return peptide,position,old_base,random_amin_acid

def send_to_prediction_as_is(peptides_file):
    """gets a file with a list of peptides and returns df of their prediction without mutating them"""
    input_file= peptides_file
    path_to_save = os.path.join(main_path,"output",'')
    output_file=os.path.join(path_to_save,"prediction.txt")
    HLA_str = 'HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B27:05,HLA-B39:01,' \
              'HLA-B40:01,HLA-B58:01,HLA-B15:01'
    command = path_to_tool + "netMHCpan " + "-p " + input_file + " -l " + "9" + " -a " + HLA_str + " >" + output_file
    raw_output = subprocess.check_output('{}'.format(command), shell=True)
    import pred_analysis
    full_df=pred_analysis.df_creator(output_file)
    return full_df

#epitope_df=send_to_prediction_as_is("/mnt/c//Users/Elinor/PycharmProjects/project_elinor/input/to_prediction.txt")


def send_pep_to_pred(peptide,desired_number_of_peeptide_for_pred=2):
    input_file= os.path.join(main_path,"input","peptides_for_pred.txt")

    #path_to_tool = "/home/sacharen/netMHCpan-4.1/"
    path_to_save = os.path.join(main_path,"output",'')
    #utput_file = path_to_save + "pred_name.txt"
    output_file=os.path.join(path_to_save,"pred_name.txt")
    mutated_list=[]
    old_AA=["no change"]
    new_AA=["no change"]
    index_changed=["no change"]
    #open('/mnt/c/Users/Elinor/PycharmProjects/pythonProject1/input/peptides_for_pred.txt', 'w').close() # clearing pred input file
    mutated_list.append(peptide)
    while len(mutated_list)!=desired_number_of_peeptide_for_pred:
        new_peptide, position, old_base, random_amin_acid =mutation_creator(peptide)
        mutated_list.append(new_peptide)
        old_AA.append(old_base)
        new_AA.append(random_amin_acid)
        index_changed.append(position)
        peptide=new_peptide
    HLA_str = 'HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B27:05,HLA-B39:01,' \
              'HLA-B40:01,HLA-B58:01,HLA-B15:01'
    #HLA_str = 'HLA-A0101,HLA-A0201,HLA-A0203,HLA-A0206,HLA-A0301,HLA-A1101,HLA-A2301,HLA-A2402,HLA-A2601,HLA-A3001,HLA-A3002,HLA-A3101,HLA-A3201,HLA-A3301,HLA-A6801,HLA-A6802,HLA-B0702,HLA-B0801,HLA-B1501,HLA-B3501,HLA-B4001,HLA-B4403,HLA-B4402,HLA-B5101,HLA-B5301,HLA-B5701,HLA-B5801'
    #HLA_first_str = HLA_str[0:999]
    #HLA_second_str = HLA_str[1000:]
    if len(mutated_list)==desired_number_of_peeptide_for_pred:
        with open(input_file,'w+') as file_object:
            for peptide in mutated_list:
                    file_object.write(peptide+'\n')
            command = path_to_tool + "netMHCpan " + "-p " + input_file + " -l " + "9" + " -a " + HLA_str + " >" + output_file
            file_object.truncate() # clearing pred input file
            file_object.close()
    #comand2 = path_to_tool + "netMHCpan " + "-p " + input_file + " -l " + "9" + " -a " + HLA_second_str + " >>" + output_file
    #raw_output = subprocess.check_output('%s' % comand1, shell=True)
    raw_output = subprocess.check_output('{}'.format(command), shell=True)
    print("mut list" +str(mutated_list))
    print("sent to prediction")
    print ('Script1 ended')
    print ('Starting script2 ,analyzing ..')

    import pred_analysis
    full_df=pred_analysis.df_creator(output_file)
    full_df["position_changed"]=index_changed
    full_df["former_AA"]=old_AA
    full_df["new_AA"]=new_AA
    return full_df



def simulation(df,function,col_contains_data):
    probabilty_res_MCMC = ["First"]
    all_data_prob = ["First"]
    Delta=["First"]
    """gets df,function, and the column that contains the data with the score im interested to check ,
    insert the data of each row into the column and return df and flag considering the simulation result"""
    delta=np.diff(df[col_contains_data]) #Calculating the delta between two values in the columns

    Delta.append(float(delta[0]))
    try:
        prob_res =params["probability_function"](float(delta[0]))

    except OverflowError:
        prob_res =0.00001
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

# df8=df8[0:2]
# df8.head(1).index[0]
# df8.set_index(["Peptide"], inplace=True)


def simulation_process(peptide,column):

    """gets a peptide and column to calculate the simulation function and return df with all the result """
    appended_data=pd.DataFrame()
    flag=False # this flag will turn to true when i reach to the desired amount of peptide that are good for my simulation
    while flag ==False:
        #df1=pd.DataFrame()
        df1=send_pep_to_pred(peptide)
        print(df1)
        print(params["probability_function"])
        the_simulation=simulation(df1,params["probability_function"],column) #saving the simulation result + flag
        print(the_simulation[0])

        if the_simulation[1]=="True":
            old_pep = the_simulation[0].tail(1).index[0]
            print("old pep"+old_pep)
            appended_data=appended_data.append(the_simulation[0])
            peptide=old_pep
        elif the_simulation[1]== "False":
            new_pep=the_simulation[0].head(1).index[0] #genearating new peptide
            print("new pep"+new_pep)
            appended_data=appended_data.append(the_simulation[0])
            peptide=new_pep

        #if len(appended_data)>=2:
        #appended_data .append(df1)
        print(appended_data)
        if len(appended_data.loc[appended_data['probabilty_res_MCMC'] == "True"])==4:  # there i decide the stop condition for this process
            flag=True
    return appended_data
#
# #
# #
res_2=simulation_process(my_peptide,"average")
df1 = send_pep_to_pred(my_peptide)
# df2= send_pep_to_pred(my_peptide)
# df1=df1.append(df2,ignore_index=True)
#
# df3=pd.DataFrame()
# df3=df3.append(df1,ignore_index=True)
###
# print(df1)