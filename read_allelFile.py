f = open("/mnt/c/Users/Elinor/PycharmProjects/pythonProject1/Weiskopf_2013_et al_HLA_class_I_plus_freq.txt", "r")
my_hla=open("/mnt/c/Users/Elinor/PycharmProjects/pythonProject1/HLA supertype MHC1.txt","r")
hla_used_for_my_pred=[]
hla_frequency={}
for line in f:
    line=line.rstrip("\n")
    print(line[0:12])
    if line.startswith("HLA"):
        hla_frequency[line[0:11]]=float(line[12:])

print(hla_frequency)