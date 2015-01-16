import sys 
import os
import shutil
import Get_mature_feature
import Get_precursor_feature
import Get_miBase_feature
import Get_RF_feature
import Get_miBase_RF_feature
import Get_arff_file

def Get_accession_family(accession_family_file):#input the accession family file get the mature-family info 
    fr = open(accession_family_file,"r")
    accession_family={}
    for line in fr:
        s = line.strip().split("\t")
        region = s[0]
        RF = s[1]
        miBase = s[2]
        accession_family[region]=(RF,miBase)
    fr.close()
    return accession_family #return the dict of the mature family

def Get_set(Dir , new_file_name , file_name,accession_family):
    new_dir = new_file_name+"/"+file_name+"/"
    os.makedirs(new_dir)#create the new_Breast/Breast NT/ the new file will save in this dir 
    files = os.listdir(Dir) 
    mature_set=set()
    ID_mature_set=set()
    ID_set=set()
    miBase_set=set()
    RF_set=set()
    miBase_RF_set=set()
    ID_miBase_set=set()
    ID_RF_set=set()
    ID_miBase_RF_set=set()
    for file in files:
        fr = open(Dir+file,"r")
        fw = open(new_dir+file,"w+")
        print "now process the "+file
        fw.write("miRNA_ID\tisoform_coords\tread_count\tmiRNA_region\tmiBase-family\tRF-family\n")
        for line in fr:
            s = line.strip().split("\t")
            if s[5].find("star")==-1 and s[5].find("mature")==-1:
               continue
            mature=s[5].split(",")[1]
            if mature not in accession_family.keys():
                fw.write(s[0]+"\t"+s[1]+"\t"+s[2]+"\t"+s[5].split(",")[1]+"\t"+"No MIR_family\tNo RF_family\n")
                miBase="No MIR_family"
                RF ="No RF_family"
            else:
                miBase = accession_family[mature][1] # the dict of accession_family is: mature:(RF,miBase)
                RF = accession_family[mature][0]
            fw.write(s[0]+"\t"+s[1]+"\t"+s[2]+"\t"+s[5].split(",")[1]+"\t"+miBase+"\t"+RF+"\n") #write the field to the new file 
            ID = s[0] 
            ID_mature = ID+"_"+mature
            ID_miBase = ID+"_"+miBase
            ID_RF = ID+"_"+RF
            ID_miBase_RF = ID+"_"+miBase+"_"+RF
            miBase_RF = miBase+"_"+RF
            # get the whole sets
            if mature not in mature_set:
               mature_set.add(mature)
            if ID not in ID_set:
               ID_set.add(ID)
            if ID_mature not in ID_mature_set:
               ID_mature_set.add(ID_mature)
            if miBase not in miBase_set and miBase!="No MIR_family":
               miBase_set.add(miBase)
            if RF not in RF_set and RF!="No RF_family":
               RF_set.add(RF)
            if miBase_RF not in miBase_RF_set and (miBase!="No MIR_family" and RF!="No RF_family"):
               miBase_RF_set.add(miBase_RF)
            if ID_miBase not in ID_miBase_set and miBase!="No MIR_family":
                ID_miBase_set.add(ID_miBase)
            if ID_RF not in ID_RF_set and RF!="No RF_family":
                ID_RF_set.add(ID_RF)
            if ID_miBase_RF not in ID_miBase_RF_set and (miBase!="No MIR_family" and RF!="No RF_family"):
                ID_miBase_RF_set.add(ID_miBase_RF)
        fr.close()
        fw.close()
    return mature_set , ID_mature_set , ID_set , miBase_set , RF_set , miBase_RF_set , ID_miBase_set ,ID_RF_set ,ID_miBase_RF_set
    
def Get_whole_info(accession_family,file_name):#get the whole set of mature, ID_mature ,ID , miBase , RF , miBase_RF
    new_file_name = "new_"+file_name
    if os.path.exists(new_file_name):
       shutil.rmtree(new_file_name)
       os.makedirs(new_file_name)
    else:
       os.makedirs(new_file_name)
    NT_mature_set, NT_ID_mature_set, NT_ID_set, NT_miBase_set, NT_RF_set,\
    NT_miBase_RF_set ,NT_ID_miBase_set , NT_ID_RF_set , NT_ID_miBase_RF_set \
        = Get_set(file_name+"/"+file_name+" NT/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",\
                               new_file_name, file_name+" NT",accession_family)   #tget the NT set 
    TN_mature_set, TN_ID_mature_set, TN_ID_set, TN_miBase_set, TN_RF_set, \
    TN_miBase_RF_set,TN_ID_miBase_set , TN_ID_RF_set , TN_ID_miBase_RF_set \
        = Get_set(file_name+"/"+file_name+" TN/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",\
                               new_file_name, file_name+" TN",accession_family)
    whole_mature = NT_mature_set | TN_mature_set #combine the NT and TN set 
    whole_ID_mature = NT_ID_mature_set | TN_ID_mature_set
    whole_ID = NT_ID_set | TN_ID_set
    whole_miBase = NT_miBase_set | TN_miBase_set
    whole_RF = NT_RF_set | TN_RF_set
    whole_miBase_RF = NT_miBase_RF_set | TN_miBase_RF_set
    whole_ID_miBase = NT_ID_miBase_set | TN_ID_miBase_set
    whole_ID_RF = NT_ID_RF_set | TN_ID_RF_set
    whole_ID_miBase_RF = NT_ID_miBase_RF_set | TN_ID_miBase_RF_set
    print "the len of mature is :"+str(len(NT_mature_set))+" | "+str(len(TN_mature_set))+" = "+str(len(whole_mature))
    print "the len of ID_mature is :"+str(len(NT_ID_mature_set))+" | "+str(len(TN_ID_mature_set))+" = "+str(len(whole_ID_mature))
    print "the len of ID is :"+str(len(NT_ID_set))+" | "+str(len(TN_ID_set))+" = "+str(len(whole_ID))
    print "the len of miBase is :"+str(len(NT_miBase_set))+" | "+str(len(TN_miBase_set))+" = "+str(len(whole_miBase))
    print "the len of RF is :"+str(len(NT_RF_set))+" | "+str(len(TN_RF_set))+" = "+str(len(whole_RF))
    print "the len of miBase_RF is :"+str(len(NT_miBase_RF_set))+" | "+str(len(TN_miBase_RF_set))+" = "+str(len(whole_miBase_RF))
    return whole_mature , whole_ID_mature , whole_ID , whole_miBase , whole_RF , whole_miBase_RF ,whole_ID_miBase, whole_ID_RF,whole_ID_miBase_RF

def main():
   accession_family_file = "Accession_Family.txt"
   file_name = sys.argv[1]
   accession_family  = Get_accession_family(accession_family_file)#get the accession family info save in dict
   whole_mature , whole_ID_mature , whole_ID , whole_miBase , whole_RF ,\
   whole_miBase_RF, whole_ID_miBase, whole_ID_RF,whole_ID_miBase_RF\
       = Get_whole_info(accession_family,file_name) #get the all whole info and creat new file
   Get_mature_feature.get_mature_features("new_"+file_name+"/"+file_name+" TN",whole_mature, whole_ID_mature,whole_ID,whole_miBase,\
                whole_RF,whole_miBase_RF,file_name+"_TN")
   Get_mature_feature.get_mature_features("new_"+file_name+"/"+file_name+" NT",whole_mature, whole_ID_mature,whole_ID,whole_miBase,\
                whole_RF,whole_miBase_RF,file_name+"_NT")
   Get_precursor_feature.get_precursor_features("new_"+file_name+"/"+file_name+" TN",whole_mature, whole_ID_mature,whole_ID,whole_miBase,\
                whole_RF,whole_miBase_RF,file_name+"_TN")
   Get_precursor_feature.get_precursor_features("new_"+file_name+"/"+file_name+" NT",whole_mature, whole_ID_mature,whole_ID,whole_miBase,\
                whole_RF,whole_miBase_RF,file_name+"_NT")
   Get_miBase_feature.get_miBase_features("new_"+file_name+"/"+file_name+" NT",whole_mature, whole_ID_mature,whole_ID,whole_miBase,\
                whole_RF,whole_miBase_RF,whole_ID_miBase, whole_ID_RF,whole_ID_miBase,file_name+"_NT")
   Get_miBase_feature.get_miBase_features("new_"+file_name+"/"+file_name+" TN",whole_mature, whole_ID_mature,whole_ID,whole_miBase,\
                whole_RF,whole_miBase_RF,whole_ID_miBase, whole_ID_RF,whole_ID_miBase,file_name+"_TN")
   Get_RF_feature.get_RF_features("new_"+file_name+"/"+file_name+" NT",whole_mature, whole_ID_mature,whole_ID,whole_miBase,\
                whole_RF,whole_miBase_RF,whole_ID_miBase, whole_ID_RF,whole_ID_miBase,file_name+"_NT")
   Get_RF_feature.get_RF_features("new_"+file_name+"/"+file_name+" TN",whole_mature, whole_ID_mature,whole_ID,whole_miBase,\
                whole_RF,whole_miBase_RF,whole_ID_miBase, whole_ID_RF,whole_ID_miBase,file_name+"_TN")
   Get_miBase_RF_feature.get_miBase_RF_features("new_"+file_name+"/"+file_name+" NT",whole_mature, whole_ID_mature,whole_ID,whole_miBase,\
                whole_RF,whole_miBase_RF ,whole_ID_miBase_RF,file_name+"_NT")
   Get_miBase_RF_feature.get_miBase_RF_features("new_"+file_name+"/"+file_name+" TN",whole_mature, whole_ID_mature,whole_ID,whole_miBase,\
                whole_RF,whole_miBase_RF ,whole_ID_miBase_RF,file_name+"_TN")
   Get_arff_file.get_arff_file(file_name,whole_mature,whole_ID_mature, whole_ID , whole_miBase,whole_RF, whole_miBase_RF, \
                            whole_ID_miBase, whole_ID_RF , whole_ID_miBase_RF)
   shutil.rmtree("new_"+file_name)
   os.remove("mature_"+file_name+"_TN")
   os.remove("mature_"+file_name+"_NT")
   os.remove("precursor_"+file_name+"_TN")
   os.remove("precursor_"+file_name+"_NT")
   os.remove("miBase_"+file_name+"_TN")
   os.remove("miBase_"+file_name+"_NT")
   os.remove("RF_"+file_name+"_TN")
   os.remove("RF_"+file_name+"_NT")
   os.remove("miBase_RF_"+file_name+"_TN")
   os.remove("miBase_RF_"+file_name+"_NT")


if __name__=="__main__":
   main()
