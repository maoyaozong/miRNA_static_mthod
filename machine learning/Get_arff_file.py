import sys

def get_arff_file(file_name,whole_mature,whole_ID_mature, whole_ID , whole_miBase,whole_RF, whole_miBase_RF, \
                            whole_ID_miBase, whole_ID_RF , whole_ID_miBase_RF):
#the feature of mature is feature_mature_number , feature_everyline_avg_in_block ,whole_mature-rates , whole_ID_mature-rates , whole_mature-IDs
#the feature of precursor is :feature_ID_number, feature_everyline_avg_in_block , whole_ID-rates ,whole_ID_mature-rates, whole_ID-matures
#the feature of miBase is :feature_family_number , feature_everyline_avg_in_block , none_miBase_family_rates , reads_sum_in_none_miBase
#                          whole_miBase-rates , whole_ID_miBase-rates , whole_miBase-IDs
#the feature of RF is :feature_family_number , feature_everyline_avg_in_block , none_RF_family_rates ,reads_sum_in_none_RF
#                          whole_RF-rates , whole_ID_RF-rates , whole_RF-IDs
#the feature of miBase_RF is : feature_family_number ,  feature_everyline_avg_in_block , none_miBase_RF_family_number , reads_sum_in_none_miBase_RF
#                           whole_miBase_RF , whole_ID_miBase_RF , whole_miBase_RF-IDs
       fw=open(file_name+".arff","w+") # this file_name comtain the TN or NT style

       # now start to create the arff file
       fw.write("@relation "+file_name+"\n")
       fw.write("@attribute mature_feature_mature_number real\n")#mature
       fw.write("@attribute mature_feature_everyline_avg_in_block real\n")
       mature_rates_list=[]  # get the mature-rates feature
       for item in whole_mature:
            mature_rates_list.append(item)
       mature_rates_list.sort()
       for item in mature_rates_list:
            fw.write("@attribute mature_"+item+"-rates real\n")
       ID_mature_rates_list=[] # get the ID_mature_rates feature
       for item in whole_ID_mature:
            ID_mature_rates_list.append(item)
       ID_mature_rates_list.sort()
       for item in ID_mature_rates_list:
            fw.write("@attribute mature_"+item+"-rates real\n")
       mature_IDs_list=[] # the the ID_mature_rates's avg
       for item in whole_mature:
            mature_IDs_list.append(item)
       mature_IDs_list.sort()
       for item in mature_IDs_list:
            fw.write("@attribute mature_"+item+"-avg real\n")
       #precursor
       fw.write("@attribute precursor_feature_precursor_number real\n")
       fw.write("@attribute precursor_feature_everyline_avg_in_block real\n")
       precursor_rates_list=[]  # get the mature-rates feature
       for item in whole_ID:
            precursor_rates_list.append(item)
       precursor_rates_list.sort()
       for item in precursor_rates_list:
            fw.write("@attribute precursor_"+item+"-rates real\n")
       ID_mature_rates_list=[] # get the ID_mature_rates feature
       for item in whole_ID_mature:
            ID_mature_rates_list.append(item)
       ID_mature_rates_list.sort()
       for item in ID_mature_rates_list:
            fw.write("@attribute precursor_"+item+"-rates real\n")
       ID_matures_list=[] # the the ID_mature_rates's avg
       for item in whole_ID:
            ID_matures_list.append(item)
       ID_matures_list.sort()
       for item in ID_matures_list:
            fw.write("@attribute precursor_"+item+"-avg real\n")
       #miBase
       fw.write("@attribute miBase_feature_number real\n")
       fw.write("@attribute miBase_feature_everyline_avg_in_block real\n")
       fw.write("@attribute mibase_none_miBase_family_rates real\n")
       fw.write("@attribute miBase_reads_sum_in_none_miBase real\n")
       miBase_rates_list=[]  # get the mature-rates feature
       for item in whole_miBase:
            miBase_rates_list.append(item)
       miBase_rates_list.sort()
       for item in miBase_rates_list:
            fw.write("@attribute miBase_"+item+"-rates real\n")
       ID_miBase_rates_list=[] # get the ID_mature_rates feature
       for item in whole_ID_miBase:
            ID_miBase_rates_list.append(item)
       ID_miBase_rates_list.sort()
       for item in ID_miBase_rates_list:
            fw.write("@attribute miBase_"+item+"-rates real\n")
       ID_miBase_list=[] # the the ID_mature_rates's avg
       for item in whole_miBase:
            ID_miBase_list.append(item)
       ID_miBase_list.sort()
       for item in ID_miBase_list:
            fw.write("@attribute miBase_"+item+"-avg real\n")
       #RF
       fw.write("@attribute RF_feature_number real\n")
       fw.write("@attribute RF_feature_everyline_avg_in_block real\n")
       fw.write("@attribute RF_none_miBase_family_rates real\n")
       fw.write("@attribute RF_reads_sum_in_none_miBase real\n")
       RF_rates_list=[]  # get the mature-rates feature
       for item in whole_RF:
            RF_rates_list.append(item)
       RF_rates_list.sort()
       for item in RF_rates_list:
            fw.write("@attribute RF_"+item+"-rates real\n")
       ID_RF_rates_list=[] # get the ID_mature_rates feature
       for item in whole_ID_RF:
            ID_RF_rates_list.append(item)
       ID_RF_rates_list.sort()
       for item in ID_RF_rates_list:
            fw.write("@attribute RF_"+item+"-rates real\n")
       ID_RF_list=[] # the the ID_mature_rates's avg
       for item in whole_RF:
            ID_RF_list.append(item)
       ID_RF_list.sort()
       for item in ID_RF_list:
            fw.write("@attribute RF_"+item+"-avg real\n")
       #miBase-RF
       fw.write("@attribute miBase_RF_feature_number real\n")
       fw.write("@attribute miBase_RF_feature_everyline_avg_in_block real\n")
       fw.write("@attribute miBase_RF_none_miBase_family_rates real\n")
       fw.write("@attribute miBase_RF_reads_sum_in_none_miBase real\n")
       miBase_RF_rates_list=[]  # get the mature-rates feature
       for item in whole_miBase_RF:
            miBase_RF_rates_list.append(item)
       miBase_RF_rates_list.sort()
       for item in miBase_RF_rates_list:
            fw.write("@attribute miBase_RF_"+item+"-rates real\n")
       ID_miBase_RF_rates_list=[] # get the ID_mature_rates feature
       for item in whole_ID_miBase_RF:
            ID_miBase_RF_rates_list.append(item)
       ID_miBase_RF_rates_list.sort()
       for item in ID_miBase_RF_rates_list:
            fw.write("@attribute miBase_RF_"+item+"-rates real\n")
       ID_miBase_RF_list=[] # the the ID_mature_rates's avg
       for item in whole_miBase_RF:
            ID_miBase_RF_list.append(item)
       ID_miBase_RF_list.sort()
       for item in ID_miBase_RF_list:
            fw.write("@attribute miBase_RF_"+item+"-avg real\n")
       fw.write("@attribute class {1,0}\n")
       fw.write("@data\n")

       fr_mature = open("mature_"+file_name+"_NT","r")
       fr_precursor = open("precursor_"+file_name+"_NT","r")
       fr_miBase = open("miBase_"+file_name+"_NT","r")
       fr_RF = open("RF_"+file_name+"_NT","r")
       fr_miBase_RF = open("miBase_RF_"+file_name+"_NT","r")
       mature_list=[]
       for line in fr_mature:
            mature_list.append(line.strip())
       fr_mature.close()
       precursor_list=[]
       for line in fr_precursor:
             precursor_list.append(line.strip())
       fr_precursor.close()
       miBase_list=[]
       for line in fr_miBase:
             miBase_list.append(line.strip())
       fr_miBase.close()
       RF_list=[]
       for line in fr_RF:
             RF_list.append(line.strip())
       fr_RF.close()
       miBase_RF_list=[]
       for line in fr_miBase_RF:
              miBase_RF_list.append(line.strip())
       fr_miBase_RF.close()
       sample_length = len(mature_list)
       if len(miBase_list)!=sample_length or len(precursor_list)!=sample_length or len(RF_list)!= sample_length\
           or len(miBase_RF_list)!=sample_length:
          print "the feature list is not equla length "
          exit()
       for i in range(sample_length):
           feature_line = mature_list[i]+"\t"+precursor_list[i]+"\t"+miBase_list[i]+"\t"+RF_list[i]+"\t"+miBase_RF_list[i]
           fw.write(",".join(feature_line.split("\t"))+","+str(1)+"\n")
       fr_mature = open("mature_"+file_name+"_TN","r")
       fr_precursor = open("precursor_"+file_name+"_TN","r")
       fr_miBase = open("miBase_"+file_name+"_TN","r")
       fr_RF = open("RF_"+file_name+"_TN","r")
       fr_miBase_RF = open("miBase_RF_"+file_name+"_TN","r")
       mature_list=[]
       for line in fr_mature:
            mature_list.append(line.strip())
       fr_mature.close()
       precursor_list=[]
       for line in fr_precursor:
             precursor_list.append(line.strip())
       fr_precursor.close()
       miBase_list=[]
       for line in fr_miBase:
             miBase_list.append(line.strip())
       fr_miBase.close()
       RF_list=[]
       for line in fr_RF:
             RF_list.append(line.strip())
       fr_RF.close()
       miBase_RF_list=[]
       for line in fr_miBase_RF:
              miBase_RF_list.append(line.strip())
       fr_miBase_RF.close()
       sample_length = len(mature_list)
       if len(miBase_list)!=sample_length or len(precursor_list)!=sample_length or len(RF_list)!= sample_length\
           or len(miBase_RF_list)!=sample_length:
          print "the feature list is not equla length "
          exit()
       for i in range(sample_length):
           feature_line = mature_list[i]+"\t"+precursor_list[i]+"\t"+miBase_list[i]+"\t"+RF_list[i]+"\t"+miBase_RF_list[i]
           fw.write(",".join(feature_line.split("\t"))+","+str(0)+"\n")
       fw.close()


