import sys 
import os 
import math
import shutil

def Get_P(ztest):
    x=abs(ztest)/math.sqrt(2)
    T=(0.0705230784,0.0422820123,0.0092705272,0.0001520143,0.0002765672,0.0000430638)
    E=1-pow((1+sum([a*pow(x,(i+1))for i,a in enumerate(T)])),-16)
    if ztest:
       p=0.5-0.5*E
    else:
       p=0.5+0.5*E
    return p

def Process_dict_isomiR(dict1,dict2,file_name):
#the dict is :(ID+iso , (avg,var,read_counts))
#the output is : "miRNA_ID\tisoform_coords\tz-test\tp-value\tread_counts"
    fw = open("Z_P_out/"+file_name+"_Z_P.txt","w+")
    fw.write("miRNA_ID\tisoform_coords\tZ-test\tP-value\tread_counts\tTN-avg\tNT-avg\tTN-sta\tNT-sta\n")
    for key in dict1.keys():
        if key not in dict2.keys():
           fw.write(key+"\t"+"No z-test"+"\t"+"No p-value"+"\t"+str(dict1[key][3])+"\t"+str(dict1[key][0])+"\tNo_NT\t"+str(dict1[key][2])+"\tNo_NT\n")
           continue 
        else: # the key in the dict1 and the dict2 
           avg1=float(dict1[key][0])
           avg2=float(dict2[key][0])
           s1=float(dict1[key][1])
           s2=float(dict2[key][1])
#as we use the Sandard Devection so need not add the sqrt
#s1=math.sqrt(float(dict1[key][1]))
#s2=math.sqrt(float(dict2[key][1]))
           n1=float(dict1[key][2])
           n2=float(dict2[key][2])
           read_counts1=float(dict1[key][3])
           read_counts2=float(dict2[key][3])
           read_counts=read_counts1+read_counts2
           fenmu=math.sqrt(s1/n1+s2/n2)
           if fenmu == 0:
              Z_test = abs(avg1-avg2)
           else:
              Z_test = abs(avg1-avg2)/math.sqrt(s1/n1+s2/n2)
           P_value = Get_P(Z_test)
           fw.write(key+"\t"+str(Z_test)+"\t"+str(P_value)+"\t"+str(read_counts)+"\t"+str(avg1)+"\t"+str(avg2)+"\t"+str(s1)+"\t"+str(s2)+"\n")
    for key in dict2.keys():
        if key not in dict1.keys():
           fw.write(key+"\t"+"No z-test"+"\t"+"No p-value"+"\t"+str(dict2[key][3])+"\tNo_TN\t"+str(dict2[key][0])+"\tNo_TN\t"+str(dict2[key][1])+"\n")
    fw.close()  
'''
def Process_dict_mature(dict1,dict2,file_name):
    fw = open("Z_P_out/"+file_name+"_Z_P.txt","w+")
    fw.write("miRNA_ID\tmiRNA_region\tZ-test\tP-value\tread_counts\taccession\tproduct\tRF-family\tmiRNA-family\n")
    for key in dict1.keys():
        if key not in dict2.keys():
           fw.write(key+"\t"+"No z-test"+"\t"+"No p-value"+"\t"+str(dict1[key][3])+"\n")
           continue 
        else:
           avg1=float(dict1[key][0])
           avg2=float(dict2[key][0])
           s1=math.sqrt(float(dict1[key][1]))
           s2=math.sqrt(float(dict2[key][1]))
           n1=float(dict1[key][2])
           n2=float(dict2[key][2])
           fenmu=math.sqrt(s1/n1+s2/n2)
           if fenmu == 0:
              Z_test = abs(avg1-avg2)
           else:
              Z_test = abs(avg1-avg2)/math.sqrt(s1/n1+s2/n2)
           P_value = Get_P(Z_test)
           fw.write(key+"\t"+str(Z_test)+"\t"+str(P_value)+"\t"+str(dict1[key][3])+"\n")
    for key in dict2.keys():
        if key not in dict1.keys():
           fw.write(key+"\t"+"No z-test"+"\t"+"No p-value"+"\t"+str(dict2[key][3])+"\n")
    fw.close()  
'''

#output is : "miRNA_ID\tZ-test\tP-value\tread_counts\n"
def Process_dict_precursor(dict1,dict2,file_name):
    fw = open("Z_P_out/"+file_name+"_Z_P.txt","w+")
    fw.write("miRNA_ID\tisoform_coords\tZ-test\tP-value\tread_counts\tTN-avg\tNT-avg\tTN-sta\tNT-sta\n")
    for key in dict1.keys():
        if key not in dict2.keys():
           fw.write(key+"\t"+"No z-test"+"\t"+"No p-value"+"\t"+str(dict1[key][3])+"\t"+str(dict1[key][0])+"\tNo_NT\t"+str(dict1[key][2])+"\tNo_NT\n")
           continue
        else:
           avg1=float(dict1[key][0])
           avg2=float(dict2[key][0])
           s1=float(dict1[key][1])
           s2=float(dict2[key][1])
    #s1=math.sqrt(float(dict1[key][1]))
    #s2=math.sqrt(float(dict2[key][1]))
           n1=float(dict1[key][2])
           n2=float(dict2[key][2])
           read_counts1=float(dict1[key][3])
           read_counts2=float(dict2[key][3])
           read_counts=read_counts1+read_counts2
           fenmu=math.sqrt(s1/n1+s2/n2)
           if fenmu == 0:
              Z_test = abs(avg1-avg2)
           else:
              Z_test = abs(avg1-avg2)/math.sqrt(s1/n1+s2/n2)
           P_value = Get_P(Z_test)
           fw.write(key+"\t"+str(Z_test)+"\t"+str(P_value)+"\t"+str(read_counts)+"\t"+str(avg1)+"\t"+str(avg2)+"\t"+str(s1)+"\t"+str(s2)+"\n")
    for key in dict2.keys():
        if key not in dict1.keys():
           fw.write(key+"\t"+"No z-test"+"\t"+"No p-value"+"\t"+str(dict2[key][3])+"\tNo_TN\t"+str(dict2[key][0])+"\tNo_TN\t"+str(dict2[key][1])+"\n")
    fw.close()  

	
def Process_dict_mature_2(dict1,dict2,file_name):
#the output is : "miRNA_ID\tmiRNA_region\tz-test\tp-value\tread_counts"
    fw = open("Z_P_out/"+file_name+"_Z_P.txt","w+")
    fw.write("miRNA_ID\tisoform_coords\tZ-test\tP-value\tread_counts\tTN-avg\tNT-avg\tTN-sta\tNT-sta\n")
    for key in dict1.keys():
        if key not in dict2.keys():
           fw.write(key+"\t"+"No z-test"+"\t"+"No p-value"+"\t"+str(dict1[key][3])+"\t"+str(dict1[key][0])+"\tNo_NT\t"+str(dict1[key][2])+"\tNo_NT\n")
           continue 
        else:
           avg1=float(dict1[key][0])
           avg2=float(dict2[key][0])
           s1=float(dict1[key][1])
           s2=float(dict2[key][1])
           read_count1=float(dict1[key][3])
           read_count2=float(dict2[key][3])
           read_counts = read_count1+read_count2
#s1=math.sqrt(float(dict1[key][1]))
#s2=math.sqrt(float(dict2[key][1]))
           n1=float(dict1[key][2])
           n2=float(dict2[key][2])		   
           fenmu=math.sqrt(s1/n1+s2/n2)
           if fenmu == 0:
              Z_test = abs(avg1-avg2)
           else:
              Z_test = abs(avg1-avg2)/math.sqrt(s1/n1+s2/n2)
           P_value = Get_P(Z_test)
           fw.write(key+"\t"+str(Z_test)+"\t"+str(P_value)+"\t"+str(read_counts)+"\t"+str(avg1)+"\t"+str(avg2)+"\t"+str(s1)+"\t"+str(s2)+"\n")
    for key in dict2.keys():
        if key not in dict1.keys():
           fw.write(key+"\t"+"No z-test"+"\t"+"No p-value"+"\t"+str(dict2[key][3])+"\tNo_TN\t"+str(dict2[key][0])+"\tNo_TN\t"+str(dict2[key][1])+"\n")
    fw.close()  

'''
def Process_dict_mirBase(dict1,dict2,file_name):
    fw = open("Z_P_out/"+file_name+"_Z_P.txt","w+")
    fw.write("miRNA_ID\tmiRNA-family\tZ-test\tP-value\tread_counts\n")
    for key in dict1.keys():
        if key not in dict2.keys():
           fw.write(key+"\t"+"No z-test"+"\t"+"No p-value"+"\t"+str(dict1[key][3])+"\n")
           continue 
        else:
           avg1=float(dict1[key][0])
           avg2=float(dict2[key][0])
           s1=math.sqrt(float(dict1[key][1]))
           s2=math.sqrt(float(dict2[key][1]))
           n1=float(dict1[key][2])
           n2=float(dict2[key][2])
           fenmu=math.sqrt(s1/n1+s2/n2)
           if fenmu == 0:
              Z_test = abs(avg1-avg2)
           else:
              Z_test = abs(avg1-avg2)/math.sqrt(s1/n1+s2/n2)
           P_value = Get_P(Z_test)
           fw.write(key+"\t"+str(Z_test)+"\t"+str(P_value)+"\t"+str(dict1[key][3])+"\n")
    for key in dict2.keys():
        if key not in dict1.keys():
           fw.write(key+"\t"+"No z-test"+"\t"+"No p-value"+"\t"+str(dict2[key][3])+"\n")
    fw.close()  
'''
#input is : "mirRNA_ID\tmiRNA-family\tread_counts\tfile_num\tAverage\tStandardDeviation\n"
#output is : "miRNA_ID\tmiRNA-family\tZ-test\tP-value\tread_counts\n"
def Process_dict_mirBase_2(dict1,dict2,file_name):
    fw = open("Z_P_out/"+file_name+"_Z_P.txt","w+")
    fw.write("miRNA_ID\tisoform_coords\tZ-test\tP-value\tread_counts\tTN-avg\tNT-avg\tTN-sta\tNT-sta\n")
    for key in dict1.keys():
        if key not in dict2.keys():
           fw.write(key+"\t"+"No z-test"+"\t"+"No p-value"+"\t"+str(dict1[key][3])+"\t"+str(dict1[key][0])+"\tNo_NT\t"+str(dict1[key][2])+"\tNo_NT\n")
           continue 
        else:
           avg1=float(dict1[key][0])
           avg2=float(dict2[key][0])
           s1=float(dict1[key][1])
           s2=float(dict2[key][1])
           #s1=math.sqrt(float(dict1[key][1]))
           #s2=math.sqrt(float(dict2[key][1]))
           n1=float(dict1[key][2])
           n2=float(dict2[key][2])
           read_count1=int(dict1[key][3])
           read_count2=int(dict2[key][3])
           read_counts = read_count1+read_count2
           fenmu=math.sqrt(s1/n1+s2/n2)
           if fenmu == 0:
              Z_test = abs(avg1-avg2)
           else:
              Z_test = abs(avg1-avg2)/math.sqrt(s1/n1+s2/n2)
           P_value = Get_P(Z_test)
           fw.write(key+"\t"+str(Z_test)+"\t"+str(P_value)+"\t"+str(read_counts)+"\t"+str(avg1)+"\t"+str(avg2)+"\t"+str(s1)+"\t"+str(s2)+"\n")
    for key in dict2.keys():
        if key not in dict1.keys():
           fw.write(key+"\t"+"No z-test"+"\t"+"No p-value"+"\t"+str(dict2[key][3])+"\tNo_TN\t"+str(dict2[key][0])+"\tNo_TN\t"+str(dict2[key][1])+"\n")
    fw.close() 
'''
def Process_dict_RF(dict1,dict2,file_name):
    fw = open("Z_P_out/"+file_name+"_Z_P.txt","w+")
    fw.write("miRNA_ID\tRF-family\tZ-test\tP-value\tread_counts\n")
    for key in dict1.keys():
        if key not in dict2.keys():
           fw.write(key+"\t"+"No z-test"+"\t"+"No p-value"+"\t"+str(dict1[key][3])+"\n")
           continue 
        else:
           avg1=float(dict1[key][0])
           avg2=float(dict2[key][0])
           s1=math.sqrt(float(dict1[key][1]))
           s2=math.sqrt(float(dict2[key][1]))
           n1=float(dict1[key][2])
           n2=float(dict2[key][2])
           fenmu=math.sqrt(s1/n1+s2/n2)
           if fenmu == 0:
              Z_test = abs(avg1-avg2)
           else:
              Z_test = abs(avg1-avg2)/math.sqrt(s1/n1+s2/n2)
           P_value = Get_P(Z_test)
           fw.write(key+"\t"+str(Z_test)+"\t"+str(P_value)+"\t"+str(dict1[key][3])+"\n")
    for key in dict2.keys():
        if key not in dict1.keys():
           fw.write(key+"\t"+"No z-test"+"\t"+"No p-value"+"\t"+str(dict2[key][3])+"\n")
    fw.close() 
'''

def Process_dict_RF_2(dict1,dict2,file_name):
#output is : "miRNA_ID\tfamily\tZ-test\tP-value\tread_counts\n"
    fw = open("Z_P_out/"+file_name+"_Z_P.txt","w+")
    fw.write("miRNA_ID\tisoform_coords\tZ-test\tP-value\tread_counts\tTN-avg\tNT-avg\tTN-sta\tNT-sta\n")
    for key in dict1.keys():
        if key not in dict2.keys():
           fw.write(key+"\t"+"No z-test"+"\t"+"No p-value"+"\t"+str(dict1[key][3])+"\t"+str(dict1[key][0])+"\tNo_NT\t"+str(dict1[key][2])+"\tNo_NT\n")
           continue 
        else:
           avg1=float(dict1[key][0])
           avg2=float(dict2[key][0])
           s1=float(dict1[key][1])
           s2=float(dict2[key][1])
           #s1=math.sqrt(float(dict1[key][1]))
           #s2=math.sqrt(float(dict2[key][1]))
           n1=float(dict1[key][2])
           n2=float(dict2[key][2])
           read_count1=int(dict1[key][3])
           read_count2=int(dict2[key][3])
           read_counts=read_count1+read_count2
           fenmu=math.sqrt(s1/n1+s2/n2)
           if fenmu == 0:
              Z_test = abs(avg1-avg2)
           else:
              Z_test = abs(avg1-avg2)/math.sqrt(s1/n1+s2/n2)
           P_value = Get_P(Z_test)
           fw.write(key+"\t"+str(Z_test)+"\t"+str(P_value)+"\t"+str(read_counts)+"\t"+str(avg1)+"\t"+str(avg2)+"\t"+str(s1)+"\t"+str(s2)+"\n")
    for key in dict2.keys():
        if key not in dict1.keys():
           fw.write(key+"\t"+"No z-test"+"\t"+"No p-value"+"\t"+str(dict2[key][3])+"\tNo_TN\t"+str(dict2[key][0])+"\tNo_TN\t"+str(dict2[key][1])+"\n")
    fw.close() 


def Process_two_file_isomiR(Fir_file,Sec_file,Dir,file_name):#the Dir is the Family 
#the input is : "miRNA_ID\tisoform_coords\tmiRNA_region\tread_counts\files_num\tAverage\t
#                StandardDeviation\taccession\tproduct\tRF-family\tmiRNA-family\n"
#the output is : "miRNA_ID\tisoform_coords\tz-test\tp-value\tread_counts\tTN-avg\tNT-avg\tTN-sta\tNT-sta"
#         notes that the read_counts must the sum of the two files 
    File1={}
    File2={}
    if Fir_file.find("TN")!=-1:
        fir_file=Fir_file
        sec_file=Sec_file
    else:
        fir_file=Sec_file
        sec_file=Fir_file
    fr1 = open(Dir+"/"+fir_file)
    fr2 = open(Dir+"/"+sec_file)
    lines1 = fr1.readlines()
    lines2 = fr2.readlines()
    for line in lines1[1:]: #handle the first file to get the dict ID-iso <-> (avg,var,)
        s=line.strip().split('\t')
        index = s[0]+"\t"+s[1]  # iso level's key is the miRNA_ID\tisoform_coords
        value_avg = s[5]
        value_var = s[6]
        value_num = s[4]  #this the num is the files_num
        value_rest = s[3] #s[2]+"\t"+s[7]+"\t"+s[8]+"\t"+s[9]+"\t"+s[10]
        value = (value_avg,value_var,value_num,value_rest)
        File1[index]=value
    for line in lines2[1:]:
        s=line.strip().split('\t')
        index = s[0]+"\t"+s[1]
        value_avg = s[5]
        value_var = s[6]
        value_num = s[4]
        value_rest = s[3] #s[2]+"\t"+s[7]+"\t"+s[8]+"\t"+s[9]+"\t"+s[10]
        value = (value_avg,value_var,value_num,value_rest)
        File2[index]=value
    print str(len(File1))+"\t"+str(len(File2))+"\t"+file_name
    Process_dict_isomiR(File1,File2,file_name)
'''
def Process_two_file_mature(Fir_file,Sec_file,Dir,file_name):
    File1={}
    File2={}
    fr1 = open(Dir+"/"+Fir_file)
    fr2 = open(Dir+"/"+Sec_file)
    lines1 = fr1.readlines()
    lines2 = fr2.readlines()
    for line in lines1[1:]:
        s=line.strip().split('\t')
        index = s[0]+"\t"+s[1]
        value_avg = s[2]
        value_var = s[3]
        value_num = s[4]
        value_rest = s[5]+"\t"+s[6]+"\t"+s[7]+"\t"+s[8]+"\t"+s[9]
        value = (value_avg,value_var,value_num,value_rest)
        File1[index]=value
    for line in lines2[1:]:
        s=line.strip().split('\t')
        index = s[0]+"\t"+s[1]
        value_avg = s[2]
        value_var = s[3]
        value_num = s[4]
        value_rest = s[5]+"\t"+s[6]+"\t"+s[7]+"\t"+s[8]+"\t"+s[9]
        value = (value_avg,value_var,value_num,value_rest)
        File2[index]=value
    print str(len(File1))+"\t"+str(len(File2))+"\t"+file_name
    Process_dict_mature(File1,File2,file_name)
'''
#input is : miRNA_ID\tread_counts\tfiles_number\tAverage\tStardarddeviation\n
def Process_two_file_precursor(Fir_file,Sec_file,Dir,file_name):
    File1={}
    File2={}
    if Fir_file.find("TN") != -1:
        fir_file=Fir_file
        sec_file=Sec_file
    else:
        fir_file=Sec_file
        sec_file=Fir_file
    fr1 = open(Dir+"/"+fir_file)
    fr2 = open(Dir+"/"+sec_file)
    lines1 = fr1.readlines()
    lines2 = fr2.readlines()
    for line in lines1[1:]:
        s=line.strip().split('\t')
        index = s[0]
        value_avg = s[3]
        value_var = s[4]
        value_num = s[2]
        value_rest = s[1]
        value = (value_avg,value_var,value_num,value_rest)
        File1[index]=value
    for line in lines2[1:]:
        s=line.strip().split('\t')
        index = s[0]
        value_avg = s[3]
        value_var = s[4]
        value_num = s[2]
        value_rest = s[1]
        value = (value_avg,value_var,value_num,value_rest)
        File2[index]=value
    print str(len(File1))+"\t"+str(len(File2))+"\t"+file_name
    Process_dict_precursor(File1,File2,file_name)

def Process_two_file_mature_2(Fir_file,Sec_file,Dir,file_name):
#input is : "miRNA_ID\tmiRNA_region\tread_counts\tfiles_number\tAverage\tStardarddeviation\n"
    File1={}
    File2={}
    if Fir_file.find("TN")!=-1:
        fir_file=Fir_file
        sec_file=Sec_file
    else:
        fir_file=Sec_file
        sec_file=Fir_file
    fr1 = open(Dir+"/"+fir_file)
    fr2 = open(Dir+"/"+sec_file)
    lines1 = fr1.readlines()
    lines2 = fr2.readlines()
    for line in lines1[1:]:
        s=line.strip().split('\t')
        index = s[0]+"\t"+s[1]
        value_avg = s[4]
        value_var = s[5]
        value_num = s[3]
        value_rest = s[2]
        value = (value_avg,value_var,value_num,value_rest)
        File1[index]=value
    for line in lines2[1:]:
        s=line.strip().split('\t')
        index = s[0]+"\t"+s[1]
        value_avg = s[4]
        value_var = s[5]
        value_num = s[3]
        value_rest = s[2]
        value = (value_avg,value_var,value_num,value_rest)
        File2[index]=value
    print str(len(File1))+"\t"+str(len(File2))+"\t"+file_name
    Process_dict_mature_2(File1,File2,file_name)

'''
def Process_two_file_mirBase(Fir_file,Sec_file,Dir,file_name):
    File1={}
    File2={}
    fr1 = open(Dir+"/"+Fir_file)
    fr2 = open(Dir+"/"+Sec_file)
    lines1 = fr1.readlines()
    lines2 = fr2.readlines()
    for line in lines1[1:]:
        s=line.strip().split('\t')
        index = s[0]+"\t"+s[5]
        value_avg = s[1] 
        value_var = s[2]
        value_num = s[3]
        value_rest = s[4]
        value = (value_avg,value_var,value_num,value_rest)
        File1[index]=value
    for line in lines2[1:]:
        s=line.strip().split('\t')
        index = s[0]+"\t"+s[5]
        value_avg = s[1]
        value_var = s[2]
        value_num = s[3]
        value_rest = s[4]
        value = (value_avg,value_var,value_num,value_rest)
        File2[index]=value
    print str(len(File1))+"\t"+str(len(File2))+"\t"+file_name
    Process_dict_mirBase(File1,File2,file_name)
'''
#input is : "mirRNA_ID\tmiRNA-family\tread_counts\tfile_num\tAverage\tStandardDeviation\n"
def Process_two_file_mirBase_2(Fir_file,Sec_file,Dir,file_name):
    File1={}
    File2={}
    if Fir_file.find("TN")!=-1:
        fir_file=Fir_file
        sec_file=Sec_file
    else:
        fir_file=Sec_file
        sec_file=Fir_file
    fr1 = open(Dir+"/"+fir_file)
    fr2 = open(Dir+"/"+sec_file)
    lines1 = fr1.readlines()
    lines2 = fr2.readlines()
    for line in lines1[1:]:
        s=line.strip().split('\t')
        index = s[0]+"\t"+s[1] 
        value_avg = s[4] 
        value_var = s[5]
        value_num = s[3]
        value_rest = s[2]
        value = (value_avg,value_var,value_num,value_rest)
        File1[index]=value
    for line in lines2[1:]:
        s=line.strip().split('\t')
        index = s[0]+"\t"+s[1]
        value_avg = s[4]
        value_var = s[5]
        value_num = s[3]
        value_rest = s[2]
        value = (value_avg,value_var,value_num,value_rest)
        File2[index]=value
    print str(len(File1))+"\t"+str(len(File2))+"\t"+file_name
    Process_dict_mirBase_2(File1,File2,file_name)

'''
def Process_two_file_RF(Fir_file,Sec_file,Dir,file_name):
    File1={}
    File2={}
    fr1 = open(Dir+"/"+Fir_file)
    fr2 = open(Dir+"/"+Sec_file)
    lines1 = fr1.readlines()
    lines2 = fr2.readlines()
    for line in lines1[1:]:
        s=line.strip().split('\t')
        index = s[0]+"\t"+s[5]
        value_avg = s[1] 
        value_var = s[2]
        value_num = s[3]
        value_rest = s[4]
        value = (value_avg,value_var,value_num,value_rest)
        File1[index]=value
    for line in lines2[1:]:
        s=line.strip().split('\t')
        index = s[0]+"\t"+s[5]
        value_avg = s[1]
        value_var = s[2]
        value_num = s[3]
        value_rest = s[4]
        value = (value_avg,value_var,value_num,value_rest)
        File2[index]=value
    print str(len(File1))+"\t"+str(len(File2))+"\t"+file_name
    Process_dict_RF(File1,File2,file_name)
'''

def Process_two_file_RF_2(Fir_file,Sec_file,Dir,file_name):
#input is : "mirRNA_ID\tmiRNA-family\tread_counts\tfile_num\tAverage\tStandardDeviation\n"
#output is : "miRNA_ID\tfamily\tZ-test\tP-value\tread_counts\n"
    File1={}
    File2={}
    if Fir_file.find("TN")!=-1:
        fir_file=Fir_file
        sec_file=Sec_file
    else:
        fir_file=Sec_file
        sec_file=Fir_file
    fr1 = open(Dir+"/"+fir_file)
    fr2 = open(Dir+"/"+sec_file)
    lines1 = fr1.readlines()
    lines2 = fr2.readlines()
    for line in lines1[1:]:
        s=line.strip().split('\t')
        index = s[0]+"\t"+s[1] 
        value_avg = s[4] 
        value_var = s[5]
        value_num = s[3]
        value_rest = s[2]
        value = (value_avg,value_var,value_num,value_rest)
        File1[index]=value
    for line in lines2[1:]:
        s=line.strip().split('\t')
        index = s[0]+"\t"+s[1] 
        value_avg = s[4]
        value_var = s[5]
        value_num = s[3]
        value_rest = s[2]
        value = (value_avg,value_var,value_num,value_rest)
        File2[index]=value
    print str(len(File1))+"\t"+str(len(File2))+"\t"+file_name
    Process_dict_RF_2(File1,File2,file_name)

#the intput is : "miRNA_ID\tiosform_coords\tmiRNA_region\tread_counts\tline_number\tAverage\t
#                 StandardDeviation\taccession\tproduct\tRF-family\tmiRNA-family\n"
def z_test_p_value_isomiR(input):
   files =os.listdir(input)# Family/Breast_acc_pro_family.txt
   Process_files=set()
   Files_list=[]
   if os.path.exists("Z_P_out"):
       shutil.rmtree("Z_P_out")
       os.makedirs("Z_P_out")
   else:
       os.makedirs("Z_P_out")
   for file in files:
       Files_list.append(file)#save all the file names
   for i in range(len(Files_list)):  # search the pair files to deal with 
        index =Files_list[i].find('_')
        Fir_file_name = Files_list[i][:index]
        Fir_file =Files_list[i]
        if Fir_file_name in Process_files:
           continue
        else:
           Process_files.add(Fir_file_name)# the file already processed
        for j in range(i+1,len(Files_list)):
            index = Files_list[j].find('_')
            Sec_file_name = Files_list[j][:index]
            Sec_file = Files_list[j]
            if Sec_file_name != Fir_file_name:
                continue
            else:
                #print "now process\t"+Fir_file+"\t"+Sec_file
                Process_two_file_isomiR(Fir_file,Sec_file,input,Sec_file_name) 
'''
def z_test_p_value_mature(input):
   files =os.listdir(input)
   Process_files=set()
   Files_list=[]
   if os.path.exists("Z_P_out"):
       shutil.rmtree("Z_P_out")
       os.makedirs("Z_P_out")
   else:
       os.makedirs("Z_P_out")
   for file in files:
       Files_list.append(file)
   for i in range(len(Files_list)):
        index =Files_list[i].find('_')
        Fir_file_name = Files_list[i][:index]
        Fir_file =Files_list[i]
        if Fir_file_name in Process_files:
           continue
        else:
           Process_files.add(Fir_file_name)
        for j in range(i+1,len(Files_list)):
            index = Files_list[j].find('_')
            Sec_file_name = Files_list[j][:index]
            Sec_file = Files_list[j]
            if Sec_file_name != Fir_file_name:
                continue
            else:
                print "now process\t"+Fir_file+"\t"+Sec_file
                Process_two_file_mature(Fir_file,Sec_file,input,Sec_file_name)
'''
#intput is : miRNA_ID\tread_counts\tfiles_number\tAverage\tStardarddeviation\n
def z_test_p_value_precursor(input):
   files =os.listdir(input)
   Process_files=set()
   Files_list=[]
   if os.path.exists("Z_P_out"):
       shutil.rmtree("Z_P_out")
       os.makedirs("Z_P_out")
   else:
       os.makedirs("Z_P_out")
   for file in files:
       Files_list.append(file)
   for i in range(len(Files_list)):
        index =Files_list[i].find('_')
        Fir_file_name = Files_list[i][:index]
        Fir_file =Files_list[i]
        if Fir_file_name in Process_files:
           continue
        else:
           Process_files.add(Fir_file_name)
        for j in range(i+1,len(Files_list)):
            index = Files_list[j].find('_')
            Sec_file_name = Files_list[j][:index]
            Sec_file = Files_list[j]
            if Sec_file_name != Fir_file_name:
                continue
            else:
                print "now process\t"+Fir_file+"\t"+Sec_file
                Process_two_file_precursor(Fir_file,Sec_file,input,Sec_file_name)

def z_test_p_value_mature_2(input):
#input is : "miRNA_ID\tmiRNA_region\tread_counts\tfiles_number\tAverage\tStardarddeviation\n"
   files =os.listdir(input)
   Process_files=set()
   Files_list=[]
   if os.path.exists("Z_P_out"):
       shutil.rmtree("Z_P_out")
       os.makedirs("Z_P_out")
   else:
       os.makedirs("Z_P_out")
   for file in files:
       Files_list.append(file)
   for i in range(len(Files_list)):
        index =Files_list[i].find('_')
        Fir_file_name = Files_list[i][:index]
        Fir_file =Files_list[i]
        if Fir_file_name in Process_files:
           continue
        else:
           Process_files.add(Fir_file_name)
        for j in range(i+1,len(Files_list)):
            index = Files_list[j].find('_')
            Sec_file_name = Files_list[j][:index]
            Sec_file = Files_list[j]
            if Sec_file_name != Fir_file_name:
                continue
            else:
                print "now process\t"+Fir_file+"\t"+Sec_file
                Process_two_file_mature_2(Fir_file,Sec_file,input,Sec_file_name)

'''
def z_test_p_value_mirBase(input):
   files =os.listdir(input)
   Process_files=set()
   Files_list=[]
   if os.path.exists("Z_P_out"):
       shutil.rmtree("Z_P_out")
       os.makedirs("Z_P_out")
   else:
       os.makedirs("Z_P_out")
   for file in files:
       Files_list.append(file)
   for i in range(len(Files_list)):
        index =Files_list[i].find('_')
        Fir_file_name = Files_list[i][:index]
        Fir_file =Files_list[i]
        if Fir_file_name in Process_files:
           continue
        else:
           Process_files.add(Fir_file_name)
        for j in range(i+1,len(Files_list)):
            index = Files_list[j].find('_')
            Sec_file_name = Files_list[j][:index]
            Sec_file = Files_list[j]
            if Sec_file_name != Fir_file_name:
                continue
            else:
                print "now process\t"+Fir_file+"\t"+Sec_file
                Process_two_file_mirBase(Fir_file,Sec_file,input,Sec_file_name)
'''
#input is : "mirRNA_ID\tmiRNA-family\tread_counts\tfile_num\tAverage\tStandardDeviation\n"
def z_test_p_value_mirBase_2(input):
   files =os.listdir(input)
   Process_files=set()
   Files_list=[]
   if os.path.exists("Z_P_out"):
       shutil.rmtree("Z_P_out")
       os.makedirs("Z_P_out")
   else:
       os.makedirs("Z_P_out")
   for file in files:
       Files_list.append(file)
   for i in range(len(Files_list)):
        index =Files_list[i].find('_')
        Fir_file_name = Files_list[i][:index]
        Fir_file =Files_list[i]
        if Fir_file_name in Process_files:
           continue
        else:
           Process_files.add(Fir_file_name)
        for j in range(i+1,len(Files_list)):
            index = Files_list[j].find('_')
            Sec_file_name = Files_list[j][:index]
            Sec_file = Files_list[j]
            if Sec_file_name != Fir_file_name:
                continue
            else:
                print "now process\t"+Fir_file+"\t"+Sec_file
                Process_two_file_mirBase_2(Fir_file,Sec_file,input,Sec_file_name)

'''
def z_test_p_value_RF(input):
   files =os.listdir(input)
   Process_files=set()
   Files_list=[]
   if os.path.exists("Z_P_out"):
       shutil.rmtree("Z_P_out")
       os.makedirs("Z_P_out")
   else:
       os.makedirs("Z_P_out")
   for file in files:
       Files_list.append(file)
   for i in range(len(Files_list)):
        index =Files_list[i].find('_')
        Fir_file_name = Files_list[i][:index]
        Fir_file =Files_list[i]
        if Fir_file_name in Process_files:
           continue
        else:
           Process_files.add(Fir_file_name)
        for j in range(i+1,len(Files_list)):
            index = Files_list[j].find('_')
            Sec_file_name = Files_list[j][:index]
            Sec_file = Files_list[j]
            if Sec_file_name != Fir_file_name:
                continue
            else:
                print "now process\t"+Fir_file+"\t"+Sec_file
                Process_two_file_RF(Fir_file,Sec_file,input,Sec_file_name)
'''

def z_test_p_value_RF_2(input):
#output is : "miRNA_ID\tfamily\tZ-test\tP-value\tread_counts\n"
   files =os.listdir(input)
   Process_files=set()
   Files_list=[]
   if os.path.exists("Z_P_out"):
       shutil.rmtree("Z_P_out")
       os.makedirs("Z_P_out")
   else:
       os.makedirs("Z_P_out")
   for file in files:
       Files_list.append(file)
   for i in range(len(Files_list)):
        index =Files_list[i].find('_')
        Fir_file_name = Files_list[i][:index]
        Fir_file =Files_list[i]
        if Fir_file_name in Process_files:
           continue
        else:
           Process_files.add(Fir_file_name)
        for j in range(i+1,len(Files_list)):
            index = Files_list[j].find('_')
            Sec_file_name = Files_list[j][:index]
            Sec_file = Files_list[j]
            if Sec_file_name != Fir_file_name:
                continue
            else:
                print "now process\t"+Fir_file+"\t"+Sec_file
                Process_two_file_RF_2(Fir_file,Sec_file,input,Sec_file_name)

