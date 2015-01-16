import sys 
import re
import os
import shutil

def Get_RF_MIR(DR):
    MIR="No MIR_family"
    RF="No RF_family"
    for item in DR:
        if item.find("RFAM")!=-1:
           f = item.strip().split(";")
           return f[1].strip(),f[2].strip()[:-1]
        else:
           item = item.lower()
           if MIR!="No MIR_family":
              continue
           else:
              pattern = re.compile('mir((-[0-9]+)|([0-9]+))[a-z]?')
              match = re.search(pattern,item)
              if match:
                 MIR=match.group()
    return RF,MIR          

def Process(DR,FT,accession_product,accession_family):
    i=0
    Accession=[]
    # the FT item must have the accession or the product
    while(i<len(FT)):
        accession = Get_the_middle_world('"','"',FT[i])
        product = Get_the_middle_world('"','"',FT[i+1])
        accession_product[accession]=product
        Accession.append(accession)
        i=i+2
    RF,MIR = Get_RF_MIR(DR)
    for acc in Accession:
        accession_family[acc]=RF+";"+MIR
    

def Get_the_middle_world(start_str,end_str,data):
    start = data.find(start_str)
    if start >=0: 
        start +=len(start_str)
        end = data.find(end_str,start)
        if end>=0:
           return data[start:end].strip()

def Process_File_isomiR(input,output,accession_product,accession_family):
#"miRNA_ID\tiosform_coords\tmiRNA_region\tread_counts\tline_number\tAverage\tStandardDeviation\n"
    fr = open(input,'r')
    fw = open(output,'w+')
    line = fr.next()
#"miRNA_ID\tiosform_coords\tmiRNA_region\tread_counts\tline_number\tAverage\t
#StandardDeviation\taccession\tproduct\tRF-family\tmiRNA-family\n"
    fw.write(line.strip()+"\taccession\tproduct\tRF-family\tmiRNA-family"+"\n")
    for line in fr:
        s = line.strip().split("\t")
        key = s[2].split(',')[1].strip() #the key comes from the miRNA_region 
        if accession_product.has_key(key): # get the accession product
           product = str(accession_product[key])[-2:]
           if product =="5p" or product =="3p":
               line = line.strip()+"\t"+str(key)+"\t"+product
           else: 
               line = line.strip()+"\t"+str(key)+"\t"+"No_product"
        else:
           line = line.strip()+"\t"+"No_accession"+"\t"+"No_product"
        if accession_family.has_key(key):#get the family
            s = str(accession_family[key]).split(";")
            line = line.strip()+"\t"+s[0]+"\t"+s[1]+"\n"
        else:
            line = line.strip()+"\t"+"No RF_family"+"\tNo MIR_family"+"\n"
        fw.write(line)
    fr.close()
    fw.close()
'''
def Process_File_mature(input,output,accession_product,accession_family):
    fr = open(input,'r')
    fw = open(output,'w+')
    line = fr.next()
    fw.write(line.strip()+"\taccession\tproduct\tRF-family\tmiRNA-family"+"\n")
    for line in fr:
        s = line.strip().split("\t")
        key = s[1].split(',')[1].strip()
        if accession_product.has_key(key):
           product = str(accession_product[key])[-2:]
           if product =="5p" or product =="3p":
               line = line.strip()+"\t"+str(key)+"\t"+product
           else: 
               line = line.strip()+"\t"+str(key)+"\t"+"No_product"
        else:
           line = line.strip()+"\t"+"No_accession"+"\t"+"No_product"
        if accession_family.has_key(key):
            s = str(accession_family[key]).split(";")
            line = line.strip()+"\t"+s[0]+"\t"+s[1]+"\n"
        else:
            line = line.strip()+"\t"+"No RF_family"+"\tNo MIR_family"+"\n"
        fw.write(line)
    fr.close()
    fw.close()
'''
def Process_File_mirBase(input,output,accession_product,accession_family):
#input is line[5]+rates 
#output is line[5]+rates+"\taccession\tproduct\tRF-family\tmiRNA-family"+"\n"
    fr = open(input,'r')
    fw = open(output,'w+')
    line = fr.next()
    fw.write(line.strip()+"\taccession\tproduct\tRF-family\tmiRNA-family"+"\n")
    for line in fr:
        s = line.strip().split("\t")
        key = s[5].split(',')[1].strip()#s[5] is the region 
        if accession_product.has_key(key):
           product = str(accession_product[key])[-2:]
           if product =="5p" or product =="3p":
               line = line.strip()+"\t"+str(key)+"\t"+product
           else: 
               line = line.strip()+"\t"+str(key)+"\t"+"No_product"
        else:
           line = line.strip()+"\t"+"No_accession"+"\t"+"No_product"
        if accession_family.has_key(key):
            s = str(accession_family[key]).split(";")
            line = line.strip()+"\t"+s[0]+"\t"+s[1]+"\n"
        else:
            line = line.strip()+"\t"+"No RF_family"+"\tNo MIR_family"+"\n"
        fw.write(line)
    fr.close()
    fw.close()

def Process_File_RF(input,output,accession_product,accession_family):
    fr = open(input,'r')
    fw = open(output,'w+')
    line = fr.next()
    fw.write(line.strip()+"\taccession\tproduct\tRF-family\miRNA-family"+"\n")
    for line in fr:
        s = line.strip().split("\t")
        key = s[5].split(',')[1].strip()
        if accession_product.has_key(key):
           product = str(accession_product[key])[-2:]
           if product =="5p" or product =="3p":
               line = line.strip()+"\t"+str(key)+"\t"+product
           else: 
               line = line.strip()+"\t"+str(key)+"\t"+"No_product"
        else:
           line = line.strip()+"\t"+"No_accession"+"\t"+"No_product"
        if accession_family.has_key(key):
            s = str(accession_family[key]).split(";")
            line = line.strip()+"\t"+s[0]+"\t"+s[1]+"\n"
        else:
            line = line.strip()+"\t"+"No RF_family"+"\tNo MIR_family"+"\n"
        fw.write(line)
    fr.close()
    fw.close()

def process_accession_product_isomiR(info_file,input):# the input is the result/Breast_NT.txt 
    fr = open(info_file,'r')
    DR=[]
    FT=[]
    accession_product={}
    accession_family={}
    for line in fr: 
        if line.strip()=="//":# a batch data sign
           Process(DR,FT,accession_product,accession_family)# get the dictionary accession-product , accession-family 
           DR=[]
           FT=[]
        else:
           if line[:2]=="DR":
              DR.append(line.strip())
           if line[:2]=="FT" and line.find("accession")!=-1:
              FT.append(line.strip())
           if line[:2]=="FT" and line.find("product")!=-1:
              FT.append(line.strip())
    fr.close()
    fw = open("Accession_Product.txt","w+")
    fw.write("Accession\tProduct\n")
    for key in accession_product:
        fw.write(str(key)+"\t"+str(accession_product[key])+"\n") 
    fw.close()
    fw = open("Accession_Family.txt","w+")
    fw.write("Accession\tRF-family\tMIR-family\n")
    for key in accession_family:
        family = str(accession_family[key]).split(";")
        fw.write(str(key)+"\t"+family[0]+"\t"+family[1]+"\n")
    fw.close()
    files = os.listdir(input) # the input is the result/Breast_NT.txt 
    if os.path.exists("Family"):
        shutil.rmtree("Family")
        os.makedirs("Family")
    else:
        os.makedirs("Family")
    for file in files:
         index = file.find('.')
         filename = file[:index]+"_acc_pro_family.txt" #this is the output file name 
         Process_File_isomiR(input+"/"+file,"Family/"+filename,accession_product,accession_family)

'''
def process_accession_product_mature(info_file,input):
    fr = open(info_file,'r')
    DR=[]
    FT=[]
    accession_product={}
    accession_family={}
    for line in fr: 
        if line.strip()=="//":
           Process(DR,FT,accession_product,accession_family)
           DR=[]
           FT=[]
        else:
           if line[:2]=="DR":
              DR.append(line.strip())
           if line[:2]=="FT" and line.find("accession")!=-1:
              FT.append(line.strip())
           if line[:2]=="FT" and line.find("product")!=-1:
              FT.append(line.strip())
    fr.close()
    fw = open("Accession_Product.txt","w+")
    for key in accession_product:
           fw.write(str(key)+"\t"+str(accession_product[key])+"\n") 
    fw.close()
    fw = open("Accession_Family.txt","w+")
    for key in accession_family:
         family = str(accession_family[key]).split(";")
         fw.write(str(key)+"\t"+family[0]+"\t"+family[1]+"\n")
    fw.close()
    files = os.listdir(input)
    if os.path.exists("Family"):
        shutil.rmtree("Family")
        os.makedirs("Family")
    else:
        os.makedirs("Family")
    for file in files:
         index = file.find('.')
         filename = file[:index]+"_acc_pro_family.txt"
         Process_File_mature(input+"/"+file,"Family/"+filename,accession_product,accession_family)
'''

def process_accession_product_mirBase(info_file,input):
    fr = open(info_file,'r')
    DR=[]
    FT=[]
    accession_product={}
    accession_family={}
    for line in fr: 
        if line.strip()=="//":
           Process(DR,FT,accession_product,accession_family)
           DR=[]
           FT=[]
        else:
           if line[:2]=="DR":
              DR.append(line.strip())
           if line[:2]=="FT" and line.find("accession")!=-1:
              FT.append(line.strip())
           if line[:2]=="FT" and line.find("product")!=-1:
              FT.append(line.strip())
    fr.close()
    fw = open("Accession_Product.txt","w+")
    for key in accession_product:
           fw.write(str(key)+"\t"+str(accession_product[key])+"\n") 
    fw.close()
    fw = open("Accession_Family.txt","w+")
    for key in accession_family:
         family = str(accession_family[key]).split(";")
         fw.write(str(key)+"\t"+family[0]+"\t"+family[1]+"\n")
    fw.close()
    files = os.listdir(input)
    if os.path.exists("Family/"+input):
        shutil.rmtree("Family/"+input)
        os.makedirs("Family/"+input)
    else:
        os.makedirs("Family/"+input)
    for file in files:
         index = file.find('.')
         filename = file[:index]+"_acc_pro_family.txt"
         Process_File_mirBase(input+"/"+file,"Family/"+input+"/"+filename,accession_product,accession_family)

def process_accession_product_RF(info_file,input):
    fr = open(info_file,'r')
    DR=[]
    FT=[]
    accession_product={}
    accession_family={}
    for line in fr: 
        if line.strip()=="//":
           Process(DR,FT,accession_product,accession_family)
           DR=[]
           FT=[]
        else:
           if line[:2]=="DR":
              DR.append(line.strip())
           if line[:2]=="FT" and line.find("accession")!=-1:
              FT.append(line.strip())
           if line[:2]=="FT" and line.find("product")!=-1:
              FT.append(line.strip())
    fr.close()
    fw = open("Accession_Product.txt","w+")
    for key in accession_product:
           fw.write(str(key)+"\t"+str(accession_product[key])+"\n") 
    fw.close()
    fw = open("Accession_Family.txt","w+")
    for key in accession_family:
         family = str(accession_family[key]).split(";")
         fw.write(str(key)+"\t"+family[0]+"\t"+family[1]+"\n")
    fw.close()
    files = os.listdir(input)
    if os.path.exists("Family/"+input):
        shutil.rmtree("Family/"+input)
        os.makedirs("Family/"+input)
    else:
        os.makedirs("Family/"+input)
    for file in files:
         index = file.find('.')
         filename = file[:index]+"_acc_pro_family.txt"
         Process_File_RF(input+"/"+file,"Family/"+input+"/"+filename,accession_product,accession_family)





