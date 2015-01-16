import sys 
import os
import shutil
import subprocess
import preprocess
import process_accession_product
import math 
def Process_Fold_isomiR(input_1,input_2,output_file):  #the process_zip will call this function with the correct level
   files = os.listdir(input_1+" "+input_2)  # find the files 
   output =output_file  # the outoput_files will mkdir and new fold named Breast_NT in the current fold to save the new_files 
   if os.path.exists("result"):  #the result will save the updata files 
       #shutil.rmtree("result")
       #os.makedirs("result")
       pass
   else:
       os.makedirs("result")
   if os.path.exists(output):# this the ouput is the Breast_NT 
       shutil.rmtree(output)
       os.makedirs(output)
   else:
       os.makedirs(output)
   print 
   print "now is process"+output
   print 
   for file in files:
       preprocess.PreProcess_isomiR(input_1,input_2+"/"+file,output+"/new"+file)  #preprocess the isomiR once handle a file a time 
   files = os.listdir(output) # get the new_file
   fw=open("result/"+output+".txt",'w+')  # as the output is Breast_NT so we will get the result/Breast_NT.txt to save the combine files' data
   fw.write("miRNA_ID\tiosform_coords\tmiRNA_region\tread_counts\tfiles_num\tAverage\tStandardDeviation\n")
   pro_list=[]
   files_num=0
   for file in files:
        is_first=1
        files_num+=1
        fr=open(output+"/"+file,'r')
        for line in fr: # iter each files 
           if is_first:
              is_first=0 #ignore the first line 
              continue
           list=line.strip().split()# " miRNA_ID[0] \t isoform_coords[1] \t read_count[2] \t miRNA_region[3] \t rates[4] " 
           pro_list.append(list)
        is_first=1
   shutil.rmtree(output) # delete the output/new_file
   pro_list.sort(key=lambda x:(x[0],x[1])) # all the files data are in the pro_list , sort as the ID and the isoform 
   pro_list2=[]  # save the same ID and the isoform data 
   is_first=1  

   def Avge_var_count_isomiR(process_list,files_num):#calculate the average in the process_list[[list],[list],[list]] 
        # each list in process_list is : " miRNA_ID[0] \t isoform_coords[1] \t read_count[2] \t miRNA_region[3] \t rates[4] " 
       sum=0.0
       read_counts=0
       for item in process_list:
           sum=sum+float(item[4])
           read_counts=read_counts+int(item[2])
       avg = sum/files_num
       var=0.0
       for item in process_list:
          var=var+(float(item[4])-avg)**2
       var=math.sqrt(var/files_num) 
	   #the output is  :  "miRNA_ID\tiosform_coords\tmiRNA_region\tread_counts\tline_number\tAverage\tStandardDeviation\n"
       fw.write(process_list[0][0]+"\t"+process_list[0][1]+"\t"+process_list[0][3]+"\t"+str(read_counts)+"\t"+str(files_num)+"\t"+str(avg)+"\t"+str(var)+"\n")

   num=1 #output data have how many lines 
   for item in pro_list: # the data in pro_list is aloso the list
       if is_first:
          is_first=0
          pre_first=item[0]
          pre_second=item[1]
          pro_list2.append(item)
          continue 
       first=item[0]
       second=item[1]
       if first!=pre_first or second != pre_second :
          Avge_var_count_isomiR(pro_list2,files_num) # call the function to calcute the average 
          num+=1
          pre_first=first
          pre_second=second
          pro_list2=[]
          pro_list2.append(item)
       else:
          pro_list2.append(item)
   Avge_var_count_isomiR(pro_list2,files_num)
   print "already process "+str(num)

'''
def Process_Fold_mature(input_1,input_2,output_file):
   files = os.listdir(input_1+" "+input_2)
   output =output_file
   if os.path.exists("result"):
       pass
       # shutil.rmtree("result")
       # os.makedirs("result")
   else:
        os.makedirs("result")
   if os.path.exists(output):
       shutil.rmtree(output)
       os.makedirs(output)
   else:
       os.makedirs(output)
   count=0
   print 
   print "now is process"+output
   print 
   for file in files:
       count+=1
       preprocess.PreProcess(input_1,input_2+"/"+file,output+"/new"+file)
   print "the totle number of files is : "+str(count)
   files = os.listdir(output)
   fw=open("result/"+output_file+".txt",'w+')
   fw.write("miRNA_ID\tmiRNA_region\tAverage\tVariance\tline_number\tread_counts\n")#this is the change part
   pro_list=[]
   for file in files:
        is_first=1
        fr=open(output+"/"+file,'r')
        for line in fr:
           if is_first:
              is_first=0
              continue
           list=line.strip().split()
           pro_list.append(list)
        is_first=1
   shutil.rmtree(output_file)
   pro_list.sort(key=lambda x:(x[0],x[5])) #this is the change part
   pro_list2=[]
   is_first=1

   def Avge_var_count_mature(process_list,count):
       count_=0
       sum=0.0
       read_counts=0
       for item in process_list:
           count_+=1
           sum=sum+float(item[6])
           read_counts=read_counts+int(item[2])
       avg = sum/count_
       var=0.0
       for item in process_list:
          var=var+(float(item[6])-avg)**2
       var=var/count_ 
       fw.write(process_list[0][0]+"\t"+process_list[0][5]+"\t"+str(avg)+"\t"+str(var)+"\t"+str(count_)+"\t"+str(read_counts)+"\n") #this is the change part

   count=0
   num=1
   for item in pro_list:
       if is_first:
          is_first=0
          pre_first=item[0]
          pre_second=item[5] # this is the change part
          pro_list2.append(item)
          continue 
       first=item[0]
       second=item[5] #this is the change part
       if first!=pre_first or second != pre_second :
          Avge_var_count_mature(pro_list2,count)
          num+=1
          count=1
          pre_first=first
          pre_second=second
          pro_list2=[]
          pro_list2.append(item)
       else:
          count=count+1
          pro_list2.append(item)
   print "already process "+str(num)
'''
def Process_Fold_precursor(input_1,input_2,output_file):
   files = os.listdir(input_1+" "+input_2)
   output =output_file
   if os.path.exists("result"):
        pass
       #shutil.rmtree("result")
       #os.makedirs("result")
   else:
       os.makedirs("result")
   if os.path.exists(output):
       shutil.rmtree(output)
       os.makedirs(output)
   else:
       os.makedirs(output)
   count=0
   print 
   print "now is process"+output
   print 
   for file in files:
       count+=1
       preprocess.PreProcess_precursor(input_1,input_2+"/"+file,output+"/new"+file)#this is the change part
   print "the totle number of files is : "+str(count)
   files = os.listdir(output)
   #the input is "miRNA_ID \t read_counts \t rates"
   fw=open("result/"+output_file+".txt",'w+')
   fw.write("miRNA_ID\tread_counts\tfiles_number\tAverage\tStandardDeviation\n")#this is the change part
   pro_list=[]
   files_num=0
   for file in files:
        is_first=1
        files_num+=1
        fr=open(output+"/"+file,'r')
        for line in fr:
           if is_first:#ignore the first line 
              is_first=0
              continue
           list=line.strip().split()
           pro_list.append(list)
        is_first=1
   shutil.rmtree(output_file)
   pro_list.sort(key=lambda x:(x[0])) #pro_list is the [[],[],[]] , and sorted by ID 
   
   def Avge_var_count_precursor(process_list,files_num):#input is : "miRNA_ID \t read_counts \t rates"
   #output is : miRNA_ID\tread_counts\tfiles_number\tAverage\tStardarddeviation\n
       sum=0.0
       read_counts=0.0
       for item in process_list:
           sum=sum+float(item[2])
           read_counts=read_counts+float(item[1])
       avg = sum/files_num
       var=0.0
       for item in process_list:
          var=var+(float(item[2])-avg)**2
       var=math.sqrt(var/files_num)
       fw.write(process_list[0][0]+"\t"+str(read_counts)+"\t"+str(files_num)+"\t"+str(avg)+"\t"+str(var)+"\n") #this is the change part

   pro_list2=[]#save the batch data to get average
   is_first=1
   num=1
   for item in pro_list: #item is a list  
       if is_first:
          is_first=0
          pre_first=item[0]
          pro_list2.append(item)
          continue 
       first=item[0]
       if first!=pre_first:
          Avge_var_count_precursor(pro_list2,files_num)
          num+=1
          pre_first=first
          pro_list2=[]
          pro_list2.append(item)
       else:
          pro_list2.append(item)
   Avge_var_count_precursor(pro_list2,files_num)
   print "already process "+str(num)

def Process_Fold_mature_2(input_1,input_2,output_file):
   files = os.listdir(input_1+" "+input_2)
   output =output_file  # the output is the Bread_NT , which saves the new_files
   if os.path.exists("result"):
        pass
       #shutil.rmtree("result")
       #os.makedirs("result")
   else:
       os.makedirs("result")
   if os.path.exists(output):
       shutil.rmtree(output)
       os.makedirs(output)
   else:
       os.makedirs(output)
   print 
   print "now is process"+output
   print 
   for file in files:
       preprocess.PreProcess_mature_2(input_1,input_2+"/"+file,output+"/new"+file)
   files = os.listdir(output)  #the output is the Bread_NT/new_files
   fw=open("result/"+output+".txt",'w+')
   fw.write("miRNA_ID\tmiRNA_region\tread_counts\tfiles_number\tAverage\tStardarddeviation\n")#this is the change part
   pro_list=[]
   files_num=0
   for file in files:#miRNA_ID[0] \t miRNA_region[1] \t read_counts[2] \t rates[3]
        is_first=1
        files_num+=1
        fr=open(output+"/"+file,'r')
        for line in fr:
           if is_first:
              is_first=0
              continue
           list=line.strip().split("\t") 
		   #this part is to delete the star and the MIMAT
           #list[1]=list[1].split(",")[1]
           pro_list.append(list)
        is_first=1
   shutil.rmtree(output_file)
   pro_list.sort(key=lambda x:(x[0],x[1])) #sorted according to the ID and region 

   def Avge_var_count_mature_2(process_list,files_num):# get the average
   #input is : miRNA_ID[0] \t miRNA_region[1] \t read_counts[2] \t rates[3]
   #output is : "miRNA_ID\tmiRNA_region\tread_counts\tfiles_number\tAverage\tStardarddeviation\n"
       sum=0.0
       read_counts=0
       for item in process_list:
           sum=sum+float(item[3])
           read_counts=read_counts+int(item[2])
       avg = sum/files_num
       var=0.0
       for item in process_list:
          var=var+(float(item[3])-avg)**2
       var=math.sqrt(var/files_num) 
       fw.write(process_list[0][0]+"\t"+process_list[0][1]+"\t"+str(read_counts)+"\t"+str(files_num)+"\t"+str(avg)+"\t"+str(var)+"\n") #this is the change part
   
   pro_list2=[]
   is_first=1
   num=1
   for item in pro_list:# pro_list = [[],[],[],[]]
       if is_first:
          is_first=0
          pre_first=item[0]  #ID 
          pre_second=item[1] #region
          pro_list2.append(item)
          continue 
       first=item[0]
       second=item[1] 
       if first!=pre_first or second!=pre_second :
          Avge_var_count_mature_2(pro_list2,files_num)
          num+=1
          pre_first=first 
          pre_second=second
          pro_list2=[]
          pro_list2.append(item)
       else:
          pro_list2.append(item)
   Avge_var_count_mature_2(pro_list2,files_num)
   print "already process "+str(num)
'''
def Process_Fold_mirBase(input_1,input_2,output_file):
   files = os.listdir(input_1+" "+input_2)
   output =output_file
   if os.path.exists("result"):
        pass
   else:
       os.makedirs("result")
   if os.path.exists(output):
       shutil.rmtree(output)
       os.makedirs(output)
   else:
       os.makedirs(output)
   count=0
   print 
   print "now is process"+output
   print 
   for file in files:
       count+=1
       preprocess.PreProcess(input_1,input_2+"/"+file,output+"/new"+file)#this is the change part
   print "the totle number of files is : "+str(count) 
   # now star to get the family information 
   process_accession_product.process_accession_product_mirBase("miRNA.dat",output)
   files = os.listdir("Family/"+output)
   fw=open("result/"+output_file+".txt",'w+')
   fw.write("miRNA_ID\tAverage\tVariance\tline_number\tread_counts\tmiRNA-family\n")#this is the change part
   pro_list=[]
   for file in files:
        is_first=1
        fr=open("Family/"+output+"/"+file,'r')
        for line in fr:
           if is_first:
              is_first=0
              continue
           list=line.strip().split("\t")
           if list[10]=="No MIR_family":
              continue
           pro_list.append(list)
        is_first=1
   shutil.rmtree(output_file)
   pro_list.sort(key=lambda x:(x[0],x[10])) #this is the change part
   pro_list2=[]
   is_first=1

   def Avge_var_count_mirBase(process_list,count):
       if process_list[0][10]!="No MIR_family":
          count_=0  # count_ is the line number
          sum=0.0
          read_counts=0
          for item in process_list:
              count_+=1
              sum=sum+float(item[6])
              read_counts=read_counts+int(item[2])
          avg = sum/count_
          var=0.0
          for item in process_list:
             var=var+(float(item[6])-avg)**2
          var=var/count_ 
          fw.write(process_list[0][0]+"\t"+str(avg)+"\t"+str(var)+"\t"+str(count_)+"\t"+str(read_counts)+"\t"+process_list[0][10]+"\n") #this is the change part

   count=0
   num=1
   for item in pro_list:
       if item[10]=="No MIR_family":
          continue
       #print item[0]+"\t"+item[10]
       if is_first:
          is_first=0
          pre_first=item[0]
          pre_second=item[10] # this is the change part
          pro_list2.append(item)
          continue 
       first=item[0]
       second=item[10] #this is the change part
       if first!=pre_first or second!=pre_second:
          Avge_var_count_mirBase(pro_list2,count)
          num+=1
          count=1
          pre_first=first
          pre_second=second
          pro_list2=[]
          pro_list2.append(item)
       else:
          count=count+1
          pro_list2.append(item)
   print "already process "+str(num)
'''
def mirBase_2_preprocess(pre_file):#pre_file is [[],[],[]]
#output is : "miRNA_ID \t family \t read_counts \t rates \n"
    pre_file.sort(key=lambda x:(x[10],x[0]))#sort the faimly first and the ID 
    is_first=1  
    pro_list=[] 
    sum_=0.0   # save the reads_counts with the same family and ID 
    new_list2=[]  # save the new data "miRNA_ID \t family \t read_counts"
    for data in pre_file:  
         if is_first:   
            is_first=0   
            pre_ID = data[0]    
            pre_family = data[10]
            pro_list.append(data)  
            continue  
         ID = data[0]   
         family = data[10]  
         if pre_ID != ID  or pre_family != family: 
               sum_=0
               for item2 in pro_list:  #pro_list is : [[],[],[]]
                   read_count2 = int(item2[2]) 
                   sum_=sum_+read_count2
               new_list2.append(item2[0]+"\t"+item2[10]+"\t"+str(sum_))
               pre_ID = ID   
               pre_family = family  
               pro_list=[]
               pro_list.append(data) 
         else :  
               pro_list.append(data)
    sum_=0
    for item2 in pro_list:  #pro_list is : [[],[],[]]
        read_count2 = int(item2[2])
        sum_=sum_+read_count2
    new_list2.append(item2[0]+"\t"+item2[10]+"\t"+str(sum_))
    is_first=1
    sum_=0.0 #same the same family's read_counts 
    pro_list=[]
    return_list=[]
    for item in new_list2:# new_list2 is : "miRNA_ID \t family \t read_counts"
        if is_first:
           is_first=0
           data=item.strip().split()
           pre_family=data[1]  #split the part according to the family 
           sum_=sum_+int(data[2])
           pro_list.append(item)
           continue
        data=item.strip().split()
        family=data[1]
        read_count=data[2]
        if pre_family!=family:
            for item2 in pro_list:# get each line's rates 
                data2=item2.strip().split("\t")
                read_count2=int(data2[2])
                rates = float(read_count2)/sum_
                return_list.append(data2[0]+"\t"+data2[1]+"\t"+data2[2]+"\t"+str(rates)+"\n")
                #output is : "miRNA_ID \t family \t read_counts \t rates \n"
            sum_ = int(read_count)
            pre_family=family
            pro_list=[]
            pro_list.append(item)
        else: 
            sum_=sum_+int(read_count)
            pro_list.append(item)
    for item2 in pro_list:# get each line's rates
        data2=item2.strip().split("\t")
        read_count2=int(data2[2])
        rates = float(read_count2)/sum_
        return_list.append(data2[0]+"\t"+data2[1]+"\t"+data2[2]+"\t"+str(rates)+"\n")
    return return_list

# as we will count the retes in the family level ,so we must get the family info first
#1) for each file we delete the data without star or mature 
#2) for each file we get get the family according to the MIMAT 
#3) for each file we delete the data with family 
#4) call the function to get the ID-family-read_counts-rates (sort family , and then ID , combine the same family and ID , get the rates for each line )
#5) put the same ID and family together 
def Process_Fold_mirBase_2(input_1,input_2,output_file):
   files = os.listdir(input_1+" "+input_2)
   output =output_file
   if os.path.exists("result"):
        pass
   else:
       os.makedirs("result")
   if os.path.exists(output):
       shutil.rmtree(output)
       os.makedirs(output)
   else:
       os.makedirs(output)
   count=0
   print 
   print "now is process"+output
   print 
   for file in files: # for each file we delete the line without star and mature and add one feature rates
       count+=1
       preprocess.PreProcess_mirBase_2(input_1,input_2+"/"+file,output+"/new"+file)#this is the change part
   print "the totle number of files is : "+str(count) 
   # now star to get the family information 
   #the input is line[5]+"rates"  , and then get the family info 
   process_accession_product.process_accession_product_mirBase("miRNA.dat",output)
   #the output is line[5]+rates+"\taccession\tproduct\tRF-family\tmiRNA-family"+"\n"
   files = os.listdir("Family/"+output)
   fw=open("result/"+output+".txt",'w+')
   fw.write("mirRNA_ID\tmiRNA-family\tread_counts\tfile_num\tAverage\tStandardDeviation\n")#this is the change part
   pro_list=[]# save all the file's data 
   files_num=0
   for file in files:
        is_first=1
        files_num+=1
        fr=open("Family/"+output+"/"+file,'r')
        pre_file=[]
        for line in fr:
           if is_first:
              is_first=0
              continue
           list=line.strip().split("\t")
           if list[10]=="No MIR_family":
              continue
           pre_file.append(list)
        pre_file= mirBase_2_preprocess(pre_file)  #call the other preprocess function 
		#output is : "miRNA_ID \t family \t read_counts \t rates \n"
        for item in pre_file:
             s = item.strip().split("\t")
             pro_list.append(s)
        is_first=1
   shutil.rmtree(output)
   pro_list.sort(key=lambda x:(x[1],x[0])) #pro_list is : [[],[],[]]

   def Avge_var_count_mirBase_2(process_list,files_num):#pro_list2 is : [[],[],[],[]]
   #[[miRNA_ID \t family \t read_counts \t rates ],[miRNA_ID \t family \t read_counts \t rates ]]
    #output is : "mirRNA_ID\tmiRNA-family\tread_counts\tfile_num\tAverage\tStandardDeviation\n"
       if process_list[0][1]!="No MIR_family":
          sum=0.0
          read_counts=0
          for item in process_list:
              sum=sum+float(item[3])
              read_counts=read_counts+int(item[2])
          avg = sum/files_num
          var=0.0
          for item in process_list:
             var=var+(float(item[3])-avg)**2
          var=math.sqrt(var/files_num) 
          fw.write(process_list[0][0]+"\t"+process_list[0][1]+"\t"+str(read_counts)+"\t"+str(files_num)+"\t"+str(avg)+"\t"+str(var)+"\n") #this is the change part

   num=1
   pro_list2=[]
   is_first=1
   for item in pro_list:
       if item[1]=="No MIR_family":
          continue
       if is_first:
          is_first=0
          pre_first=item[1]
          pre_second=item[0] 
          pro_list2.append(item)#pro_list2 is : [[],[],[],[]]
          continue 
       first=item[1]
       second=item[0] 
       if first!=pre_first or second != pre_second:
          Avge_var_count_mirBase_2(pro_list2,files_num)
          num+=1 
          pre_first=first
          pre_second=second
          pro_list2=[]
          pro_list2.append(item)
       else: 
          pro_list2.append(item)
   Avge_var_count_mirBase_2(pro_list2,files_num)
   print "already process "+str(num)

'''
def Process_Fold_RF(input_1,input_2,output_file):
   files = os.listdir(input_1+" "+input_2)
   output =output_file
   if os.path.exists("result"):
        pass
   else:
       os.makedirs("result")
   if os.path.exists(output):
       shutil.rmtree(output)
       os.makedirs(output)
   else:
       os.makedirs(output)
   count=0
   print 
   print "now is process"+output
   print 
   for file in files:
       count+=1
       preprocess.PreProcess(input_1,input_2+"/"+file,output+"/new"+file)#this is the change part
   print "the totle number of files is : "+str(count) 
   # now star to get the family information 
   process_accession_product.process_accession_product_RF("miRNA.dat",output)
   files = os.listdir("Family/"+output)
   fw=open("result/"+output_file+".txt",'w+')
   fw.write("miRNA_ID\tAverage\tVariance\tline_number\tread_counts\tRF-family\n")#this is the change part
   pro_list=[]
   for file in files:
        is_first=1
        fr=open("Family/"+output+"/"+file,'r')
        for line in fr:
           if is_first:
              is_first=0
              continue
           list=line.strip().split("\t")
           if list[10]=="No RF_family":
              continue
           pro_list.append(list)
        is_first=1
   shutil.rmtree(output_file)
   pro_list.sort(key=lambda x:(x[0],x[9])) #this is the change part
   pro_list2=[]
   is_first=1

   def Avge_var_count_RF(process_list,count):
       if process_list[0][9]!="No RF_family":
          count_=0  # count_ is the line number
          sum=0.0
          read_counts=0
          for item in process_list:
              count_+=1
              sum=sum+float(item[6])
              read_counts=read_counts+int(item[2])
          avg = sum/count_
          var=0.0
          for item in process_list:
             var=var+(float(item[6])-avg)**2
          var=var/count_ 
          fw.write(process_list[0][0]+"\t"+str(avg)+"\t"+str(var)+"\t"+str(count_)+"\t"+str(read_counts)+"\t"+process_list[0][9]+"\n") #this is the change part

   count=0
   num=1
   for item in pro_list:
       if item[10]=="No RF_family":
          continue
       #print item[0]+"\t"+item[9]
       if is_first:
          is_first=0
          pre_first=item[0]
          pre_second=item[9] # this is the change part
          pro_list2.append(item)
          continue 
       first=item[0]
       second=item[9] #this is the change part
       if first!=pre_first or second!=pre_second:
          Avge_var_count_RF(pro_list2,count)
          num+=1
          count=1
          pre_first=first
          pre_second=second
          pro_list2=[]
          pro_list2.append(item)
       else:
          count=count+1
          pro_list2.append(item)
   print "already process "+str(num)
'''
def RF_2_preprocess(pre_file):
    pre_file.sort(key=lambda x:(x[9],x[0]))
    is_first=1  
    pro_list=[] 
    sum_=0.0   
    new_list2=[]
    for data in pre_file:  
         if is_first:   
            is_first=0  
            pre_ID = data[0]    
            pre_family = data[9]
            pro_list.append(data)  
            continue   
         ID = data[0]   
         family = data[9]  
         if pre_ID != ID  or pre_family != family: 
               sum_=0
               for item2 in pro_list: 
                   read_count2 = int(item2[2]) 
                   sum_=sum_+read_count2
               new_list2.append(item2[0]+"\t"+item2[9]+"\t"+str(sum_)) 
               pre_ID = ID   
               pre_family = family  
               pro_list=[]
               pro_list.append(data) 
         else :  
               pro_list.append(data)
    sum_=0
    for item2 in pro_list:
        read_count2 = int(item2[2])
        sum_=sum_+read_count2
    new_list2.append(item2[0]+"\t"+item2[9]+"\t"+str(sum_))
    is_first=1
    sum_=0.0
    pro_list=[]
    return_list=[]
    for item in new_list2:
        if is_first:
           is_first=0
           data=item.strip().split() 
           pre_family=data[1]
           sum_=sum_+int(data[2])
           pro_list.append(item)
           continue
        data=item.strip().split() 
        family=data[1]
        read_count=data[2]
        if pre_family!=family:
            for item2 in pro_list:
                data2=item2.strip().split()
                read_count2=int(data2[2])
                rates = float(read_count2)/sum_
                return_list.append(data2[0]+"\t"+data2[1]+"\t"+data2[2]+"\t"+str(rates)+"\n")
                #output is : "miRNA_ID \t family \t read_counts \t rates \n"
            sum_ = int(read_count) 
            pre_family=family
            pro_list=[]
            pro_list.append(item)
        else: 
            sum_=sum_+int(read_count)
            pro_list.append(item)
    for item2 in pro_list:
        data2=item2.strip().split()
        read_count2=int(data2[2])
        rates = float(read_count2)/sum_
        return_list.append(data2[0]+"\t"+data2[1]+"\t"+data2[2]+"\t"+str(rates)+"\n")
    return return_list

def Process_Fold_RF_2(input_1,input_2,output_file):
   files = os.listdir(input_1+" "+input_2)
   output =output_file
   if os.path.exists("result"):
        pass
   else:
       os.makedirs("result")
   if os.path.exists(output):
       shutil.rmtree(output)
       os.makedirs(output)
   else:
       os.makedirs(output)
   count=0
   print 
   print "now is process"+output
   print 
   for file in files:
       count+=1
       preprocess.PreProcess_RF_2(input_1,input_2+"/"+file,output+"/new"+file)#this is the change part
   print "the totle number of files is : "+str(count) 
   # now star to get the family information 
   process_accession_product.process_accession_product_RF("miRNA.dat",output)
   files = os.listdir("Family/"+output)
   fw=open("result/"+output_file+".txt",'w+')
   fw.write("mirRNA_ID\tmiRF-family\tread_counts\tfile_num\tAverage\tStandardDeviation\n")#this is the change part
   pro_list=[]
   files_num=0
   for file in files:
        is_first=1
        files_num+=1
        fr=open("Family/"+output+"/"+file,'r')
        pre_file=[]
        for line in fr:
           if is_first:
              is_first=0
              continue
           list=line.strip().split("\t")
           if list[9]=="No RF_family":
              continue
           pre_file.append(list)
        pre_file= RF_2_preprocess(pre_file)  #call the other preprocess function 
        for item in pre_file:
             s = item.strip().split("\t")
             pro_list.append(s)
        is_first=1
   shutil.rmtree(output_file)
   pro_list.sort(key=lambda x:(x[1],x[0])) #this is the change part

   def Avge_var_count_RF_2(process_list,files_num):
   #input is : "mirRNA_ID\tmiRNA-family\tread_counts\tfile_num\tAverage\tStandardDeviation\n"
       if process_list[0][1]!="No RF_family": 
          sum=0.0
          read_counts=0
          for item in process_list: 
              sum=sum+float(item[3])
              read_counts=read_counts+int(item[2])
          avg = sum/files_num
          var=0.0
          for item in process_list:
             var=var+(float(item[3])-avg)**2
          var=math.sqrt(var/files_num)
          fw.write(process_list[0][0]+"\t"+process_list[0][1]+"\t"+str(read_counts)+"\t"+str(files_num)+"\t"+str(avg)+"\t"+str(var)+"\n") #this is the change part

   num=1
   pro_list2=[]
   is_first=1
   for item in pro_list:
       if item[1]=="No RF_family":
          continue 
       if is_first:
          is_first=0
          pre_first=item[1]
          pre_second=item[0] 
          pro_list2.append(item)
          continue 
       first=item[1]
       second=item[0]  
       if first!=pre_first or second != pre_second:
          Avge_var_count_RF_2(pro_list2,files_num)
          num+=1 
          pre_first=first
          pre_second=second
          pro_list2=[]
          pro_list2.append(item)
       else: 
          pro_list2.append(item)
   Avge_var_count_RF_2(pro_list2,files_num)
   print "already process "+str(num)




