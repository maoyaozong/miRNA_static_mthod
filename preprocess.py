import sys 

#preprocess isomiR once handle a file a time 
# 1)iter the file if the line hast not the start or the mature delete it 
# 2)calculate each line's rate , sum is the reads_counds of the same ID and region with MIMAT or star 
# 3)output is the " miRNA_ID[0] \t isoform_coords[1] \t read_count[2] \t miRNA_region[3] \t rates[4] " 
def PreProcess_isomiR(input_1,input_2,output):
    fr = open(input_1+" "+input_2,'r')  
    fw = open(output,"w+") # output is the Breast_NT/new_files this will be delete in the process_fold
    is_first=1 
    new_list=[]  #the new_list save the data without first line
    for line in fr: 
        if is_first:  
            is_first=0  
            fw.write("miRNA_ID\tisoform_coords\tread_count\tmiRNA_region\trates\n")
            continue
        data = line.strip().split()  
        if data[5].lower().find("star")==-1 and data[5].lower().find("mature")==-1: # no found the star or the mature  
            continue  
        new_list.append(line.strip()) 
    fr.close()
    is_first=1  
    pro_list=[] # the data with the same ID and the region 
    read_count_sum=0.0   
    for item in new_list:  
         if is_first:   
            is_first=0  
            data=item.strip().split() 
            pre_ID = data[0]    
            pre_region = data[5][0]
            read_count = data[2]
            read_count_sum = read_count_sum + int(read_count) 
            pro_list.append(item)  
            continue 
         data = item.strip().split() 
         ID = data[0]   # get the ID 
         miRNA_region = data[5][0]  # get the region's first char
         read_count = data[2]  
         if pre_ID != ID or  pre_region != miRNA_region : # ID and region have one is different 
               for item2 in pro_list:  
                   data2= item2.strip().split() # item2 is the line of each data
                   read_count2 = int(data2[2]) 
                   rate = float(read_count2)/read_count_sum   # get the rates
				   #" miRNA_ID[0] \t isoform_coords[1] \t read_count[2] \t miRNA_region[3] \t rates[4] "
                   fw.write(data2[0]+"\t"+data2[1]+"\t"+data2[2]+"\t"+data2[5]+"\t"+str(rate)+"\n")
               read_count_sum = int(read_count)   
               pre_ID = ID   
               pre_region = miRNA_region  
               pro_list=[]
               pro_list.append(item) 
         else :  
               read_count_sum =read_count_sum + int(read_count)
               pro_list.append(item)
    for item2 in pro_list:
        data2= item2.strip().split() # item2 is the line of each data
        read_count2 = int(data2[2])
        rate = float(read_count2)/read_count_sum   # get the rates
		#" miRNA_ID[0] \t isoform_coords[1] \t read_count[2] \t miRNA_region[3] \t rates[4] "
        fw.write(data2[0]+"\t"+data2[1]+"\t"+data2[2]+"\t"+data2[5]+"\t"+str(rate)+"\n")
'''
def PreProcess_mature_2_(input_1,input_2,output):
    fr = open(input_1+" "+input_2,'r')  
    fw = open(output,"w+") 
    is_first=1 
    new_list=[] 
    for line in fr: 
        if is_first:  
            is_first=0  
            fw.write("miRNA_ID\tmiRNA_region\tread_counts"+"\t"+"reads_rate\n")
            continue
        data = line.strip().split()  
        if data[5].lower().find("star")==-1 and data[5].lower().find("mature")==-1:  
            continue  
        new_list.append(line.strip()) 
    fr.close()
    is_first=1  
    pro_list=[] 
    sum_=0.0   
    new_list2=[]
    for item in new_list:  
         if is_first:   
            is_first=0  
            data=item.strip().split() 
            pre_ID = data[0]    
            pre_region = data[5]
            pro_list.append(item)  
            continue 
         data = item.strip().split() 
         ID = data[0]   
         miRNA_region = data[5]  
         if pre_ID != ID  or pre_region != miRNA_region: 
               sum_=0
               for item2 in pro_list:  
                   data2= item2.strip().split()
                   read_count2 = int(data2[2]) 
                   sum_=sum_+read_count2
               new_list2.append(data2[0]+"\t"+data2[5]+"\t"+str(sum_))
               #print new_list2
               pre_ID = ID   
               pre_region = miRNA_region  
               pro_list=[]
               pro_list.append(item) 
         else :  
               pro_list.append(item)
    is_first=1
    sum_=0.0
    pro_list=[]
    for item in new_list2:
        if is_first:
           is_first=0
           data=item.strip().split()
           pre_ID=data[0]
           pre_region=data[1][0]
           sum_=sum_+int(data[2])
           pro_list.append(item)
           continue
        data=item.strip().split()
        ID=data[0]
        miRNA_region=data[1][0]
        read_count=data[2]
        if pre_ID!=ID or pre_region!=miRNA_region:
            for item2 in pro_list:
                data2=item2.strip().split()
                read_count2=int(data2[2])
                rates = float(read_count2)/sum_
                fw.write(data2[0]+"\t"+data2[1]+"\t"+data2[2]+"\t"+str(rates)+"\n")
            sum_ = int(read_count)
            pre_ID = ID
            pre_region=miRNA_region
            pro_list=[]
            pro_list.append(item)
        else: 
            sum_=sum_+int(read_count)
            pro_list.append(item)
'''
#handle the mature level 
#1)delete the data without star and mature 
#2)sort the data according to the region , so the same MIMAT will come to together
#3)count the totle read_counts of the same MIMAT ----sum 
#4)combine the same ID and same MIMAT , get the read_counts , and then get the rates=read_counts/sum
#5)so the output is the miRNA_ID[0] \t miRNA_region[1] \t read_counts[2] \t rates[3] 
def PreProcess_mature_2(input_1,input_2,output):
#the input is : miRNA_ID\tisoform_coords\tread_count\treads_per_million_miRNA_mapped\tcross-mapped\tmiRNA_region
    fr = open(input_1+" "+input_2,'r')  
    fw = open(output,"w+") 
    is_first=1 
    new_list=[]  
    for line in fr: 
        if is_first:  
            is_first=0  
            fw.write("miRNA_ID\tmiRNA_region\tread_counts\treads_rate\n")
            continue
        data = line.strip().split("\t")  
        if data[5].lower().find("star")==-1 and data[5].lower().find("mature")==-1:  
            continue 
        list=line.strip().split("\t") 
        new_list.append(list) 
    fr.close() 
    new_list.sort(key=lambda x:(x[5],x[0]))#first sorted by the MIMAT , and then is the ID 
    is_first=1  
    pro_list=[] 
    read_count_sum=0.0   
    new_list2=[]
    for data in new_list: # we combine the same ID and same MIMAT get the new read_counts , new_list=[[],[],[]]
         if is_first:   
            is_first=0 
            pre_ID = data[0]    
            pre_region = data[5]
            pro_list.append(data)  
            continue   
         ID = data[0]   
         miRNA_region = data[5]  
         if pre_ID != ID  or pre_region != miRNA_region: 
               read_count_sum=0
               for item2 in pro_list: 
                   read_count2 = int(item2[2]) 
                   read_count_sum=read_count_sum+read_count2
               new_list2.append(item2[0]+"\t"+item2[5]+"\t"+str(read_count_sum)) #combine the same Id and MIMAT to get the new read_count
               pre_ID = ID   
               pre_region = miRNA_region  
               pro_list=[]
               pro_list.append(data) 
         else :  
               pro_list.append(data)
    read_count_sum=0
    for item2 in pro_list:
        read_count2 = int(item2[2])
        read_count_sum=read_count_sum+read_count2
    new_list2.append(item2[0]+"\t"+item2[5]+"\t"+str(read_count_sum)) #combine the same Id and MIMAT to get the new read_coun
    is_first=1
    sum_=0.0 # the sum_ is the MIMAT's sum
    pro_list=[]
    for item in new_list2:#the new_list2 is the "ID \t region \t read_counts " , new_list2 is the line
	#now is the one ID face one MIMAT , we just differ the region 
        if is_first:
           is_first=0
           data=item.strip().split() 
           pre_region=data[1]
           sum_=sum_+int(data[2])
           pro_list.append(item)
           continue
        data=item.strip().split() 
        miRNA_region=data[1]
        read_count=data[2]
        if pre_region!=miRNA_region:
            for item2 in pro_list:
                data2=item2.strip().split()
                read_count2=int(data2[2])
                rates = float(read_count2)/sum_
                fw.write(data2[0]+"\t"+data2[1]+"\t"+data2[2]+"\t"+str(rates)+"\n")
				#the output is : ID + region + read_counts + rates
            sum_ = int(read_count) 
            pre_region=miRNA_region
            pro_list=[]
            pro_list.append(item)
        else: 
            sum_=sum_+int(read_count)
            pro_list.append(item)
    for item2 in pro_list:
        data2=item2.strip().split()
        read_count2=int(data2[2])
        rates = float(read_count2)/sum_
        fw.write(data2[0]+"\t"+data2[1]+"\t"+data2[2]+"\t"+str(rates)+"\n")

def PreProcess_mirBase_2(input_1,input_2,output):
    fr = open(input_1+" "+input_2,'r')  
    fw = open(output,"w+")
    is_first=1
    for line in fr:
        if is_first:
          is_first=0
          fw.write(line.strip()+"\t"+"reads_rate\n")
          continue
        data = line.strip().split()  
        if data[5].lower().find("star")==-1 and data[5].lower().find("mature")==-1:  
            continue  
        fw.write(line.strip()+"\t"+"rates\n")
    fr.close()
    fw.close()
    '''  
    is_first=1 
    new_list=[] 
    for line in fr: 
        if is_first:  
            is_first=0  
            fw.write(line.strip()+"\t"+"reads_rate\n")
            continue
        data = line.strip().split()  
        if data[5].lower().find("star")==-1 and data[5].lower().find("mature")==-1:  
            continue  
        new_list.append(line.strip()) 
    fr.close()
    is_first=1  
    pro_list=[] 
    sum_=0.0   
    for item in new_list:  
         if is_first:   
            is_first=0  
            data=item.strip().split() 
            pre_ID = data[0]    
            pre_region = data[5][0]
            read_count = data[2]
            sum_ = sum_ + int(read_count) 
            pro_list.append(item)  
            continue 
         data = item.strip().split() 
         ID = data[0]   
         miRNA_region = data[5][0]  
         read_count = data[2]  
         if pre_ID != ID or  pre_region != miRNA_region : 
               for item2 in pro_list:  
                   data2= item2.strip().split()
                   read_count2 = int(data2[2]) 
                   rate = float(read_count2)/sum_   
                   fw.write(item2+"\t"+str(rate)+"\n")
               sum_ = int(read_count)   
               pre_ID = ID   
               pre_region = miRNA_region  
               pro_list=[]
               pro_list.append(item) 
         else :  
               sum_ =sum_ + int(read_count)
               pro_list.append(item)
     '''

def PreProcess_RF_2(input_1,input_2,output):
    fr = open(input_1+" "+input_2,'r')  
    fw = open(output,"w+")
    is_first=1
    for line in fr:
        if is_first:
          is_first=0
          fw.write(line.strip()+"\t"+"reads_rate\n")
          continue
        data = line.strip().split()  
        if data[5].lower().find("star")==-1 and data[5].lower().find("mature")==-1:  
            continue  
        fw.write(line.strip()+"\t"+"rates\n")
    fr.close()
    fw.close()

#1) count the first row one ID have how much rates ,combine same ID 
#2) calculate the sum of the reads ----reads_sum
#3) for each same ID get the rates 
#the output is "miRNA_ID \t read_counts \t rates"
def PreProcess_precursor(input_1,input_2,output):
    fr = open(input_1+" "+input_2,'r')  
    fw = open(output,"w+") 
    is_first=1 
    new_list=[] 
    for line in fr: 
        if is_first:  
            is_first=0  
            fw.write("miRNA_ID\tread_counts\treads_rate\n")
            continue
        data = line.strip().split()  
        if data[5].lower().find("star")==-1 and data[5].lower().find("mature")==-1:  
            continue  
        new_list.append(data) #all the data in the new_list , [[],[],[]]
    fr.close()
    is_first=1  
    pro_list=[] 
    sum_=0.0 #one part's sum 
    sum_all=0.0#whole files's sum
    new_list.sort(key=lambda x:(x[0])) #combine all the ID in one file 
    for item in new_list: 
        sum_all = sum_all + int(item[2])# get the whole file's reads sum 
    for data in new_list:
        if is_first:
           is_first=0 
           pre_ID=data[0]
           read_count = data[2]
           sum_ = sum_+int(read_count)
           pro_list.append(data)
           continue
        #data = item.strip().split()
        ID=data[0]
        read_count=data[2]
        if pre_ID!=ID:
           #ata2=pro_list[0].strip().split("\t")
           fw.write(pro_list[0][0]+"\t"+str(sum_)+"\t"+str(float(sum_)/sum_all)+"\n")
           sum_ =int(read_count)
           pre_ID =ID 
           pro_list=[]
           pro_list.append(data)
        else:
           sum_ = sum_+int(read_count)
           pro_list.append(data)
    fw.write(pro_list[0][0]+"\t"+str(sum_)+"\t"+str(float(sum_)/sum_all)+"\n")
