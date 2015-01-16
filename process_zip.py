import sys 
import os 
import subprocess
import shutil
import process_fold
import process_accession_product
import process_ztest_pvalue
from optparse import OptionParser

# this part show the help of our codes 
usage = "usage: python %prog [options] arg1[-l zip] arg2[-l level]"
parser = OptionParser(usage=usage)
parser.add_option("-f","--file",dest="zip_files",help="the fold of the zip files")
parser.add_option("-l","--level",dest="calculation_level",default="misomiR",help="isomiR , mature , precursor , mirBase , RF")
(options,args)=parser.parse_args()

level=options.calculation_level
if len(sys.argv)!=5:
   parser.error("please input the correct number of the args")

if level!="isomiR" and level!="mature" and level!="precursor" and level!="mirBase" and level!="RF":
   parser.error("please input the correct -l level args : isomiR , mature , precursor ,mirBase , RF")

def main():
   files = os.listdir(options.zip_files)  # iter the zip fold 
   for file in files:
           index = file.find('.') # get the file name index
           child=subprocess.Popen(["unzip",options.zip_files+"/"+file,"-d",file[:index]])
           print "unzip the file:"+file
           child.wait()
           if level=="isomiR":
           #   pass
              if os.path.exists(file[:index]+"/"+file[:index]+" NT/"): # if the dir Breast NT is exists and then handle this dir
                 process_fold.Process_Fold_isomiR(file[:index]+"/"+file[:index],"NT/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",file[:index]+"_NT")
              if os.path.exists(file[:index]+"/"+file[:index]+" TN/"):
                 process_fold.Process_Fold_isomiR(file[:index]+"/"+file[:index],"TN/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",file[:index]+"_TN")
              shutil.rmtree(file[:index]) #delete the Breast NT and the Breast TN , this two folds store the new_files 
           if level=="mature":
        #     pass
              if os.path.exists(file[:index]+"/"+file[:index]+" NT/"):
                 process_fold.Process_Fold_mature_2(file[:index]+"/"+file[:index],"NT/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",file[:index]+"_NT")
              if os.path.exists(file[:index]+"/"+file[:index]+" TN/"):
                 process_fold.Process_Fold_mature_2(file[:index]+"/"+file[:index],"TN/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",file[:index]+"_TN")
              shutil.rmtree(file[:index])
           if level=="precursor":
         #   pass
              if os.path.exists(file[:index]+"/"+file[:index]+" NT/"):
                 process_fold.Process_Fold_precursor(file[:index]+"/"+file[:index],"NT/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",file[:index]+"_NT")
              if os.path.exists(file[:index]+"/"+file[:index]+" TN/"):
                 process_fold.Process_Fold_precursor(file[:index]+"/"+file[:index],"TN/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",file[:index]+"_TN")
              shutil.rmtree(file[:index])
           if level=="mirBase":
        #   pass
              if os.path.exists(file[:index]+"/"+file[:index]+" NT/"):
                 process_fold.Process_Fold_mirBase_2(file[:index]+"/"+file[:index],"NT/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",file[:index]+"_NT")
              if os.path.exists(file[:index]+"/"+file[:index]+" TN/"):
                 process_fold.Process_Fold_mirBase_2(file[:index]+"/"+file[:index],"TN/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",file[:index]+"_TN")
              shutil.rmtree(file[:index])
           if level =="RF":
              if os.path.exists(file[:index]+"/"+file[:index]+" NT/"):
                 process_fold.Process_Fold_RF_2(file[:index]+"/"+file[:index],"NT/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",file[:index]+"_NT")
              if os.path.exists(file[:index]+"/"+file[:index]+" TN/"):
                 process_fold.Process_Fold_RF_2(file[:index]+"/"+file[:index],"TN/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",file[:index]+"_TN")
              shutil.rmtree(file[:index])
           #as some history reason we have the mature-2 and mirBase-2 RF-2 selections but we may not use them
		   #if level =="mature-2":
           #  #    pass
           #   if os.path.exists(file[:index]+"/"+file[:index]+" NT/"):
           #      process_fold.Process_Fold_mature_2(file[:index]+"/"+file[:index],"NT/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",file[:index]+"_NT")
           #   if os.path.exists(file[:index]+"/"+file[:index]+" TN/"):
           #      process_fold.Process_Fold_mature_2(file[:index]+"/"+file[:index],"TN/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",file[:index]+"_TN")
           #   shutil.rmtree(file[:index])
           #if level =="mirBase-2":
           #   #    pass
           #   if os.path.exists(file[:index]+"/"+file[:index]+" NT/"):
           #      process_fold.Process_Fold_mirBase_2(file[:index]+"/"+file[:index],"NT/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",file[:index]+"_NT")
           #   if os.path.exists(file[:index]+"/"+file[:index]+" TN/"):
           #      process_fold.Process_Fold_mirBase_2(file[:index]+"/"+file[:index],"TN/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",file[:index]+"_TN")
           #   shutil.rmtree(file[:index])
           #if level =="RF-2":
           #   #    pass
           #   if os.path.exists(file[:index]+"/"+file[:index]+" NT/"):
           #      process_fold.Process_Fold_RF_2(file[:index]+"/"+file[:index],"NT/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",file[:index]+"_NT")
           #   if os.path.exists(file[:index]+"/"+file[:index]+" TN/"):
           #      process_fold.Process_Fold_RF_2(file[:index]+"/"+file[:index],"TN/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/isomiR/",file[:index]+"_TN")
           #   shutil.rmtree(file[:index])  

   if level=="isomiR":   # after the average and the standard deviation we will calculate the z-test and the p-value 
   #the input is  : "miRNA_ID\tiosform_coords\tmiRNA_region\tread_counts\files_num\tAverage\tStandardDeviation\n"  save in the result 
      process_accession_product.process_accession_product_isomiR("miRNA.dat","result")
   #the output is : "miRNA_ID\tiosform_coords\tmiRNA_region\tread_counts\files_num\tAverage\tStandardDeviation\taccession\tproduct\tRF-family\tmiRNA-family\n"
      process_ztest_pvalue.z_test_p_value_isomiR("Family") 
   #the output is : "miRNA_ID\tisoform_coords\tz-test\tp-value\tread_counts"
   elif level =="mature":
   #input is : "miRNA_ID\tmiRNA_region\tread_counts\tfiles_number\tAverage\tStardarddeviation\n"
       process_ztest_pvalue.z_test_p_value_mature_2("result")
	#the output is : "miRNA_ID\tmiRNA_region\tz-test\tp-value\tread_counts"
	#  process_accession_product.process_accession_product_mature("miRNA.dat","result")
    #  process_ztest_pvalue.z_test_p_value_mature("Family") 
   elif level =="precursor":
    #input is : miRNA_ID\tread_counts\tfiles_number\tAverage\tStardarddeviation\n
      process_ztest_pvalue.z_test_p_value_precursor("result")
	#output is : "miRNA_ID\tZ-test\tP-value\tread_counts\n"
   elif level=="mirBase":
   #input is : "mirRNA_ID\tmiRNA-family\tread_counts\tfile_num\tAverage\tStandardDeviation\n"
      process_ztest_pvalue.z_test_p_value_mirBase_2("result")
	#output is : "miRNA_ID\tfamily\tZ-test\tP-value\tread_counts\n"
   elif level=="RF":
   #input is : "mirRNA_ID\tmiRNA-family\tread_counts\tfile_num\tAverage\tStandardDeviation\n"
      process_ztest_pvalue.z_test_p_value_RF_2("result")
   #output is : "miRNA_ID\tfamily\tZ-test\tP-value\tread_counts\n"
   #elif level=="mature-2":
   #    process_ztest_pvalue.z_test_p_value_mature_2("result")
   #elif level=="mirBase-2":
   #    process_ztest_pvalue.z_test_p_value_mirBase_2("result")
   #elif level=="RF-2":
   #    process_ztest_pvalue.z_test_p_value_RF_2("result")

if __name__=="__main__":
     main()
