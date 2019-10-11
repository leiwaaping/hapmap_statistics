#############################################################################
############################## Header #######################################
# Date: Jul 3 2018
# File Name: nucleobase_count.py 
# Functions: 
#        1. data processing
# command: python2 *.py testCase_file testControl_file  output_file output_filename_of_matrix

#############################################################################
############################## Nucleobase_Count #############################
# Purpose: Process and analyze test data
# output format: SNPID  Main_Base   Other_Base  Case_A   Case_T   Ctr_A Ctr_T
# Tol_AA Tol_AT   Tol_TT   Tol_NN   Num_of_Sample
#
# NOTE: 
# 1, If case SNPID diff from control SNPID, only that line will print out
#    an error messgae. No other lines will be affected
# 2, All base combinations are ordered as MainMain MainMinor MinorMinor



# import package
import re
from itertools import izip
import sys

# take in file name 
caseFileName = sys.argv[1]
ctrFileName = sys.argv[2]
resultFileName = sys.argv[3]
resultFileName2 = sys.argv[4]

# open two files and create result file
fileCase = open(caseFileName,"r")
fileCtr = open(ctrFileName, "r")
fileResult = open(resultFileName, "w+")

# write header for the result file
fileResult.write("SNPID\tMain_Base\tMinor_Base\tCase_Main\tCase_Minor\t" +
"Ctr_Main\tCtr_Minor\tTol_#1_AA\tTol_#2_AT\tTol_#3_TT\tTol_NN\t"
+"Num_of_Sample\tCase_#1_AA\tCase_#2_AT\tCase_#3_TT\t" 
+"Contrl_#1_AA\tContrl_#2_AT\tContrl_#3_TT\n")

minorBF2 = [" "]
majorBF2 = [" "]
# initialize 
header = 0

# read in case file line by line
for lineCase, lineCtr in izip (fileCase, fileCtr):
  
   # initialize all variables for each line
   SNPID = ""
   mainBase = ""
   otherBase = ""
   Num_case_A = 0
   Num_case_T = 0
   Num_ctr_A =0
   Num_ctr_T = 0
   Num_Tol_AA = 0
   Num_Tol_AT = 0
   Num_Tol_TT = 0
   Num_Tol_NN = 0
   Num_sample= 0

   Num_case_AA = 0
   Num_case_AT = 0
   Num_case_TT = 0
   Num_ctr_AA = 0
   Num_ctr_AT = 0
   Num_ctr_TT = 0

  
   if header == 0:
      header = 1
   
   else:
      # seperate all lines
      splCase = lineCase.split()
      splCtr = lineCtr.split()
      if splCase[0] != splCtr[0] :
         print("Case ID: " + str(splCase[0]) + " different from control ID : " +
         str(splCtr[0]))
         minorBF2.append("DNE")
         majoeBF2.append("DNE")
      else:

         # get base name for the data
         SNPID = splCase[0]
         A = splCase[1][0]
         T = splCase[1][2]
         AA = str(A) + str(A)
         TT = str(T) + str(T)
      
         # loop through all lines
         for i in range(2,int(len(splCase))): 
            Num_sample += 1

            # calc case file
            if splCase[i] == AA:
               Num_Tol_AA += 1
               Num_case_A += 2
               Num_case_AA += 1
            elif splCase[i] == TT:
               Num_Tol_TT += 1
               Num_case_T += 2
               Num_case_TT += 1
            elif splCase[i] == "NN":
               Num_Tol_NN += 1
            else:
               Num_Tol_AT += 1
               Num_case_A += 1
               Num_case_T += 1
               Num_case_AT += 1

            # calc ctr file
         for i in range(2,int(len(splCtr))):
            Num_sample += 1

            if splCtr[i] == AA:
               Num_Tol_AA += 1
               Num_ctr_A += 2
               Num_ctr_AA += 1
            elif splCtr[i] == TT:
               Num_Tol_TT += 1
               Num_ctr_T += 2
               Num_ctr_TT += 1
            elif splCtr[i] == "NN":
               Num_Tol_NN += 1
            else:
               Num_Tol_AT += 1
               Num_ctr_A += 1
               Num_ctr_T += 1
               Num_ctr_AT += 1

         # calc main base
         if (Num_case_A + Num_ctr_A) > (Num_case_T + Num_case_T):
            mainBase = A
            otherBase = T
            minorBF2.append(T)
            majorBF2.append(A)
            # print message for Total, case, and ctr
            msgTol = (str(Num_Tol_AA) + "\t" + str(Num_Tol_AT)
            + "\t" + str(Num_Tol_TT))
            msgCase = (str(Num_case_AA) + "\t" + str(Num_case_AT)
            + "\t" + str(Num_case_TT))
            msgCtr = (str(Num_ctr_AA) + "\t" + str(Num_ctr_AT)
            + "\t" + str(Num_ctr_TT))
            msgCase1 = (str(Num_case_A) + "\t" + str(Num_case_T))
            msgCtr1 = (str(Num_ctr_A) + "\t" + str(Num_ctr_T))


         else: 
            mainBase = T
            otherBase = A
            minorBF2.append(A)
            majorBF2.append(T)
            # print message for Total, case, and ctr
            msgTol = (str(Num_Tol_TT) + "\t" + str(Num_Tol_AT)
            + "\t" + str(Num_Tol_AA))
            msgCase = (str(Num_case_TT) + "\t" + str(Num_case_AT)
            + "\t" + str(Num_case_AA))
            msgCtr = (str(Num_ctr_TT) + "\t" + str(Num_ctr_AT)
            + "\t" + str(Num_ctr_AA))
            msgCase1 = (str(Num_case_T) + "\t" + str(Num_case_A))
            msgCtr1 = (str(Num_ctr_T) + "\t" + str(Num_ctr_A))


         # write into file
         print ("SNPID = " + SNPID)
         fileResult.write(SNPID + "\t" + str(mainBase) + "\t" + str(otherBase)
         + "\t" + msgCase1 + "\t" + msgCtr1 + "\t" 
         + msgTol + "\t" +  str(Num_Tol_NN)+ "\t" + str(Num_sample) + "\t" + msgCase + "\t"
         + msgCtr + "\n") 

# close file
fileResult.close()


#############################################################################
############################## Nucleobase_Matrix #############################
# Purpose: Create a Matrix with number of minor base in the data
# output format: 
#        Sample_Name Case(1)/Control(0)   SNP1  SNP2  SNP3  SNP4  ...
#        sample1
#        sample2
#        sample3
#        ...

# NOTE: 
fileCase = open(caseFileName,"r")
fileCtr = open(ctrFileName, "r")
fileResult2 = open(resultFileName2, "w+")
fileCount = open(resultFileName, "r")




# open a new file
# create a list of strings
# one string per Sample
strList = ['\tCase(1)/Control(0)']

firstTime = 0
count =1 
# loop through all lines in both files
for lineCase in fileCase:
   # split the line
   splCase = lineCase.split()
         
   # first time:
   # create strings and append all sample names into all strings
   # append 1 to all strings
   if firstTime == 0 :
      for j in range (0,int(len(splCase))-2):
         newStr = str(splCase[j+2]) + "\t" + "1" + "\t"
         strList.append(newStr)
         
      firstTime = 1
   # print number of minor bases
   else :
      strList[0] = str(strList[0]) + "\t" + str(splCase[0])
      minorBF = minorBF2[count]
      MajorBF = majorBF2[count]
         
      for j in range (0,int(len(splCase))-2):
         mm = str(minorBF) + str(minorBF)
         mM = str(minorBF) + str(MajorBF)
         Mm = str(MajorBF) + str(minorBF)
         
         if str(splCase[j+2]) == str(mm) :
            strList[j+1] = str(strList[j+1]) + "2\t"
         elif str(splCase[j+2]) == str(mM) :
            strList[j+1] = str(strList[j+1]) + "1\t"
         elif str(splCase[j+2]) == str(Mm) :
            strList[j+1] = str(strList[j+1]) + "1\t"
         else:
            strList[j+1] = str(strList[j+1]) + "0\t"   
      count += 1

# get number of existing item in list
listSize = len(strList) - 1
# reset numbers
firstTime = 0
count =1 

# loop through all lines in control files
for lineCtr in fileCtr:
   # split the line
   splCtr = lineCtr.split()
         
   # first time:
   # create strings and append all sample names into all strings
   # append 1 to all strings
   if firstTime == 0 :
      for j in range (0,int(len(splCtr))-2):
         newStr = str(splCtr[j+2]) + "\t" + "0" + "\t"
         strList.append(newStr)
      firstTime = 1

   # access minor base   
   else :
      minorBF = minorBF2[count]
      MajorBF = majorBF2[count]
      
      if minorBF == "DNE":
         print("skip for diff SNPID")
      else:   
         # count number of minor bases   
         for j in range (0,int(len(splCtr))-2):
            mm = str(minorBF) + str(minorBF)
            mM = str(minorBF) + str(MajorBF)
            Mm = str(MajorBF) + str(minorBF)
            if str(splCtr[j+2]) == str(mm) :
               strList[j+1+listSize] = str(strList[j+1 +listSize]) + "2\t"
            elif str(splCtr[j+2]) == str(mM) :
               strList[j+1+listSize] = str(strList[j+1 + listSize]) + "1\t"
            elif str(splCtr[j+2]) == str(Mm) :
               strList[j+1+listSize] = str(strList[j+1+listSize]) + "1\t"
            elif str(splCtr[j+2]) == "NN" :
               strList[j+1+listSize] = str(strList[j+1+listSize]) + "3\t"
	    else:
               strList[j+1+listSize] = str(strList[j+1+listSize]) + "0\t"   
      count += 1
# write all items into the second file
for i in range (0, int(len(strList))):

   fileResult2.write(strList[i])
   fileResult2.write("\n")

# close the second result file
fileResult2.close()

























