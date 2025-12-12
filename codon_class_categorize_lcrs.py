#!/usr/bin/python3

import re
from sys import argv



#splits a DNA sequence into a list of its respective codons
def Split_Into_Codons(DNAseq):

    codon_list = []
    for i in range(0, len(DNAseq)-2, 3):
        codon = DNAseq[i] + DNAseq[i+1] + DNAseq[i+2]
        codon_list.append(codon)
    return codon_list


#takes a DNA sequences and returns a dictionary of the number of each codon type (c1,c2,c3) in the sequence
def countcodonclasscategories(DNAseq):
    
    codon_list = Split_Into_Codons(DNAseq)
    d = {"c1":0, "c2":0, "c3":0}

    totalcod = 0  #for checking the numbers add up
    for cod in codon_list:
        class1 = (cod[0] == cod[1]) and (cod[0] == cod[2])     #defining a class1 codon
        class2 = (cod[0] == cod[1] and cod[0] != cod[2]) or (cod[0] == cod[2]\
                  and cod[0] != cod[1]) or (cod[1] == cod[2] and cod[0] != cod[1]) #defining a class 2 codon
        class3 = cod[0] != cod[1] and cod[0] != cod[2] and cod[1] != cod[2]  #defining a class3 codon

        if class1 == True:
            d["c1"] += 1
        elif class2 == True:
            d["c2"] += 1
        elif class3 == True:
            d["c3"] += 1
        totalcod += 1

    if d["c1"]+d["c2"]+d["c3"] != totalcod:    #to confirm that you are counting all of the codons
        raise ValueError("the codon types do not add up to the total number of codons.\
                c1: "+str(d["c1"])+"\tc2: "+str(d["c2"])+"\tc3: "+str(d["c3"])+ "\ttotalcod: "+str(totalcod))

    return d


#takes the dictionary of codon_classes and their counts and returns a list of the predominant codon class(es)
def predominantcodonclass(dictofcodclasscount):

    MAX = 0
    for key, val in dictofcodclasscount.items():
        if val > MAX:
            MAX = val
    predomcodclass = []
    for key, val in dictofcodclasscount.items():
        if val == MAX:
            predomcodclass.append(key)
    return predomcodclass


#def over50():
#    pass

#takes _Seq file and returns 4 files: LCRs where c1 codons are the most prevalent, LCRs where c2 codons are the most prevalent, LCRs where C3 codons are the most prevalent,
# and a file just with the number of LCRs in each category for ease... :)
def codonclass_categorizeLCRs(fileR, fileW=""):

    print("started running")
    totalLCRs = 0
    totalcods = 0
    c1_seq_info = {}
    c2_seq_info = {}
    c3_seq_info = {}
    totalcodoncount = {"c1":0, "c2":0, "c3":0}
    totalLCRcount   = {"c1":0, "c2":0, "c3":0} 
    with open(fileR, "r") as fh:
        f = fh.readlines()[4:]
        for line in range(len(f)):
            if line%7 == 0 :
                totalLCRs += 1
                ID = f[line][:-1]
                DNAseq = f[line+1][:-1]  #removes '\n'
               # Proseq = f[line+4][:-1]

                #call functions
                dictofcodonclasscounts = countcodonclasscategories(DNAseq)   #exp. {'c1': 9, 'c2': 44, 'c3': 3}
                listofpredomcodclass = predominantcodonclass(dictofcodonclasscounts) #exp. ['c2'] OR ['c2', 'c3']

                #get info from functions
                c1codnum = dictofcodonclasscounts["c1"]
                c2codnum = dictofcodonclasscounts["c2"]
                c3codnum = dictofcodonclasscounts["c3"]
                totalcodoncount["c1"] += c1codnum
                totalcodoncount["c2"] += c2codnum
                totalcodoncount["c3"] += c3codnum
                totalcods += c1codnum + c2codnum + c3codnum

                #store info in dictionary
                seqinfo = ["c1: "+str(c1codnum)+"\tc2: "+str(c2codnum)+"\tc3: "+str(c3codnum), DNAseq]
                #print(seqinfo)
                for codtype in listofpredomcodclass:
                    add = round(1/len(listofpredomcodclass),2)  #adds a fraction if 2 or 3 codon class types are in equal proportions
                    if codtype == "c1":
                        c1_seq_info[ID]= seqinfo
                        totalLCRcount["c1"]+= add  
                    elif codtype == "c2":
                        c2_seq_info[ID]= seqinfo
                        totalLCRcount["c2"]+= add
                    elif codtype == "c3":
                        c3_seq_info[ID]= seqinfo
                        totalLCRcount["c3"]+= add

    #print(c1_seq_info)
    #print(c2_seq_info)
    #print(c3_seq_info)
    #print(totalcodoncount["c1"], totalcodoncount["c2"], totalcodoncount["c3"])
    #print(totalLCRcount["c1"], totalLCRcount["c2"], totalLCRcount["c3"])

    #write to files
    #get some info to generate filenames
    pd = re.search("_pro_dna", fileR) 
    dp = re.search("_dna_pro", fileR) 
    if pd:
        direction = pd.group()
    elif dp:
        direction = dp.group()
    else:
        direction=""
    sc = re.search("Saccharomyces_cerevisiae", fileR) 
    hs = re.search("Homo_sapiens", fileR) 
    dm = re.search("Drosophila_melanogaster", fileR) 
    ce = re.search("Caenorhabditis_elegans", fileR) 
    at = re.search("Arabidopsis_thaliana", fileR) 
    if sc:
        org = "_sc" 
    elif hs:
        org = "_hs"
    elif dm:
        org = "_dm"
    elif ce:
        org = "_ce"
    elif at:
        org = "_at"
    else:
        org=""
    #can redirect files by specifying a fileW if needed
    c1filename = fileW+"predomc1"+org+direction
    c2filename = fileW+"predomc2"+org+direction
    c3filename = fileW+"predomc3"+org+direction
    countfilename = fileW+"codclassLCRcount"+org+direction
   # print(c1filename, c2filename, c3filename, countfilename)

    #gather the info into a string
    filenote = "#LCRs were categorized based on their predominant codon type\n#Info came from file: "+fileR+"\n"

    to_c1file = filenote+"#This file contains ID and codon counts for all LCRs encoded predominantly by class1 codons\n\n"
    for ID, seqInfo in c1_seq_info.items():
        to_c1file += ID + "\t" + seqInfo[0] + "\n" + seqInfo[1] + "\n\n"
    to_c2file = filenote+"#This file contains ID and codon counts for all LCRs encoded predominantly by class2 codons\n\n"
    for ID, seqInfo in c2_seq_info.items():
        to_c2file += ID + "\t" + seqInfo[0] + "\n" + seqInfo[1] + "\n\n"
    to_c3file = filenote+"#This file contains ID and codon counts for all LCRs encoded predominantly by class3 codons\n\n"
    for ID, seqInfo in c3_seq_info.items():
        to_c3file += ID + "\t" + seqInfo[0] + "\n" + seqInfo[1] + "\n\n"
    to_totcountfile = filenote+"#This file contains the number of LCRs which were encoded predominantly by each respective \
    codon class as well as the total number of codon types\ntotLCR\tLCRc1\tLCRc2\tLCRc3\tpLCRc1\tpLCRc2\tpLCRc3\ttotcod\ttotc1\ttotc2\ttotc3\tptotc1\tptotc2\tptotc3\n"
    to_totcountfile += str(totalLCRs)+"\t"  #total LCRs
    for c in totalLCRcount.values():    #number of LCRs with predominantly   c1,c2,c3
        to_totcountfile += str(c)+"\t"
    #sumto100 = 0
    for c in totalLCRcount.values():    #proportion of LCRs with predominanty c1,c2,c3
        lcrprop = str(round(c/totalLCRs,3))
        to_totcountfile += lcrprop + "\t"
        #sumto100 += float(lcrprop)  
    to_totcountfile += str(totalcods)+"\t"
    for c in totalcodoncount.values():  #total number of c1,c2,c3
        to_totcountfile += str(c)+"\t"
    for c in totalcodoncount.values():  #proportion of c1,c2,c3
        codprop = str(round(c/totalcods,3))
        to_totcountfile += codprop +"\t"
        #sumto100+= float(codprop)


    with open(c1filename, "w") as fh:
        fh.write(to_c1file)
    with open(c2filename, "w") as fh:
        fh.write(to_c2file)
    with open(c3filename, "w") as fh:
        fh.write(to_c3file)
    with open(countfilename, "w") as fh:
        fh.write(to_totcountfile)
    print("files are written")
    #print(sumto100)

#codonclass_categorizeLCRs("/home/JosLaptop/proj1/S.cer_Output/Pro_DNA/ran0722/Seq/Saccharomyces_cerevisiae151.92.2_pro_dna_Seq", fileW="/home/JosLaptop/test_codclasscatLCRs/")
#codonclass_categorizeLCRs("/home/JosLaptop/proj1/Human_Output/Pro_DNA/ran0729/Seq/Homo_sapiens151.92.2_pro_dna_Seq", fileW="/home/JosLaptop/test_codclasscatLCRs/")
#codonclass_categorizeLCRs("/home/JosLaptop/proj1/d.melanogaster_Output/Pro_DNA/ran0722/Seq/Drosophila_melanogaster151.92.2_pro_dna_Seq", fileW="/home/JosLaptop/test_codclasscatLCRs/")
#codonclass_categorizeLCRs("/home/JosLaptop/proj1/c.elegans_Output/Pro_DNA/ran0722/Seq/Caenorhabditis_elegans151.92.2_pro_dna_Seq", fileW="/home/JosLaptop/test_codclasscatLCRs/")
#codonclass_categorizeLCRs("/home/JosLaptop/proj1/A.thaliana_Output/Pro_DNA/ran0722/Seq/Arabidopsis_thaliana151.92.2_pro_dna_Seq", fileW="/home/JosLaptop/test_codclasscatLCRs/")

#codonclass_categorizeLCRs("/home/JosLaptop/proj1/S.cer_Output/DNA_Pro/ran0722/Seq/Saccharomyces_cerevisiae451.31.5_dna_pro_Seq", fileW="/home/JosLaptop/test_codclasscatLCRs/")
#codonclass_categorizeLCRs("/home/JosLaptop/proj1/Human_Output/DNA_Pro/ran0729/Seq/Homo_sapiens451.31.5_dna_pro_Seq", fileW="/home/JosLaptop/test_codclasscatLCRs/")
#codonclass_categorizeLCRs("/home/JosLaptop/proj1/d.melanogaster_Output/DNA_Pro/ran0722/Seq/Drosophila_melanogaster451.31.5_dna_pro_Seq", fileW="/home/JosLaptop/test_codclasscatLCRs/")
#codonclass_categorizeLCRs("/home/JosLaptop/proj1/c.elegans_Output/DNA_Pro/ran0722/Seq/Caenorhabditis_elegans451.31.5_dna_pro_Seq", fileW="/home/JosLaptop/test_codclasscatLCRs/")
#codonclass_categorizeLCRs("/home/JosLaptop/proj1/A.thaliana_Output/DNA_Pro/ran0722/Seq/Arabidopsis_thaliana451.31.5_dna_pro_Seq", fileW="/home/JosLaptop/test_codclasscatLCRs/")
if len(argv) == 3:
    codonclass_categorizeLCRs(argv[1], argv[2])
else:
    print("Seq file (from CombinedUsingSysCall_Codon3 to read from), fileW (mostly just for path redirection)")
        
