#!/usr/bin/python3

import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from datetime import datetime
from sys import argv
import copy
import subprocess
from math import log
import tempfile

#preference coefficient:
#    c1    , c2  , c3    , c1c2    , c1c3   , c2c3  , c1c2c3
#sc: 0.0511 , 1 , 0.6159 ,  0.5256 , 0.3335 , 0.808 ,  0.5557
#hs: 0.0401 , 1, 0.362 ,  0.5201,   0.2011,   0.681 ,  0.4674
#dm: 0.0104, 1 , 0.7067 , 0.5052 , 0.3586 ,   0.8534 , 0.5724
#ce: 0.037 , 1 , 0.4573 , 0.5185 , 0.2472 ,   0.7287 , 0.4981
#at: 0.0297 , 1 , 0.4197 , 0.5149 , 0.2247 ,  0.7099 , 0.4831


#takes the file of codon proportions and returns it as a list in order
def parsecodonpropfile(codonpropfile):
        
    with open(codonpropfile) as fh:
        file = fh.readlines()
        for line in file:
            if "prop" in line:
                list_props = line.split("\t")
                weights = []  #this will have codon proportions appended to it
                for i in range(1,len(list_props)-1):  #avoid 'props' and '' at start and end of list_props
                    weights.append(float(list_props[i]))
    return weights


def simulateCDS(num_iterations, codonbias, len_protein, slope, species):

    print("generating sequences")
    codons = ["TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG",\
              "TAT", "TAC", "TGT", "TGC", "TGG", "CTT", "CTC", "CTA", "CTG",\
              "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT",\
              "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC",\
              "ACA", "ACG", "AAT", "AAC","AAA", "AAG", "AGT", "AGC", "AGA",\
              "AGG","GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",\
               "GAT","GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"]


    allCDS = []
    j = 0
    while j < num_iterations: #this is how many proteins you want in your fasta and genbank files
        w2 = copy.copy(codonbias) ##reset it back to being the initial proportions for each new sequences time  
                                  #makes a non-aliased copy of codonbias
        codons100 = []            #this is where all the codons will get added
        cod1 = random.choices(codons, weights=w2, k=1)[0]  #choose initial codon with equal weights   
        cod_index = codons.index(cod1)   #this is so I change the correct proportion later
        codons100.append(cod1)     #add initial codon to start of sequence

        i = 1
        dna_sequence = ""
        while i < len_protein:      #sample 99 more codons to make a length of 100 (or however many)
            cod = random.choices(codons, weights=w2, k=1)[0]   #w2 is the array (61 long) of weights for each codon
            codons100.append(cod)                    #add the codon to the growing sequence
            cod_index = codons.index(cod)
            samplingbias = w2[cod_index]
            if codons100[i] == codons100[i-1]:   #compare the codon you just added to the previous one to see if they are the same
                samplingbias += slope       #increase probability of resampling linearly until the repeat is broken
            else:
                samplingbias = codonbias[cod_index]    #sampling weight goes back to all the same if the repeat is broken
            w2 = copy.copy(codonbias)   #reset it to the start and then change probability for only the one codon (if it is a repeat)
            w2[cod_index] = samplingbias  #and adjust that one proportion  (it may be the original weighting or greater depending)
            i += 1
#
        #string codons together
        for c in codons100:
            dna_sequence += c
        allCDS.append(dna_sequence)
        j += 1

    #write all coding sequences to fasta file and return file name
    to_file_DNA = ""
    ID = 0
    for sequence in allCDS:
        to_file_DNA += ">" + str(ID) + "\n"
        to_file_DNA += sequence + "TAA" + "\n\n"  #each coding sequence needs a stop codon because in the next program this is run through it removes the stop codon
        ID += 1
    #tfDNA2 = "/home/JosLaptop/temp/incslipsim_speciesspecific_codonclassbias."+str(num_iterations)+"."+str(slope)+"."+species+"."+"fasta"
    tfDNA2 = tempfile.NamedTemporaryFile(delete = False).name  #make new tempfile2
    with open(tfDNA2, "w") as fh:
        fh.write(to_file_DNA)
    return tfDNA2



def MakeTranslatedFasta(CDSfasta):
    print("making translated fasta")
    to_file_protein = ""
    for record in SeqIO.parse(CDSfasta, "fasta"):
        ID = record.id
        seq = record.seq
        to_file_protein += ">"+ID+"\n"+ str(Seq(seq).translate()[:-1])+"\n\n"
    filename = CDSfasta + ".translate"
    with open(filename, "w") as fh:
        fh.write(to_file_protein)

    return(filename)


def MakeGenBankFile(CDSfasta, num_iter, seed, len_protein, species, fileW):


    print("making genbank file")

    #concatenate all DNA seqs and make it the chromosome sequence
    chromosome_string = ""
    for rec in SeqIO.parse(CDSfasta, "fasta"):
        chromosome_string += rec.seq
    chromosome_object = Seq(chromosome_string)
    
    #create a record
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    name = "ImprovedCodClassSeqSim_"+species
    description = "ImprovedCodClassSeqSim " +species+"    Sequences were generated by generating random sequences using the increasing slippage\
    model and species specific rates of slippage, as well as species specific codon bias and protein length. Sequences were classified based on\
    predominant codon class type and were either accepted or rejected according to their preference coefficient.\nNum_iter = "+ str(num_iter)+\
    ", Seed = " + str(seed) + "   len_protein = " + str(len_protein) + "   species = "+species
    organism = "ImprovedCodClassSeqSim_"+species
    accessions = "NC_"+str(seed)+"_"+species
    genbankname = fileW+"improvedcodclassLCRsim_"+ str(num_iter)+"."+str(seed)+"."+species+".gb"
    record = SeqRecord(chromosome_object, id = dt_string, name = name,\
            description = description, annotations = {"molecule_type": "DNA", "accessions": [accessions], "organism": organism})
    
    for rec in SeqIO.parse(CDSfasta, "fasta"):
        #print("adding features")
        ID = rec.id
        DNAseq = rec.seq
        start = int(ID)*len(DNAseq)
        end = start + len(DNAseq)
        translation = Seq(DNAseq).translate()[:-1]
       
        #create feature annotations
        feature = SeqFeature(FeatureLocation(start=start, end=end), type = "CDS", strand=1,\
                  qualifiers= {"protein_id":[ID], "translation":[translation], "db_xref":["GeneID:"+ID]})
        record.features.append(feature)
    
    #save as genbank file
    with open(genbankname, "w") as fh:
        SeqIO.write(record, fh, "genbank")
    print("genbank file is written")



def findLCRs(tfName1, W, K1, K2):
    #run SEG on protein and pass all LCRs into new Fasta format file
    print("seg")
    SEG = subprocess.check_output(["segA", tfName1, W, K1, K2, "pro", "-l",]) #default alphabet size = 20
    tfName2 = tempfile.NamedTemporaryFile(delete = False).name  #make new tempfile2
    #tfName2 = "checkLCRsfromseg"
    with open(tfName2, "w") as fh:
        fh.write(SEG.decode())

    return tfName2

#can do a shorter version of the entropy calculator becuase we do not have to worry about ambiguous letters!!
def EntropyCalculator(proteinseq):

    #print("lcrs were the same length")
    count = {}  #will count all the different AAs in the sequence
    for let in proteinseq:
        if let in count:
            count[let] += 1
        else:
            count[let] = 1
    E = 0
    for freq in count.values():
        if freq != 0:
            Pi = freq/len(proteinseq)
            E += Pi*(log(Pi,2))
    E = -1*E
    return round(E,3)

def longestLCR(segfile):
    
    LCRs = {}
    for record in SeqIO.parse(segfile, "fasta"):
        index = record.id.index("(")
        ID = record.id[:index]
        Pos = record.id[index+1:-1].split("-")
        start = int(Pos[0])  #1 indexed
        end = int(Pos[1])    #1 indexed
        length = end+1-start   #+1 to count the AA it ends on too
        seq = record.seq
        if ID not in LCRs:  #if there is only 1 LCR in the protein sequence
            LCRs[ID] = [start, end, length, seq]
        elif length > LCRs[ID][2]:  #compare lengths and replace it if it is longer, or leave the previous one if that one is longer
            LCRs[ID] = [start, end, length, seq]
        elif length == LCRs[ID][2]:    #if the LCRs are equal length check which has the highest entropy
            newlcr = EntropyCalculator(seq)
            oldlcr = EntropyCalculator(LCRs[ID][3])
            if newlcr > oldlcr:       #if current one has a higher entropy replace it, otherwise oldlcr stays
                LCRs[ID] = [start, end, length, seq]

    #now that we have only the longest LCRs from each protein, will return a dictionary just of the LCRs ID:[start, end]
    return LCRs


#takes a fasta file of DNA coding sequences, and a dictionary containing the ID and the position that the LCR starts and stops at in the protein
#will determine the corresponding positions in the CDS which encode the LCR, and will return a dictionary of the ID, with the corresponding sequence
def maptoDNAseq(lcrbiasedfasta, dictoflcrpos):

    dnapos = {}
    for record in SeqIO.parse(lcrbiasedfasta, "fasta"):
        ID = record.id
        if ID in dictoflcrpos:
            dnaseq = record.seq
            protstart = dictoflcrpos[ID][0]  #1 indexed
            protend = dictoflcrpos[ID][1]
            dnastart = protstart*3-2     #1 indexed
            dnaend   = protend*3         #1 indexed
            lcrdnaseq = str(dnaseq[dnastart-1:dnaend])
            dnapos[ID] = lcrdnaseq
    return dnapos 

#splits a DNA sequence into a list of its respective codons
def Split_Into_Codons(DNAseq):

    codon_list = []
    for i in range(0, len(DNAseq)-2, 3):
        codon = DNAseq[i] + DNAseq[i+1] + DNAseq[i+2]
        codon_list.append(codon)
    return codon_list


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


#takes a dictionary of DNA sequences, counts the codon classes, feeds to predominantcodonclass to figure out 
#what the predominant codon type is --> returns dictionary of ID:"codonclasstype"
def countcodonclass(dnaseqs):

    print("classifying LCR according to predominant codon class")
    codtype = {}
    for ID, dnaseq in dnaseqs.items():
        d = {"c1":0, "c2":0, "c3":0}
        codon_list = Split_Into_Codons(dnaseq)

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

        predomcodclass = predominantcodonclass(d)   # a list of string(s) exp. ["c1"] or ["c1", "c2"]
        predomcodclassnew = ""      #if codons classes are equal in number string them together to make the combined type. exp. "c1c2"
        for ele in predomcodclass:
            predomcodclassnew += ele
        codtype[ID] = predomcodclassnew

    return codtype

def keeporrejectlcr(dictofcodclass, species):

    print("keep or reject")
    pref_coefficients = {"sc":{"c1":0.0511, "c2":1, "c3":0.6159, "c1c2":0.5256, "c1c3":0.3335, "c2c3":0.808 , "c1c2c3":0.5557},
                         "hs":{"c1":0.0401, "c2":1, "c3":0.362 , "c1c2":0.5201, "c1c3":0.2011, "c2c3":0.681 , "c1c2c3":0.4674},
                         "dm":{"c1":0.0104, "c2":1, "c3":0.7067, "c1c2":0.5052, "c1c3":0.3586, "c2c3":0.8534, "c1c2c3":0.5724},
                         "ce":{"c1":0.037 , "c2":1, "c3":0.4573, "c1c2":0.5185, "c1c3":0.2472, "c2c3":0.7287, "c1c2c3":0.4981},
                         "at":{"c1":0.0297, "c2":1, "c3":0.4197, "c1c2":0.5149, "c1c3":0.2247, "c2c3":0.7099, "c1c2c3":0.4831},
                         }

    keeporreject = {}
    for ID, codclass in dictofcodclass.items():
        randnum = random.random()    #random float 0 <= x < 1.0
        if randnum < pref_coefficients[species][codclass]:   #if its below the preference coefficient keep the sequence
            keeporreject[ID] = "keep"
        else:
            keeporreject[ID] = "reject"

    return keeporreject

#takes a fasta of all the sequences, as well as a file from the output of seg with only the sequences that contain LCRs
#determines which sequence IDs do not contain LCRs, and returns the ID as well as the sequences wihtout LCRs
def findseqswithoutlcrs(dnaseqfasta, segfilefasta):

    seqswithlcrs = {}
    seqswithoutlcrs = {}
    for record in SeqIO.parse(segfilefasta, "fasta"):
        ID = record.id
        endindex = ID.index("(")
        ID = ID[:endindex]
        seqswithlcrs[ID]="lcr"   #just using a dictionary cuz dictionary lookup is faster than a list
    for record in SeqIO.parse(dnaseqfasta, "fasta"):
        ID = record.id
        if ID not in seqswithlcrs:
            seqswithoutlcrs[ID]=str(record.seq)
    return seqswithoutlcrs

#this takes a file of secies specific codon bias proportions, and species type
# will simulate LCR sequences using increasing slippage model that also takes into account the codonclass proportions (class 1,2 or 3)
# codon class proportions are based off of precalculated preference coefficient
# returns fasta file of coding sequences, fasta file of corresponding protein sequences, genbank file
def CallonFunctions(fileW="", codonpropfile="", species="", existingfasta="", num_iterations = 100000, len_protein=100, slope=0.3, seed=10, W="15", K1="1.9", K2="2.2"):

    random.seed(seed)
    print("started running")
    if codonpropfile != "":
        codonbias = parsecodonpropfile(codonpropfile)   #this will give a list of weights in the same order as the 'codons' global list
    else:
        codonbias = ""

    #species can be 'hs' (humans), 'sc' (yeast), 'dm' (drosophila), 'at' (arabidopsis), 'ce' (c.elegans)
    if species == "hs":
      #  len_protein = 317
       # slope = 0.271       #trial and error (must test to figure out) 
        slope = 0.218
        len_protein = 584
    elif species == "sc":
      #  len_protein = 467
       # slope = 0.391
        slope = 0.384
        len_protein = 488
    elif species == "dm":
      #  len_protein = 500
       # slope = 0.336
        slope = 0.328
        len_protein = 533
    elif species == "ce":
     #   len_protein = 435
      #  slope = 0.341
        slope = 0.349
        len_protein = 409
    elif species == "at":
     #   len_protein = 410
        slope = 0.339
        len_protein = 406
     #   slope = 0.337

    if existingfasta != "":
        lcrbiasedfasta = existingfasta

    seqstokeep = []
    num_iterations2 = copy.copy(num_iterations)  #make another one that you can decrease the number of sequences generated to you don't have to generate more than 100 000 and then remove some to make 100 000
    while len(seqstokeep) < num_iterations:   #until you generate 100 000 LCRs with biological codon class proportions
        num_iterations2 = num_iterations - len(seqstokeep)
        print("kept sequences: ", len(seqstokeep))
        print("num_iter2 ", num_iterations2)
        lcrbiasedfasta = simulateCDS(num_iterations2, codonbias, len_protein, slope, species)
        translatedfasta = MakeTranslatedFasta(lcrbiasedfasta)   #translate these sequences into protein sequences
        segfile = findLCRs(translatedfasta, W, K1, K2)   #run translatedfasta through seg to ID lcrs
        dictlongestlcrs = longestLCR(segfile)        #dictionary in the form [ID]:[start, end, length, seq]   #for all the LCRs
        dictdnaseq = maptoDNAseq(lcrbiasedfasta, dictlongestlcrs)   #dictionary [ID]:"lcrdnaseq"
        dictcodclass = countcodonclass(dictdnaseq)   #dictionary [ID]:"c1", [ID]:"c2c3"
        dictkeeporreject = keeporrejectlcr(dictcodclass, species)    #dictionary in the form:  [ID]:"reject", ID:"keep"


        #### a quick test ####
        #for key in dictkeeporreject.keys():       #remove this test later
        #    prolcrseq = dictlongestlcrs[key][3]
        #    dnalcrseq = dictdnaseq[key]
        #    if prolcrseq != Seq(dnalcrseq).translate().lower():
        #        raise ValueError("The translated DNA sequence and the protein sequence do not correspond. There is a bug in your code.")
        ########################

        #find the sequences which do not contain LCRs and add them to the list as well.
        seqswithoutlcrs = findseqswithoutlcrs(lcrbiasedfasta, segfile)    #dictionary in the form {[ID]: "ACGTTTGCTCA", [ID]: "CCGTTGTTGCAAAAA"}
        for seq in seqswithoutlcrs.values():
            seqstokeep.append(seq)
        # add the sequences with LCRs which we are keeping to the list
        for record in SeqIO.parse(lcrbiasedfasta, "fasta"):
            ID = record.id
            dnaseq = record.seq
            if ID in dictkeeporreject and dictkeeporreject[ID]=="keep":
                seqstokeep.append(str(dnaseq))

    #seqstokeep = seqstokeep[:num_iterations]


    #make new fasta based off of all the LCRs in seqstokeep
    to_fasta = ""
    ID = 0
    for seq in seqstokeep:
        to_fasta += ">"+str(ID)+"\n"+seq+"\n\n"
        ID += 1
    filename = fileW + "incslipsim_speciesspecific_improvedcodclassbias."+str(num_iterations)+"."+str(seed)+"."+species+"."+"fasta"
    with open(filename, "w") as fh:
        fh.write(to_fasta)
    #maketranslatedfasta
    MakeTranslatedFasta(filename)
    #makegenbankfile
    MakeGenBankFile(filename, num_iterations, seed, len_protein, species, fileW) #codingseqfasta

    #should probably remove the intermediate fasta files
    subprocess.call(["rm", lcrbiasedfasta, translatedfasta, segfile])  #remove tempfile
    print("finished")



if len(argv) == 0:
    print("codon proportion file, species (exp. 'sc'), fileName (to write to)")
else:
    CallonFunctions(codonpropfile=argv[1], species=argv[2], fileW=argv[3])
#CallonFunctions(codonpropfile="/home/JosLaptop/proj1/nt_cod_aa_proportions/genomewidecounts/Saccharomyces_cerevisiae_codon_count", species="sc", fileW = "/home/JosLaptop/improvedcodclass_newprotlengths_files/")
#CallonFunctions(codonpropfile="/home/JosLaptop/proj1/nt_cod_aa_proportions/genomewidecounts/Homo_sapiens_codon_count"            , species="hs", fileW = "/home/JosLaptop/improvedcodclass_newprotlengths_files/")
#CallonFunctions(codonpropfile="/home/JosLaptop/proj1/nt_cod_aa_proportions/genomewidecounts/Drosophila_melanogaster_codon_count" , species="dm", fileW = "/home/JosLaptop/improvedcodclass_newprotlengths_files/")
#CallonFunctions(codonpropfile="/home/JosLaptop/proj1/nt_cod_aa_proportions/genomewidecounts/Caenorhabditis_elegans_codon_count"  , species="ce", fileW = "/home/JosLaptop/improvedcodclass_newprotlengths_files/")
#CallonFunctions(codonpropfile="/home/JosLaptop/proj1/nt_cod_aa_proportions/genomewidecounts/Arabidopsis_thaliana_codon_count"    , species="at", fileW = "/home/JosLaptop/improvedcodclass_newprotlengths_files/")
#

