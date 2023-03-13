#!/usr/bin/python3

import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from datetime import datetime
from sys import argv
import copy

#global variable of list of codons
codons = ["TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG",\
          "TAT", "TAC", "TGT", "TGC", "TGG", "CTT", "CTC", "CTA", "CTG",\
          "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT",\
          "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC",\
          "ACA", "ACG", "AAT", "AAC","AAA", "AAG", "AGT", "AGC", "AGA",\
          "AGG","GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",\
           "GAT","GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"]


## parses a file of codon proportions from nt_pro_cod_genomebias
#returns list of the proportions 'weights' in the same order that corresponds to 'codons'
def parsecodonpropfile(codonpropfile):
        
    with open(codonpropfile) as fh:
        file = fh.readlines()
        for line in file:
            if "prop" in line:
                list_props = line.split("\t")
                weights = []  #this will have the nt or codon proportions appended to it
                for i in range(1,len(list_props)-1):  #avoid 'props' and '' at start and end of list_props
                    weights.append(float(list_props[i]))
    return weights



#this code generates a proteome of 100 000 proteins or genes
# by randomly selecting codons and stringing them together to make a protein 100 AAs long
#prints sequence in Fasta format to random file
def RandomSeqSimulation(num_iterations, seed, fileW):

    print("getting random sequences")
    random.seed(seed)
    allCDS = []
    i = 0
    while i < num_iterations:
        j = 0
        DNASeq = ""
        while j < 100:
            DNASeq += random.choice(codons)
            j += 1
        allCDS.append(DNASeq)
        i += 1
    to_file_DNA = ""
    ID = 0
    for sequence in allCDS:
        to_file_DNA += ">" + str(ID) + "\n"
        to_file_DNA += sequence + "TAA" + "\n\n"
        ID += 1
    
    tfDNA1 = fileW+"RandomSeqSimulation.fasta"+ "." + str(num_iterations)+"." + str(seed)
    with open(tfDNA1, "w") as fh:
        fh.write(to_file_DNA)

    return tfDNA1

#this is the second part to the simulation.
#100 codon long sequences are created with a bias towards creating LCRs 
#(more likely to sample the same codon you are already on)
#p is the probability of resampling the same codon
def UnequalCodonSampling(num_iterations, seed, p, len_protein, species, codonbias, codontypefactor, fileW):

    print("getting biased sequences")
    random.seed(seed)

    if codonbias == "":
        filenote = "_initpropsequal"
        codonbias = []
        for prop in range(61):
            codonbias.append(1)  #start it off as equal 
    else:
        filenote = "_initbioprops_"+species
        for i in range(len(codonbias)):
            codonbias[i] = codonbias[i]
    if (codontypefactor[0] == 1) and (codontypefactor[1] == 1) and (codontypefactor[2] == 1):
        codbiasnote = ""
    else:
        codbiasnote = "_codonclassbias"



    allCDS = []
    j = 0
    while j < num_iterations:  #this is how many proteins you want in your genbank file
        w2 = copy.copy(codonbias)
        codons100 = ""
        cod = random.choices(codons, weights=w2, k=1)[0]  #choose initial codon with equal weights or weights of initial codon bias
        cod_index = codons.index(cod)
        codons100 += cod

        #establish what constitutes a class1 versus class2 versus class3 codon
        class1 = (cod[0] == cod[1]) and (cod[0] == cod[2])     #defining a class1 codon
        class2 = (cod[0] == cod[1] and cod[0] != cod[2]) or (cod[0] == cod[2]\
                 and cod[0] != cod[1]) or (cod[1] == cod[2] and cod[0] != cod[1]) #defining a class 2 codon
        class3 = cod[0] != cod[1] and cod[0] != cod[2] and cod[1] != cod[2]  #defining a class3 codon
        if class1 == True:
            codclassbias = codontypefactor[0]
        elif class2 == True:
            codclassbias = codontypefactor[1]
        elif class3 == True:
            codclassbias = codontypefactor[2]

        w2[cod_index] = p*codclassbias       #adjust the proportion of the codon you just chose   #multiply by the bias factor depending on the codon type
        
        i = 0
        while i < (len_protein-1):       #choose 99 more codons with different weights depending on the codon you are currently on
            cod = random.choices(codons, weights=w2, k=1)[0]  #w2 is the array (61 long) of weights for each codon
            cod_index = codons.index(cod)
            codons100 += cod
            w2 = copy.copy(codonbias)  #reset each time

            #establish what constitutes a class1 versus class2 versus class3 codon
            class1 = (cod[0] == cod[1]) and (cod[0] == cod[2])     #defining a class1 codon
            class2 = (cod[0] == cod[1] and cod[0] != cod[2]) or (cod[0] == cod[2]\
                     and cod[0] != cod[1]) or (cod[1] == cod[2] and cod[0] != cod[1]) #defining a class 2 codon
            class3 = cod[0] != cod[1] and cod[0] != cod[2] and cod[1] != cod[2]  #defining a class3 codon
            if class1 == True:
                codclassbias = codontypefactor[0]
            elif class2 == True:
                codclassbias = codontypefactor[1]
            elif class3 == True:
                codclassbias = codontypefactor[2]

            w2[cod_index] = p*codclassbias  #increase probability of codon you just chose
            i += 1
        allCDS.append(codons100)
        #print(j)
        j += 1

    #write all coding sequences to fasta file and return file name, seed and p value
    to_file_DNA = ""
    ID = 0
    for sequence in allCDS:
        to_file_DNA += ">" + str(ID) + "\n"
        to_file_DNA += sequence + "TAA" + "\n\n"
        ID += 1
    tfDNA2 = fileW+"UnequalCodonSamplingSimulation.fasta"+ "."+str(num_iterations)+"."+str(seed)+"."+str(p)+filenote+codbiasnote
    with open(tfDNA2, "w") as fh:
        fh.write(to_file_DNA)

    return tfDNA2

#this is the third part of the simulation.
#DNA is hammered with synonomous substitutions (only at 3rd codon position)
def ForceSynonomousMutations(Unequal_Sampling_fasta, num_mutations, seed, species, fileW):
   
    print("mutating seq")
    random.seed(seed)
    nuc = ["A", "C", "G", "T"]   #a list of nt to randomly sample from
    #dictionary of a codon and its corresponding synonomous codons
    synon_codons = {"TTT":["TTC"], "TTC":["TTT"], "TTA":["TTG","CTT","CTC", "CTA", "CTG"],\
                    "TTG":["TTA","CTT","CTC","CTA","CTG"], "CTT":["TTA","TTG", "CTC", "CTA", "CTG"],\
                    "CTC":["TTA","TTG","CTT","CTA","CTG"], "CTA":["TTA","TTG","CTT","CTC","CTG"],\
                    "CTG":["TTA","TTG","CTT","CTC","CTA"], "ATT":["ATC","ATA"], "ATC":["ATT","ATA"], "ATA":["ATT","ATC"], "ATG":[],\
                    "GTT":["GTC", "GTA", "GTG"], "GTC":["GTT","GTA", "GTG"], "GTA":["GTT","GTC","GTG"],\
                    "GTG":["GTC", "GTA", "GTT"], "TCT":["TCC","TCA","TCG", "AGT", "AGC"], "TCC":["TCT","TCA","TCG", "AGT", "AGC"],\
                    "TCA":["TCT","TCC","TCG", "AGT", "AGC"], "TCG":["TCT","TCC","TCA", "AGT", "AGC"], "CCT":["CCC","CCA","CCG"],\
                    "CCC":["CCT","CCA","CCG"], "CCA":["CCT","CCC","CCG"], "CCG":["CCT","CCC","CCA"],\
                    "ACT":["ACC","ACA","ACG"], "ACC":["ACT","ACA","ACG"], "ACA":["ACT","ACC","ACG"],\
                    "ACG":["ACT","ACC","ACA"], "GCT":["GCC","GCA","GCG"], "GCC":["GCT", "GCA", "GCG"],\
                    "GCA":["GCT", "GCC", "GCG"], "GCG":["GCT","GCC", "GCA"], "TAT":["TAC"], "TAC":["TAT"],\
                    "CAT":["CAC"], "CAC":["CAT"], "CAA":["CAG"], "CAG":["CAA"], "AAT":["AAC"], "AAC":["AAT"],\
                    "AAA":["AAG"], "AAG":["AAA"], "GAT":["GAC"], "GAC":["GAT"], "GAA":["GAG"], "GAG":["GAA"],\
                    "TGT":["TGC"], "TGC":["TGT"], "TGG":[], "CGT":["CGC", "CGA", "CGG", "AGA", "AGG"],\
                    "CGC":["CGT", "CGA", "CGG", "AGA", "AGG"], "CGA":["CGT", "CGC", "CGG", "AGA", "AGG"],\
                    "CGG":["CGT", "CGC", "CGA", "AGA", "AGG"], "AGT":["AGC", "TCT","TCC", "TCA", "TCG"],\
                    "AGC":["AGT", "TCT","TCC", "TCA", "TCG"], "AGA":["CGT", "CGC", "CGA", "CGG", "AGG"],\
                    "AGG":["CGT", "CGC", "CGA", "CGG", "AGA"], "GGT":["GGC", "GGA", "GGG"],\
                    "GGC":["GGT", "GGA", "GGG"], "GGA":["GGT", "GGC", "GGG"], "GGG":["GGT", "GGC", "GGA"]}

    to_new_fasta = ""
    count = 0
    for record in SeqIO.parse(Unequal_Sampling_fasta, "fasta"):   #takes from a previously created fasta of sequences biased towards LCR formation
        ID = record.id
        seq = record.seq[:-3]      #-3 so the TAA stop is removed (don't want to mutate)

        indv_codons = []
        #pick out indvl codons from sequence (str --> list)
        cod = ""
        for i in range(0, len(seq), 3):                   
            cod += seq[i] + seq[i+1] + seq[i+2]
            indv_codons.append(cod)        #add codons to list (of 100)
            cod = ""
        #now do while loop to randomly mutate a random codon (1--> numberofcodons)) at a nt position 1 or 3, with either: A, C, G, or T
        mutations = 0
        while mutations < num_mutations:
            nt = random.choice(nuc)               #mutate it with a random nt
            codon_num = random.randint(0, len(indv_codons)-1)    #len(seq)-1 because we will use this number for indexing, and 0 indexing is a thing
            codon_from_seq = indv_codons[codon_num]    #a random codon from the sequence
            index_codon = random.randrange(0,3,2)   #number will be only 0 or two (to mutate only first or third nt, as the second nt will never return a synonomous codon with 1 mutation)

            #now peice the codon back together, with the new mutation
            if index_codon == 0:
                mutant_codon = nt + codon_from_seq[1:]
            elif index_codon == 2:
                mutant_codon = codon_from_seq[0:2] + nt

            #if a synonomous mutation was made, replace the original codon with the mutant codon 
            if mutant_codon in synon_codons[codon_from_seq]:  
                indv_codons[codon_num] = mutant_codon
                mutations += 1
        count += 1
        #print(count)
        #now put mutated sequence back into string format and write new fasta
        mutated_seq = ""
        for c in indv_codons:
            mutated_seq += c
        mutated_seq += "TAA"
        to_new_fasta += ">" + ID + "\n" + mutated_seq + "\n\n"

    if species != "":
        species = "_"+species

    tfDNA3 = fileW+"MutatedLCRSimulation.fasta"+"."+str(num_mutations)+"."+str(seed)+species
    with open(tfDNA3, "w") as fh:
        fh.write(to_new_fasta)

    return tfDNA3

#ForceSynonomousMutations("randomfastafiles/biased/UnequalCodonSamplingSimulation.fasta.100000.10.27.8", 1000, 10, 1)


#this just makes another fasta file with corresponding protein sequences
#handy if you just want to run it on SEG and not the entire entropy program
def MakeTranslatedFasta(CDSfasta):
    to_file_protein = ""
    for record in SeqIO.parse(CDSfasta, "fasta"):
        ID = record.id
        seq = record.seq
        to_file_protein += ">"+ID+"\n"+ str(Seq(seq).translate()[:-1])+"\n\n"
    filename = CDSfasta + ".translate"
    with open(filename, "w") as fh:
        fh.write(to_file_protein)

#now generate genbank file from these random DNA sequences
#simulation type indicates whether the sequences are coming from truly random sequences or biased sequences
#CDSfasta is the fasta file that is converted to a genbank file
#CDSfasta2 is only for mutated type (to show which file the sequences that were mutated came from)
#simulation type is random, biased, or mutate
#num_iter is the number of sequences in the fasta file (argument is used to name file)
#seed (argument is used to name gb file)
#p (argument is used to name gb file)
def MakeGenBankFile(CDSfasta, simulation_type, num_iter, seed, p, num_mutations, len_protein, species, codonclassbias, codontypefactor, fileW, CDSfasta2=""):



    print("making genbank file")

    if species != "":
        filenote = "_initbioprops_"+species
    else:
        filenote = "_initequalprops"
   # if codonclassbias == "True":
   #     codbiasnote = "\tcodon classes 1, 2, and 3 were multiplied by the factors "+str(codontypefactor[0])+", "+str(codontypefactor[1])+", and "+str(codontypefactor[2])+", respectively"
   #     codbiasname = "_codonclassbias"
   # else:
   #     codbiasnote = ""
   #     codbiasname = ""



    #concatenate all DNA seqs and make it the chromosome sequence
    chromosome_string = ""
    for rec in SeqIO.parse(CDSfasta, "fasta"):
        chromosome_string += rec.seq
    chromosome_object = Seq(chromosome_string)
    
    #create a record
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    if simulation_type == "random":
        name = "Random_Sequences_Simulation"
        description = "Randomly Simulated genbank file containing 100 codon long\
        DNA sequences and their amino acid sequences. Num_iter = "+ str(num_iter) + ", Seed = " + str(seed)
        organism = "Random Sequences"
        accessions = "NC_"+str(seed)
        genbankname = fileW+"RandomSeqSimulation.genbank."+ str(num_iter)+"."+str(seed)
    elif simulation_type == "biased":
        name = "Biased_Sequences_Simulation"
        description = "Unequal Sampling of codons given a certain probability to\
        simulate random sequences with a bias towards LCR generation. Num_iter = "+ str(num_iter)+", Seed = " + str(seed) + ", p = " + str(p)+"  "+filenote 
        organism = "Biased Sequences"
        accessions = "NC_"+str(seed)+"_"+str(p)
        genbankname = fileW+"UnequalCodonSamplingSimulation.genbank."+ str(num_iter)+"."+str(seed)+"."+str(p)+filenote
    elif simulation_type == "mutate":
        name = "synonomous_mutations_Simulation"
        description = "Synonomous mutations are applied to sequences with a predisposition to form LCRs. From file: "\
        + CDSfasta2 + ", Seed = " + str(seed) + ", num_mutations = "+ str(num_mutations)
        organism = "Synonomous mutations"
        accessions = "NC_"+str(seed)+str(num_mutations)
        genbankname = fileW+"SynonomousMutationsSimulation.genbank."+ str(seed) + "." + str(num_mutations)+filenote

    record = SeqRecord(chromosome_object, id = dt_string, name = name,\
            description = description, annotations = {"molecule_type": "DNA", "accessions": [accessions], "organism": organism})
    
    for rec in SeqIO.parse(CDSfasta, "fasta"):
        print("adding features")
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

#this function takes arguments from the command line, and calls the proper functions to generate a genbank file
#if OnlyGenbank == T, existing_fasta must be specified
#if OnlyGenbank == F, existing_fasta can still be specified if simulation_type == mutate
#codonbias can be "True" or "False", if true it will
def CallonFunctions(simulation_type, fileW="", OnlyGenbank="F", existing_fasta="", species="", codonpropfile="", codonclassbias="True", len_protein="100", num_iterations = "100000", seed = "10", p = "27.8", num_mutations = "1000"):

    print("started running")
    #check to make sure arguments are correct
    if simulation_type != "random" and simulation_type != "biased" and simulation_type != "mutate":
        raise ValueError("simulation_type must be 'random', 'biased', or 'mutate'")
    if OnlyGenbank == "":
        OnlyGenbank == "F"
    elif OnlyGenbank != "F" and OnlyGenbank != "T":
        raise ValueError("Only Genbank must be 'T' or 'F'")
    if num_iterations == "":
        num_iterations = 100000
    elif type(eval(num_iterations)) != int:
        raise ValueError("num_iterations must be integer, not " + str(type(eval(num_iterations))))
    elif type(num_iterations) == str:
        num_iterations = int(num_iterations)
    if seed == "":
        seed = 10
    elif type(eval(seed)) != int:
        raise ValueError("seed must be integer, not " + str(type(eval(seed))))
    elif type(seed) == str:
        seed = int(seed) 
    if p == "":
        p = 27.8
    elif type(eval(p)) != int and type(eval(p)) != float:
        raise ValueError("p value must be int or float, not " + str(type(eval(p))))
    elif type(p) == str:
        try:
            p = int(p)
        except ValueError:
            p = float(p)
    if num_mutations == "":
        num_mutations = 1000
    elif type(eval(num_mutations)) != int:
        raise VaueError("num_mutations must be integer, not " + str(type(eval(num_mutations))))
    elif type(num_mutations) == str:
        num_mutations = int(num_mutations)
    if len_protein == "" and species == "":
        len_protein = 100
    elif type(eval(len_protein)) != int:
        raise ValueError("len_protein must be integer, not " + str(type(eval(num_iterations))))
    elif type(len_protein) == str:
        len_protein = int(len_protein)
    if species == "" and codonpropfile != "":
        raise ValueError("need to specify a species associated with the codon proportions file")
    #if species != "" and codonpropfile == "":
    #    raise ValueError("need to specify a codonpropfile associated with the given species")
    #END of argument check

    #now if a species and initial codon props were given, specify p to give a proportion of LCRs as seen in that organism
    # also change the length of the protein to be the average length seen in that organism
    # also going to need to get the codon bias as a list from the nt_cod_aa_props file
    if codonpropfile != "":
        codonbias = parsecodonpropfile(codonpropfile)
    else:
        codonbias = ""

    if species == "hs":
        len_protein = 317
        p = 0.156
        codontypefactor = [0.0572, 1.4248, 0.5158]
    elif species == "sc":
        len_protein = 467
        p = 0.236
        codontypefactor = [0.0647, 1.2661, 0.7798]
    elif species == "dm":
        len_protein = 500
        p = 0.214
        codontypefactor = [0.0124, 1.1904, 0.8413]
    elif species == "ce":
        len_protein = 435
        p = 0.215
        codontypefactor = [0.0494, 1.3341, 0.6101]
    elif species == "at":
        len_protein = 410
        p = 0.197
        codontypefactor = [0.0401, 1.3495, 0.5664]


    if codonclassbias == "False":
        codontypefactor = [1, 1, 1]  #if you do not want to change the proportions based on codontype, then you don't have to.  #only needed for UneqSamp

    if OnlyGenbank == "T":
        MakeGenBankFile(existing_fasta, simulation_type, num_iterations, seed, p, num_mutations, len_protein, species, codonclassbias, codontypefactor, fileW) 
    else:
        if simulation_type == "random":
            CDSfasta = RandomSeqSimulation(num_iterations, seed, fileW)
        elif simulation_type == "biased":
            CDSfasta = UnequalCodonSampling(num_iterations, seed, p, len_protein, species, codonbias, codontypefactor, fileW)
        elif simulation_type == "mutate":
            if existing_fasta == "":
                CDSfasta0 = UnequalCodonSampling(num_iterations, seed, p, len_protein, species, codonbias, codontypefactor, fileW)      #if the corresponding biased fasta does not exist you must make it first
                MakeTranslatedFasta(CDSfasta0)                                 #also make its translated fasta
                MakeGenBankFile(CDSfasta0, "biased", num_iterations, seed, p, num_mutations, len_protein, species, codonclassbias, codontypefactor, fileW)  #make its genbank file too just for comparison
                CDSfasta = ForceSynonomousMutations(CDSfasta0, num_mutations, seed, species, fileW)  #now make the actual mutated fasta
            else:   #if an exisiting fasta file is specified to make mutated seq from                                        
                CDSfasta = ForceSynonomousMutations(existing_fasta, num_mutations, seed, species, fileW)
        #make final translated fasta and genbank files
        MakeTranslatedFasta(CDSfasta)                                             
        MakeGenBankFile(CDSfasta, simulation_type, num_iterations, seed, p, num_mutations, len_protein, species, codonclassbias, codontypefactor, fileW)  #existing fasta argument is just for 'mutate'

#CallonFunctions(input("simulation_type: "), input("OnlyGenbank(d='F'):"), input("existing_fasta(d=''): "),\
#input("num_iterations(d=100000): "), input("seed(d=10): "), input("p(d=27.8): "), input("num_mutations(d=1000): "))
if len(argv) < 4:
    print("arguments are: simulation_type ('mutate' or 'random'), species ('sc'), existing_fasta (fastafile to mutate), fileW (name of file to write to)")
else:
    CallonFunctions(argv[1], species=argv[2], existing_fasta = argv[3], fileW = argv[4])
#CallonFunctions("mutate", species="sc", codonpropfile = "proj1/nt_cod_aa_proportions/genomewidecounts/Saccharomyces_cerevisiae_codon_count")
#CallonFunctions("mutate", species="sc", existing_fasta = "/home/JosLaptop/randomfastafiles/increasingslippage/IncreasingSlippageSimulation.fasta.100000.10.0.358.initbioprops_sc_codonclassbias")
#CallonFunctions("mutate", species="hs", existing_fasta = "/home/JosLaptop/randomfastafiles/increasingslippage/IncreasingSlippageSimulation.fasta.100000.10.0.238.initbioprops_hs_codonclassbias")
#CallonFunctions("mutate", species="dm", existing_fasta = "/home/JosLaptop/randomfastafiles/increasingslippage/IncreasingSlippageSimulation.fasta.100000.10.0.319.initbioprops_dm_codonclassbias")
#CallonFunctions("mutate", species="ce", existing_fasta = "/home/JosLaptop/randomfastafiles/increasingslippage/IncreasingSlippageSimulation.fasta.100000.10.0.31.initbioprops_ce_codonclassbias")
#CallonFunctions("mutate", species="at", existing_fasta = "/home/JosLaptop/randomfastafiles/increasingslippage/IncreasingSlippageSimulation.fasta.100000.10.0.305.initbioprops_at_codonclassbias")
#
#CallonFunctions("mutate", species="sc", existing_fasta = "/home/JosLaptop/randomfastafiles/biased/UnequalCodonSamplingSimulation.fasta.100000.10.0.236_initbioprops_sc_codonclassbias")
#CallonFunctions("mutate", species="hs", existing_fasta = "/home/JosLaptop/randomfastafiles/biased/UnequalCodonSamplingSimulation.fasta.100000.10.0.156_initbioprops_hs_codonclassbias")
#CallonFunctions("mutate", species="dm", existing_fasta = "/home/JosLaptop/randomfastafiles/biased/UnequalCodonSamplingSimulation.fasta.100000.10.0.214_initbioprops_dm_codonclassbias")
#CallonFunctions("mutate", species="ce", existing_fasta = "/home/JosLaptop/randomfastafiles/biased/UnequalCodonSamplingSimulation.fasta.100000.10.0.215_initbioprops_ce_codonclassbias")
#CallonFunctions("mutate", species="at", existing_fasta = "/home/JosLaptop/randomfastafiles/biased/UnequalCodonSamplingSimulation.fasta.100000.10.0.197_initbioprops_at_codonclassbias")
#CallonFunctions("mutate", existing_fasta = "/home/JosLaptop/randomfastafiles/biased/UnequalCodonSamplingSimulation.fasta.100000.10.0.236_initbioprops_sc_codonclassbias")


#CallonFunctions("mutate", species="sc", existing_fasta = "/home/JosLaptop/randomfastafiles/biased/thegoodones/UnequalCodonSamplingSimulation.fasta.100000.10.0.236_initbioprops_sc", fileW = "uneqsampinitbioprops_")
#CallonFunctions("mutate", species="hs", existing_fasta = "/home/JosLaptop/randomfastafiles/biased/thegoodones/UnequalCodonSamplingSimulation.fasta.100000.10.0.156_initbioprops_hs", fileW = "uneqsampinitbioprops_")
#CallonFunctions("mutate", species="ce", existing_fasta = "/home/JosLaptop/randomfastafiles/biased/thegoodones/UnequalCodonSamplingSimulation.fasta.100000.10.0.215_initbioprops_ce", fileW = "uneqsampinitbioprops_")
#CallonFunctions("mutate", species="dm", existing_fasta = "/home/JosLaptop/randomfastafiles/biased/thegoodones/UnequalCodonSamplingSimulation.fasta.100000.10.0.214_initbioprops_dm", fileW = "uneqsampinitbioprops_")
#CallonFunctions("mutate", species="at", existing_fasta = "/home/JosLaptop/randomfastafiles/biased/thegoodones/UnequalCodonSamplingSimulation.fasta.100000.10.0.197_initbioprops_at", fileW = "uneqsampinitbioprops_")
#
#CallonFunctions("mutate", species="sc", existing_fasta = "/home/JosLaptop/randomfastafiles/increasingslippage/thegoodones/IncreasingSlippageSimulation.fasta.100000.10.0.358.initbioprops_sc", fileW = "incslipinitbioprops_")
#CallonFunctions("mutate", species="hs", existing_fasta = "/home/JosLaptop/randomfastafiles/increasingslippage/thegoodones/IncreasingSlippageSimulation.fasta.100000.10.0.238.initbioprops_hs", fileW = "incslipinitbioprops_")
#CallonFunctions("mutate", species="ce", existing_fasta = "/home/JosLaptop/randomfastafiles/increasingslippage/thegoodones/IncreasingSlippageSimulation.fasta.100000.10.0.31.initbioprops_ce", fileW = "incslipinitbioprops_")
#CallonFunctions("mutate", species="dm", existing_fasta = "/home/JosLaptop/randomfastafiles/increasingslippage/thegoodones/IncreasingSlippageSimulation.fasta.100000.10.0.319.initbioprops_dm", fileW = "incslipinitbioprops_")
#CallonFunctions("mutate", species="at", existing_fasta = "/home/JosLaptop/randomfastafiles/increasingslippage/thegoodones/IncreasingSlippageSimulation.fasta.100000.10.0.305.initbioprops_at", fileW = "incslipinitbioprops_")


#CallonFunctions("mutate", species="sc", existing_fasta = "/home/JosLaptop/incslipsim_speciesspecific_improvedcodclassbias.100000.10.sc.fasta" , fileW = "/home/JosLaptop/improvedcodclass_initbioprops_")
#CallonFunctions("mutate", species="hs", existing_fasta = "/home/JosLaptop/incslipsim_speciesspecific_improvedcodclassbias.100000.10.hs.fasta" , fileW = "/home/JosLaptop/improvedcodclass_initbioprops_")
#CallonFunctions("mutate", species="dm", existing_fasta = "/home/JosLaptop/incslipsim_speciesspecific_improvedcodclassbias.100000.10.dm.fasta" , fileW = "/home/JosLaptop/improvedcodclass_initbioprops_")
#CallonFunctions("mutate", species="ce", existing_fasta = "/home/JosLaptop/incslipsim_speciesspecific_improvedcodclassbias.100000.10.ce.fasta" , fileW = "/home/JosLaptop/improvedcodclass_initbioprops_")
#CallonFunctions("mutate", species="at", existing_fasta = "/home/JosLaptop/incslipsim_speciesspecific_improvedcodclassbias.100000.10.at.fasta" , fileW = "/home/JosLaptop/improvedcodclass_initbioprops_")



#CallonFunctions("biased", species=argv[1], codonpropfile=argv[2])
#CallonFunctions("biased", species=argv[1], codonpropfile=argv[2])
