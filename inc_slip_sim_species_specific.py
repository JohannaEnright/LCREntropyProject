#!/usr/bin/python3

import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from datetime import datetime
import re
import copy
from sys import argv

#global variable of list of codons
codons = ["TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG",\
          "TAT", "TAC", "TGT", "TGC", "TGG", "CTT", "CTC", "CTA", "CTG",\
          "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT",\
          "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC",\
          "ACA", "ACG", "AAT", "AAC","AAA", "AAG", "AGT", "AGC", "AGA",\
          "AGG","GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",\
           "GAT","GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"]



#takes a file of codon proportions and returns a list of their proportions (in the order corresponding to "codons" above)
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



#this is the second part to the simulation.
#100 codon long sequences are created with a bias towards creating LCRs 
#(more likely to sample the same codon you are already on)
#p is the probability of resampling the same codon
# len_protein is the length you want you protein to be
def increasingchanceofslippage(num_iterations, codonbias, len_protein, slope, seed, species, fileW):  #will need a samplingweightparameter once I figure out what to make it

    print("generating sequences with increasing chance of slippage")
    random.seed(seed)  #set the seed, default 10

    if codonbias == "":   #build an array of all equal proportions if the initial proportions were not specified
        filenote = "initpropsequal"
        codonbias = []
        for prob in range(61):
            codonbias.append(1)    #build the initial w2 as all equal weights to start (even after the 1st codon sample)  (cuz the array takes relative numbers)
    else:
        filenote= "initbioprops_"+species

    allCDS = []
    j = 0
  #  print(codonbias, end="\n\n")
    while j < num_iterations: #this is how many proteins you want in your fasta and genbank files
        w2 = copy.copy(codonbias) ##reset it back to being the initial proportions each time  #makes a non-aliased copy of codonbias
        codons100 = []  #this is where all the codons will get added
        cod1 = random.choices(codons, weights=w2, k=1)[0]  #choose initial codon with equal weights   
        cod_index = codons.index(cod1)   #this is so I change the correct proportion later
        codons100.append(cod1)     #add initial codon to start of sequence

        i = 1
        dna_sequence = ""
        while i < len_protein:      #sample 99 more codons to make a length of 100
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
        #string codons together
        for c in codons100:
            dna_sequence += c
        allCDS.append(dna_sequence)
       # print(j)
        j += 1

    #write all coding sequences to fasta file and return file name

    to_file_DNA = ""
    ID = 0
    for sequence in allCDS:
        to_file_DNA += ">" + str(ID) + "\n"
        to_file_DNA += sequence + "TAA" + "\n\n"    #each coding sequence needs a stop codon because in the next program this is run through it removes the stop codon
        ID += 1
    tfDNA2 = fileW+"IncreasingSlippageSimulation.fasta"+ "."+str(num_iterations)+"."+str(seed)+"."+str(slope)+"."+filenote
    with open(tfDNA2, "w") as fh:
        fh.write(to_file_DNA)
    return tfDNA2

#this is the third part of the simulation.
#DNA is hammered with synonomous substitutions (only at 3rd codon position)
# mutprob is a list of the probability of each type of mutation, these probabilities are relatively weighted (with the probability of transitions being higher than transversions)
# exp.  mutprob = [transition, transversion],  mutprob = [5, 1]  #relative weights
def ForceSynonomousMutations(Unequal_Sampling_fasta, num_mutations, mutprob, seed):
   
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

    allCDS = []
    to_new_fasta = ""
    counter = 0
    for record in SeqIO.parse(Unequal_Sampling_fasta, "fasta"):   #takes from a previously created fasta of sequences biased towards LCR formation
        counter+=1
        print(counter)
        ID = record.id
        seq = record.seq[:-3]      #-3 so the TAA stop is removed (don't want to mutate)
        #print(len(seq))

        indv_codons = []
        #pick out indvl codons from sequence (str --> list)
        cod = ""
        for i in range(0, len(seq)-2, 3):                   
            cod += seq[i] + seq[i+1] + seq[i+2]  #here we are just breaking the sequence up into indvl codons
            indv_codons.append(cod)             #add codons to list (of 100)
            cod = ""
        #now do while loop to randomly mutate a random codon (1-100) at a nt position 3, with either: A, C, G, or T
        mutations = 0
        while mutations < num_mutations:
            codon_num = random.randint(0, (len(seq)/3)-1)    #len(seq)-1 because we will use this number for indexing, and 0 indexing is a thing
            codon_from_seq = indv_codons[codon_num]    #a random codon from the sequence
            index_codon = random.randrange(0,3,2)   #number will be only 0 or two (to mutate only first or third nt, as the second nt will never return a synonomous codon with 1 mutation)
            nuctomutate = codon_from_seq[index_codon]  #will index randomly at either the 0th or 2nd position to choose which nt to mutate (mutating the second will never work as serine needs 2 mutations to change codons anyway)
            #now decide what the probability array will be based on what 3rdnuc is
            prob_transition = mutprob[0] 
            prob_transversion = mutprob[1]      
            if nuctomutate == "A":      # nuc = ["A", "C", "G", "T"]  just to remember the order
                nucprobs = [0, prob_transversion, prob_transition, prob_transversion]    ##made A 0, just because it will continue iterating until it successfully mutates anyway
            elif nuctomutate == "C":
                nucprobs = [prob_transversion, 0, prob_transversion, prob_transition]
            elif nuctomutate == "G":
                nucprobs = [prob_transition, prob_transversion, 0, prob_transversion]
            elif nuctomutate == "T":
                nucprobs = [prob_transversion, prob_transition, prob_transversion, 0]

            nt = random.choices(nuc, weights=nucprobs, k=1)[0]   #nucprobs is an array of the probabilities of choosing each nt based on transition and transversion, order corresponds to 'nuc'
            #now peice the codon back together
            if index_codon == 0:
                mutant_codon = nt + codon_from_seq[1:]
            elif index_codon == 2:
                mutant_codon = codon_from_seq[0:2] + nt

            #if a synonomous mutation was made, replace the original codon with the mutant codon 
            if (mutant_codon in synon_codons[codon_from_seq]): 
                indv_codons[codon_num] = mutant_codon
                mutations += 1  #this should be under if statement so that it is only increased when the change is actually made
        
        #now put mutated sequence back into string format and write new fasta
        mutated_seq = ""
        for c in indv_codons:
            mutated_seq += c
        mutated_seq += "TAA"
        #print(mutated_seq)

        to_new_fasta += ">" + ID + "\n" + mutated_seq + "\n\n"

    tfDNA3 = "ti"+str(mutprob[0])+"tv"+str(mutprob[1])+"_MutatedLCRSimulation.fasta"+"."+str(num_mutations)+"."+str(seed)
    with open(tfDNA3, "w") as fh:
        fh.write(to_new_fasta)

    return tfDNA3

##ForceSynonomousMutations("randomfastafiles/biased/UnequalCodonSamplingSimulation.fasta.100000.10.27.8", 1000, 10, 1)


#this just makes another fasta file with corresponding protein sequences
#handy if you just want to run it on SEG and not the entire entropy program
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

#now generate genbank file from these random DNA sequences
#simulation type indicates whether the sequences are coming from truly random sequences or biased sequences
#CDSfasta is the fasta file that is converted to a genbank file
#CDSfasta2 is only for mutated type (to show which file the sequences that were mutated came from)
#simulation type is random, biased, or mutate
#num_iter is the number of sequences in the fasta file (argument is used to name file)
#seed (argument is used to name gb file)
#p (argument is used to name gb file)
def MakeGenBankFile(CDSfasta, simulation_type, num_iter, seed, slope, num_mutations, mutprob, species, fileW, CDSfasta2=""):

    print("making genbank file")
    #concatenate all DNA seqs and make it the chromosome sequence
    chromosome_string = ""
    for rec in SeqIO.parse(CDSfasta, "fasta"):
        chromosome_string += rec.seq
    chromosome_object = Seq(chromosome_string)
    
    if species != "":
        filenote = "initbioprops_"+species
    else:
        filenote = "initequalprops"
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
        name = "Increasing_chance_of_slippage_Simulation"
        description = "Increasing slippage change as the repeat length increase. Num_iter = "+ str(num_iter)+", Seed = " + str(seed) + ", slope = " + str(slope) + "  "+filenote
        organism = "Increasing Slippage"+" "+filenote
        accessions = "NC_"+str(seed)+"_"+str(slope)
        genbankname = fileW+"IncreasingSlippageSimulation.genbank."+ str(num_iter)+"."+str(seed)+"."+str(slope)+"_"+filenote
    elif simulation_type == "mutate":
        name = "synonomous_mutations_Simulation"
        description = "Synonomous mutations are applied to sequences with a predisposition to form LCRs from file "\
                + CDSfasta2 + ", Seed = " + str(seed) + ", num_mutations = "+ str(num_mutations) + " transition/transversion = "+ str(mutprob[0])+"/"+str(mutprob[1])+"  "+filenote
        organism = "Synonomous mutations"
        accessions = "NC_"+str(seed)+str(num_mutations)
        genbankname = fileW+"Ti"+str(mutprob[0])+"Tv"+str(mutprob[1])+"_SynonomousMutationsSimulation.genbank."+ str(seed) + "." + str(num_mutations)

    record = SeqRecord(chromosome_object, id = dt_string, name = name,\
            description = description, annotations = {"molecule_type": "DNA", "accessions": [accessions], "organism": organism})
    
    for rec in SeqIO.parse(CDSfasta, "fasta"):
        print("adding features")
        ID = rec.id
        DNAseq = rec.seq
        cds_indexfactor = len(DNAseq)  #so that it gives the correct start and end positions as the cds  (len(proseq)*3)+3
        start = int(ID)*cds_indexfactor   
        end = start + cds_indexfactor
        translation = Seq(DNAseq).translate()[:-1]
       
        #create feature annotations
        feature = SeqFeature(FeatureLocation(start=start, end=end), type = "CDS", strand=1,\
                  qualifiers= {"protein_id":[ID], "translation":[translation], "db_xref":["GeneID:"+ID]})
        record.features.append(feature)
    
    #save as genbank file
    with open(genbankname, "w") as fh:
        SeqIO.write(record, fh, "genbank")
    print("genbank file is written")


##this function takes arguments from the command line, and calls the proper functions to generate a genbank file
##if OnlyGenbank == T, existing_fasta must be specified
##if OnlyGenbank == F, existing_fasta can still be specified if simulation_type == mutate
#len_protein is the length you want each protein sequence to be
#slope will be how quickly you want the LCR to slip (amount by which resampling a codon increases by if there is a repeat)
# it is set to 54.5 because this generates the same proportion of LCRs as would be seen in a yeast proteome
# mutprob is the relative weighting of transition to transversions 
# codonpropfile is going to be a file from nt_cod_aa_bias output which contains all of the codon proportions
# input the prop file which corresponds to a particular species proportions
def CallonFunctions(simulation_type, codonpropfile="", species="", existingfasta="", fileW="", num_iterations = 100000, len_protein=100, slope=54.5, mutprob=[4,1], num_mutations = 1000, seed=10):

    if codonpropfile != "":
        codonbias = parsecodonpropfile(codonpropfile)   #this will give a list of weights in the same order as the 'codons' global list
    else:
        codonbias = ""

    #species can be 'hs' (humans), 'sc' (yeast), 'dm' (drosophila), 'at' (arabidopsis), 'ce' (c.elegans)
    if species == "hs":
       # len_protein = 317
        len_protein = 584
       # slope = 0.238       #trial and error (must test to figure out)
        slope = 0.192       #trial and error (must test to figure out)
    elif species == "sc":
       # len_protein = 467
        len_protein = 488
        #slope = 0.358
        slope = 0.352
    elif species == "dm":
      #  len_protein = 500
        len_protein = 533
        #slope = 0.319
        slope = 0.311
    elif species == "ce":
      #  len_protein = 435
        len_protein = 409
        #slope = 0.31
        slope = 0.316
    elif species == "at":
      #  len_protein = 410
        len_protein = 406
       # slope = 0.305
        slope = 0.307

    if simulation_type == "biased":
        if existingfasta != "":
            lcrbiasedfasta = existingfasta
        else:
            lcrbiasedfasta = increasingchanceofslippage(num_iterations, codonbias, len_protein, slope, seed, species, fileW)
            MakeTranslatedFasta(lcrbiasedfasta)   #translate these sequences into protein sequences
#    elif simulation_type == "mutate":
#        lcrbiasedfasta = ForceSynonomousMutations("/home/JosLaptop/randomfastafiles/increasingslippage/IncreasingSlippageSimulation.fasta.100000.10.54.5.1", num_mutations, mutprob, seed)

    MakeGenBankFile(lcrbiasedfasta, simulation_type, num_iterations, seed, slope, num_mutations, mutprob, species, fileW)       #make the Genbank file!


#CallonFunctions("biased", species="sc", codonpropfile="")
#CallonFunctions("biased", fileW = "/home/JosLaptop/incslip_newproteinlengthfiles/", species="sc", codonpropfile="/home/JosLaptop/proj1/nt_cod_aa_proportions/genomewidecounts/Saccharomyces_cerevisiae_codon_count")
#CallonFunctions("biased", fileW = "/home/JosLaptop/incslip_newproteinlengthfiles/", species="hs", codonpropfile="/home/JosLaptop/proj1/nt_cod_aa_proportions/genomewidecounts/Homo_sapiens_codon_count")
#CallonFunctions("biased", fileW = "/home/JosLaptop/incslip_newproteinlengthfiles/", species="dm", codonpropfile="/home/JosLaptop/proj1/nt_cod_aa_proportions/genomewidecounts/Drosophila_melanogaster_codon_count")
#CallonFunctions("biased", fileW = "/home/JosLaptop/incslip_newproteinlengthfiles/", species="ce", codonpropfile="/home/JosLaptop/proj1/nt_cod_aa_proportions/genomewidecounts/Caenorhabditis_elegans_codon_count")
#CallonFunctions("biased", fileW = "/home/JosLaptop/incslip_newproteinlengthfiles/", species="at", codonpropfile="/home/JosLaptop/proj1/nt_cod_aa_proportions/genomewidecounts/Arabidopsis_thaliana_codon_count")
if len(argv) < 4:
    print("arguments are: 'biased', species='at', codonpropfile=''")
else:
    CallonFunctions(argv[1], species=argv[2], codonpropfile= argv[3])
#simulation_type:  'biased' or 'mutate'

