#!/usr/bin/python3

from sys import argv
from Bio import SeqIO
import re
from Bio.Seq import Seq

#this takes a genbank file and returns the name of the organism
#"genus species"
#this is literally so I can see it at the top of the file and be
#less confused
def OrganismName(fileR):

    with open(fileR, "r") as fh:
        for line in fh:
            if re.search(r"^DEFINITION", line):
                c = 0
                name = ""
                for let in line:
                    if let == " ":
                         c += 1
                    if c == 4:
                        break
                    if c > 1:
                        if let == " ":
                            name += "_"
                        else:
                            name += let
                break
    return name[1:]


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


#takes genbank file to read from 
def countcodtypesCDS(fileR, fileW=""):

    #print("started running")
    totalcodoncount = {"c1":0, "c2":0, "c3":0}
    totalcods = 0
    AllGeneIDs = {}  #no duplicate GeneIDs!

    ##Get the GeneSeqeunce to count the codons from
    for record in SeqIO.parse(fileR, "gb"):
        try:
            Acc = record.annotations["accessions"][0]
            full_chrom_seq = str(record.seq)       
          #  print("Acc ", Acc)
        except: #some sort of Bio.SeqError
          #  print("getting other seqs for humans and drosophila")
            Organism = OrganismName(fileR)
            if Organism == "Homo_sapiens":
                filepath = "/scratch/Human/Unplaced.fna"
            elif Organism == "Drosophila_melanogaster":
                filepath = "/scratch/d.melanogaster/drosophilaUnplaced.fna"
            else:
                raise FileNotFoundError("need a fasta file path for "+ Organism)

                for newrecord in SeqIO.parse(filepath, "fasta"):
                    if newrecord.id == ID:
                        full_chrom_seq = str(newrecord.seq)

        for i in range(len(record.features)):
            P = record.features[i]
            if P.type == "CDS":
                strand = P.strand                            #will be +1 or -1
                GeneStart = P.location.nofuzzy_start + 1     #1 indexed
                GeneEnd = P.location.nofuzzy_end  
                intron = False
                exonPos = {}

                try:
                    for j in P.qualifiers["db_xref"]:
                        find = re.search(r"GeneID:\d+", j)
                        if find != None:
                            GeneID = find.group()
                            break
                except KeyError:
                    raise KeyError("It is not finding the GeneID")
                try:
                    ProID = P.qualifiers["protein_id"][0]
                    aaSeq = P.qualifiers["translation"][0]
                except KeyError:
                    continue

                if GeneID in AllGeneIDs: #if same Gene, then check aaSeq length.
                    a = AllGeneIDs[GeneID][-1]  #this indexes the aaSeq
                    b = aaSeq
                    if len(a) >= len(b): #if previous aaSeq(withLCR)is longer, keep it and continue, otherwise check if it has LCR
                        continue     

                if P.location_operator == "join": #otherwise add in new information for that geneID
                    intron = True                 #or get info for that unique gene

                if intron == True:
                    for exon in P.location.parts:
                        ex = re.search(r"(\d+):[><]*(\d+)", str(exon)) #< or < cuz some exons are partial sequences
                        exstrand = re.search(r"\([+-]\)", str(exon))
                        exonPos[(int(ex.group(1))+1, int(ex.group(2)))] = exstrand.group() #1 indexed

                #get geneSeq
                if intron == True:    #concatenate the exons together into GeneSeq variable
                    GeneSeq = ""
                    for tup, s in exonPos.items():
                        forward = full_chrom_seq[tup[0]-1:tup[1]]
                        if s == "(+)":
                            GeneSeq += forward
                        elif s == "(-)":
                            GeneSeq += str(Seq(forward).reverse_complement())
                    #if there are no introns 
                else:
                    GeneSeq = full_chrom_seq[GeneStart-1: GeneEnd]
                    if strand == -1:
                        GeneSeq = str(Seq(GeneSeq).reverse_complement())   
                GeneSeq = GeneSeq[:-3].lower()  #remove stop codon and ensure all in lower case

                ##Don't consider sequences with a stop codon or a which are not 3X length of the protein sequence
                #Ensure DNAseq is 3X length of aaSeq
                if len(aaSeq)*3 != len(GeneSeq):
                    continue  #these will not be checked for LCRs
                #check for stop codons in ORF, these sequences will not be included
                codon_list = Split_Into_Codons(GeneSeq)
                if "tag" in codon_list or "tga" in codon_list or "taa" in codon_list:
                    continue  #has a stop codon in the middle
                ##you now have the same GeneSeqs which were used in CombinedUsingSysCall_Codon3

                ##count the codon types
                dictofcodonclasscounts = countcodonclasscategories(GeneSeq)  #{"c1":7, "c2":3, "c3":9}
                #print("counting codon types") 

                #get info from functions
                c1codnum = dictofcodonclasscounts["c1"]
                c2codnum = dictofcodonclasscounts["c2"]
                c3codnum = dictofcodonclasscounts["c3"]
                totalcodoncount["c1"] += c1codnum
                totalcodoncount["c2"] += c2codnum
                totalcodoncount["c3"] += c3codnum
                totalcods += c1codnum + c2codnum + c3codnum

   # print(totalcodoncount)
   # print(totalcods)
    ##Write info to files
    #countfilename = fileW+"codclasscount_allCDS"

    #gather the info into a string
    filenote = "#File contains number and proportion of codon class in all CDS (same as CombindUsingSysCall_Codon3)\n#Info came from file: "+fileR+"\n"
    header = "totcod\tc1\tc2\tc3\tpropc1\tpropc2\tpropc3\n"
   # print(filenote, header)

    to_totcountfile = filenote + header + str(totalcods) + "\t"
    #print(to_totcountfile)
    for codnum in totalcodoncount.values():
        to_totcountfile += str(codnum)+"\t"
    for codnum in totalcodoncount.values():
        to_totcountfile += str(round(codnum/totalcods,4)) + "\t"
    to_totcountfile += "\n\n"

    print(to_totcountfile)
    #with open(countfilename, "w") as fh:
    #    fh.write(to_totcountfile)
    #print("files are written")



#countcodtypesCDS("/scratch/S.cer_gen_all")
if len(argv) < 2:
    print("takes file containing all the genbank files from an organism")
else:
    countcodtypesCDS(argv[1])
