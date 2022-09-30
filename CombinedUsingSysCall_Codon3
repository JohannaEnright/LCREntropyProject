#!/usr/bin/python3

import sys
import warnings
from Bio import SeqIO
import subprocess
from Bio.Seq import Seq
import re
from math import log
import tempfile

print(sys.argv[1:])
print("welcome")



#will take gb file and find LCRs within every gene (by converting to fasta and running
#through seg)
#then will find corresponding LCR positions within the protein
#will calulate the entropy of this subsegment of protein
#map it to a dictionary geneID:(information)
#write this dictionary to a fil


#takes the start position of DNAseq (a), the end position of DNAseq (b)
#and the length of DNA seq (c)
#returns the new a, b,  and c when trimmed to line up with the start and end of a codon
def trim_dna_ends(a, b, c):
    if (a-1)%3 == 0 and b%3 == 0: #both are fine (first and third)
        pass
    elif (a-1)%3 == 0 and (b+1)%3 == 0: #start in 1st pos, end in 2nd pos
        b = b-2
        c = c-2
    elif (a-1)%3 == 0 and (b+2)%3 == 0: #start in 1st pos, end is in 1st pos
        b = b-1
        c = c-1
    elif (a-2)%3 == 0 and b%3 == 0: #start in 2nd pos, end in 3rd pos
        a = a+2
        c = c-2
    elif (a-2)%3 == 0 and (b+1)%3 == 0: #start in 2nd pos, end in 2nd pos
        a = a+2
        b = b-2
        c = c-4
    elif (a-2)%3 == 0 and (b-1)%3 == 0: #start in 2nd pos, end in 1st pos
        a = a+2
        b = b-1
        c = c-3
    elif a%3 == 0 and b%3 == 0: #start in 3rd pos, end in 3rd pos
        a = a+1
        c = c-1
    elif a%3 == 0 and (b+1)%3 == 0: #start in 3rd pos, end in 2nd pos
        a = a+1
        b = b-2
        c = c-3
    else:  #start in 3rd pos, end in 1st pos
        a = a+1
        b = b-1
        c = c-2
    return (a,b,c)


#takes a Fasta file (output from seg) and returns the positions of the longest LCR in the file
#as a list of 1 tuple
#if 2 LCRs are the same length, it returns a list of >1 tuples
def LongestLCR(FastaFile, KeepEnds):
    
    LCRs = {}
    for record in SeqIO.parse(FastaFile, "fasta"):
        Index = record.id.index("(")
        Pos = record.id[Index+1:-1]
        Pos = Pos.split("-")
        a = int(Pos[0])         #start pos of LCR
        b = int(Pos[1])         #end pos of LCR
        c = len(str(record.seq))  #length of LCR

        if KeepEnds == False:  #make sure LCR starts at 1st nt at first codon and ends at 3rd nt of last codon
            a = trim_dna_ends(a,b,c)[0]
            b = trim_dna_ends(a,b,c)[1]
            c = trim_dna_ends(a,b,c)[2]
        LCRs[(a,b)] = c  
            
    result = []
    Max = 0
    for num in LCRs.values():
        if num > Max:
            Max = num
    if Max == 0:
        return "noLCR"
    else:
        for key, value in LCRs.items():
            if value == Max:
                result.append(key)
        return result   #returns a list of 1 tuple or more tuples (only >1 if LCRs are equal length)
                        #if 2 LCRs are the same length you must later compare their entropy
                        #positions could be of geneLCR, codonLCR or proteinLCR


#this is to take a DNA string and return it as a list of codons within the string
def Split_Into_Codons(DNAseq):
    codon_list = []
    for i in range(0, len(DNAseq)-2, 3):
        codon = DNAseq[i] + DNAseq[i+1] + DNAseq[i+2]
        codon_list.append(codon)
    return codon_list

#takes list of codons and converts to string of ascii characters to run through seg
def convert_to_ascii(ProID, Acc, codon_list):

    cod_ascii = {"ttt":'"', "ttc":"#", "tta":"$", "ttg":"%", "tct":"&", "tcc":"'", "tca":"(", "tcg":")",\
                "tat":"*", "tac":"+", "tgt":",", "tgc":"-", "tgg":".", "ctt":"/", "ctc":"0", "cta":"1", "ctg":"2",\
                "cct":"3", "ccc":"4", "cca":"5", "ccg":"6", "cat":"7", "cac":"8", "caa":"9", "cag":":", "cgt":";",\
                "cgc":"<", "cga":"=", "cgg":"?", "att":"@", "atc":"A", "ata":"B", "atg":"C", "act":"D", "acc":"E",\
                "aca":"F", "acg":"G", "aat":"H", "aac":"I","aaa":"J", "aag":"K", "agt":"L", "agc":"M", "aga":"N",\
                "agg":"O","gtt":"P", "gtc":"Q", "gta":"R", "gtg":"S", "gct":"T", "gcc":"U", "gca":"V", "gcg":"W",\
                "gat":"X","gac":"Y", "gaa":"Z", "gag":"[", "ggt":"\\", "ggc":"]", "gga":"^", "ggg":"_"}   
    CodonSeq = ""
    for cod in codon_list:
        try:
            CodonSeq += cod_ascii[cod]
        except KeyError:
            warnings.warn("ambiguous amino acids, "+cod+", in codon.")
            CodonSeq = "NA"  #means there is a stop codon in the ORF or ambig
                             #nts in ORF, cannot properly be converted to ascii
    return CodonSeq

#searches Fastafile (with all LCRs for a protein) and returns the Entropy 
# for the longest LCR (based on matching positions from LongestLCR function)
def Entropy(ProID, Acc, seqString, Type):

    seqString = seqString.lower() #only for dna and protein

    count = {}
    ambig_nts = {"a":["a"], "c":["c"], "g":["g"], "t":["t"], "n":["a","c","g","t"], "m":["a","c"],\
                "r":["a","g"], "w":["a","t"], "s":["c","g"], "y":["c","t"], "k":["g","t"],\
                "v":["a","c","g"], "h":["a","c","t"], "d":["a","g","t"], "b":["c","g","t"]}
    ambig_aas = {"x":["a","r","n","d","c","q","e","g","h","i","l","k","m","f","p","s","t","w","y","v"],\
                 "b":["d","n"], "z":["e","q"], "j":["i","l"], "u":["c"], "o":["k"]}

    if Type == "dna":
        num_let = len(seqString)  #this is for Shannon's E calc later
        count["a"] = 0; count["c"] = 0; count["g"] = 0; count["t"] = 0

        for let in seqString:
            if let in count:
                count[let] += 1
            #ambiguous nts
            else:
                warnings.warn(ProID + " " + Acc + "ambiguous nucleotide: " + let + ", in DNA sequence")
                list_pos_nts = ambig_nts[let]   #if let is not a,c,g,t or ambig, then it will key error
                for nt in list_pos_nts:
                    count[nt] += 1/len(list_pos_nts)

    elif Type == "protein":
        num_let = len(seqString)  #for Shannon's E calc
        count["a"] = 0; count["c"] = 0; count["d"] = 0; count["e"] = 0
        count["f"] = 0; count["g"] = 0; count["h"] = 0; count["i"] = 0
        count["k"] = 0; count["l"] = 0; count["m"] = 0; count["n"] = 0
        count["p"] = 0; count["q"] = 0; count["r"] = 0; count["s"] = 0
        count["t"] = 0; count["v"] = 0; count["w"] = 0; count["y"] = 0

        for let in seqString:
            num_let = len(seqString) 
            if let in count:
                count[let] += 1
            else:
                warnings.warn(ProID + " " + Acc + "ambiguous AA: " + let + ", in protein sequence.")
                list_pos_aas = ambig_aas[let]
                for aa in list_pos_aas:
                    count[aa] += 1/len(list_pos_aas) 
    elif Type == "codon":
        num_let = len(seqString)/3
        #all of the codons (not including stop codons)
        count["ttt"] = 0; count["ttc"] = 0; count["tta"] = 0; count["ttg"] = 0
        count["tct"] = 0; count["tcc"] = 0; count["tca"] = 0; count["tcg"] = 0
        count["tat"] = 0; count["tac"] = 0; count["tgt"] = 0; count["tgc"] = 0
        count["tgg"] = 0; count["ctt"] = 0; count["ctc"] = 0; count["cta"] = 0
        count["ctg"] = 0; count["cct"] = 0; count["ccc"] = 0; count["cca"] = 0
        count["ccg"] = 0; count["cat"] = 0; count["cac"] = 0; count["caa"] = 0
        count["cag"] = 0; count["cgt"] = 0; count["cgc"] = 0; count["cga"] = 0
        count["cgg"] = 0; count["att"] = 0; count["atc"] = 0; count["ata"] = 0
        count["atg"] = 0; count["act"] = 0; count["acc"] = 0; count["aca"] = 0
        count["acg"] = 0; count["aat"] = 0; count["aac"] = 0; count["aaa"] = 0
        count["aag"] = 0; count["agt"] = 0; count["agc"] = 0; count["aga"] = 0
        count["agg"] = 0; count["gtt"] = 0; count["gtc"] = 0; count["gta"] = 0
        count["gtg"] = 0; count["gct"] = 0; count["gcc"] = 0; count["gca"] = 0
        count["gcg"] = 0; count["gat"] = 0; count["gac"] = 0; count["gaa"] = 0
        count["gag"] = 0; count["ggt"] = 0; count["ggc"] = 0; count["gga"] = 0
        count["ggg"] = 0

        codon_list = Split_Into_Codons(seqString)  #split the DNA sequence into respective codons
        #now make a string of ascii characters representing the codons in the original string
        for cod in codon_list:
            if cod in count:  #sequences with stop codons were removed earlier
                count[cod] += 1
            #if it is an ambiguous nucleotide need to make proper codons
            else:
                print(codon_list)
                warnings.warn(ProID + " " + Acc + " ambiguous nt(s): " + cod + ", in codons.")
                possible_cods = []
                possib_nt_1 = ambig_nts[cod[0]]   #list of all possible nts for 1st pos in codon
                possib_nt_2 = ambig_nts[cod[1]]   #list of all possible nts for 2nd pos in codon
                possib_nt_3 = ambig_nts[cod[2]]   #list of all possible nts for 3rd pos in codon
                for char1 in possib_nt_1:
                    for char2 in possib_nt_2:
                        for char3 in possib_nt_3:
                            possible_cods.append(char1 + char2 + char3)
                if "taa" in possible_cods:
                    possible_cods.remove("taa")
                if "tag" in possible_cods:
                    possible_cods.remove("tag")
                if "tga" in possible_cods:
                    possible_cods.remove("tga")
                #now add the fraction for proportion of each codon observed
                for new_cod in possible_cods:
                    count[new_cod] += 1/len(possible_cods)

    E = 0
    for freq in count.values():
        if freq != 0:
            pi = freq/num_let    #len seqString /3 
            E += (pi*(log(pi, 2)))
    E = -1*E
    return round(E, 3)


#spits out corresponding DNA positions of the LCR within the chromosome
#1 indexed
#if 3*aaLength != DNA sequence, then it returns NA because proper mapping cannot take place.
def MapToChrom(direction, strand, longestlcr_pos, GeneStart, GeneEnd, intron, exonPos):

    #all aaseq are 3X the length of the DNA seq
    ans = {}
    if direction == "pro_dna" or direction == "pro_cod"\
        or direction == "cod_pro" or direction == "cod_dna":
        LS = (longestlcr_pos[0][0]*3)-2 #these are AA positions or ascii positions converted to gene positions
        LE = (longestlcr_pos[0][1]*3)
    elif direction == "dna_pro" or direction == "dna_cod":
        LS = longestlcr_pos[0][0]  #these are gene positions
        LE = longestlcr_pos[0][1]

    #no introns
    if intron == False:
        if strand == 1:
            start = GeneStart + LS -1 #lower chrom pos
            end = GeneStart + LE - 1  #higher chrom pos
            st = "(+)"
        elif strand == -1:
            start = GeneEnd - LE + 1 #lower chrom pos
            end = GeneEnd - LS + 1  #higher chrom pos
            st = "(-)"
        ans[(start,end)] = st  #ans is a dictionary (4056, 6078):"(+)"

    #introns
    else: #works for 1, -1, or None strand 
        ChromGene = {}    #will have (chromS,chromE):(geneS,geneE)              
        Next = 1
        for chrom, st in exonPos.items():
            key = (chrom, st)              #looks like: ((234, 55667), "(+)")
            diff = chrom[1] - chrom[0]
            ChromGene[key] = (Next, Next+diff)  #these are the corresponding gene Pos   key : (1,20)
            Next += diff+1

        indexS = 0
        for key, val in ChromGene.items():
            #LCR is within same exon
            if LS >= val[0] and LS < val[1] and LE > val[0] and LE <= val[1]:
                st = key[1]
                if st == "(+)":
                    start = LS-val[0] + key[0][0]
                    end = LE-val[0] + key[0][0]
                elif st == "(-)":
                    end = key[0][1] - (LS-val[0])   #end as in highest chromPos
                    start = key[0][1] - (LE-val[0])
                ans[(start, end)] = st  #a tuple of 2 ints (lowest to highest):"(-)"

            #LCR spans 1 or more exons
            elif LS >= val[0] and LS <= val[1] and LE > val[1]:
                begin = indexS #the index (0 ind) in ChromGene that firstEx is at
                if st == "(+)":
                    firstEx = (key[0][0] + LS-val[0], key[0][1]) #add position onto start, end
                elif st == "(-)":
                    firstEx = (key[0][0], key[0][1] - (LS-val[0])) #start, substract position from end
                ans[firstEx] = st

                indexE = 0
                for k, v in ChromGene.items():  #for loop to find which exon the end of the LCR lies in 
                    #if 1 junction is spanned
                    if LE >= v[0] and LE <= v[1]:
                        finish = indexE  #the index (0 ind) in ChromGene that lastEx is at
                        s = k[1]
                        if s == "(+)":
                            lastEx = (k[0][0], k[0][0] + LE-v[0]) #start, add position onto end
                        elif s == "(-)":
                            lastEx = (k[0][1] - (LE-v[0]), k[0][1]) #subtr pos from end, end
                    indexE += 1                                     #tuple of ints from lowest to highest

                #if more than 1 junction is spanned
                if ((firstEx[1]+1 - firstEx[0]) + (lastEx[1]+1 - lastEx[0])) != (LE+1 - LS):
                    addEx = 0
                    for key2 in ChromGene.keys(): #add the exons in between
                        if (addEx > begin) and (addEx < finish):
                            ans[key2[0]] = key2[1]
                        addEx += 1
                ans[lastEx] = s  #add the end exon
            indexS += 1  #in what exon (indexed by nums)

    return ans  
    #ans will be a dictionary with (int,int):"(+)"   (or "(-)")
    #if could be a length of 1 if there are no introns, or if LCR is within a single exon
    #or it could be a length >1 if LCR spans 1 or more exon junctions


#takes LCR_pos within the gene (tuple) 
#returns LCR positions within the protein (tuple) (1 indexed)
def MapToProtein(LCR_pos):

    AAstart = ((LCR_pos[0][0] - 1) // 3) + 1
    AAend = ((LCR_pos[0][1] -1) // 3) + 1
    return (AAstart, AAend)


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


#this is the main function which calls on the other smaller functions
#takes gb file to read, name of file you want to write to, window length, seg highcut and lowcut parameters, 
#and an optional parameter in case you want to look at all LCRs in a protein instead of choosing the longest one
#stores info on LCR positioning, entropy, corresponding chromsome positioning, entropy for that sequence, parameters used
#writes all info for a protein to a line in a file.

def CombinedSysCall(fileR, direction, W, K1, K2, trim = "10000", AllLCR = False, KeepEnds = False):

    if direction != "pro_dna" and direction != "dna_pro" and direction != "pro_cod"\
        and direction != "cod_pro" and direction != "cod_dna" and direction != "dna_cod" :
        raise NameError("must specify direction: 'dna_pro' or 'pro_dna' or 'pro_cod'\
        or 'cod_pro' or dna_cod or cod_dna")
    if (type(eval(W)) != int and type(eval(W)) != float) or\
        (type(eval(K1)) != int and type(eval(K1)) != float) or\
        (type(eval(K2)) != int and type(eval(K2)) != float):
        raise TypeError("paramters for seg must be int or float")
    if trim == "":
        trim = "10000"
    elif type(eval(trim)) != int and type(eval(trim)) != float:
        raise TypeError("trim must be int or float")
    if AllLCR == "" or AllLCR == "False" or AllLCR == "F":
        AllLCR = False
    elif AllLCR == "True" or AllLCR == "T":
        AllLCR = True
    elif AllLCR != False:
        raise TypeError("AllLCR must be boolean: True, T, False, F")
    if (KeepEnds == "" and direction == "dna_pro") or\
        (KeepEnds == "" and direction == "dna_cod") or\
        KeepEnds == "False" or KeepEnds == "F":  #only applies if going from DNA to other
        KeepEnds = False
    elif KeepEnds == "True" or KeepEnds == "T" or direction == "pro_dna"\
        or direction == "pro_cod" or direction == "cod_pro" or direction == "cod_dna":
        KeepEnds = True  #dont want to drop ends (only for DNA) 
    elif KeepEnds != False: #if they put it as anything else
        raise TypeError("KeepEnds must be boolean: True, T, False, F")

    print("Started running")
    Organism = OrganismName(fileR)
    print(Organism)
    AllGeneIDs = {}  #one copy of each gene ID will be stored here
    ExtraGeneIDs = {}  #extra gene ID copies will be stored here (this is so I can count them up later for checking)
    stopcod_OR_ambignts = []
    unequal_lengths = []
    NoLCR = {}

    #parse genbank file for indv proten sequences and convert to fasta
    dupnotreplaced = 0  
    dup = 0
    total = 0
    replacedORAdded = 0
    hadstopcod = 0
    nolcr_count = 0
    for record in SeqIO.parse(fileR, "gb"):
        try:
            Acc = record.annotations["accessions"][0]
            print("On new gb file: ", Acc)
            full_chrom_seq = str(record.seq)       
        except: #some sort of Bio.SeqError
            print("cannot get full chrom seq")
            #now we need to pull from the fasta files for the unplaced genes
            ID = record.id
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
                total += 1
                strand = P.strand                            #will be +1 or -1
                GeneStart = P.location.nofuzzy_start + 1     #1 indexed
                GeneEnd = P.location.nofuzzy_end  
                intron = False
                exonPos = {}

                #get protein description information
                try:
                    Description = P.qualifiers["note"][0]
                except KeyError:
                    Description = "no description"
                try:
                    Product = P.qualifiers["product"][0]
                except KeyError:
                    Product = "no product info"

                #now get CDS information
                try:
                    for j in P.qualifiers["db_xref"]:
                        find = re.search(r"GeneID:\d+", j)
                        if find != None:
                            GeneID = find.group()
                            break
                except KeyError:
                    GeneID = "LocusTag:" + P.qualifiers["locus_tag"][0]
                    warnings.warn(Acc + "  i: " + str(i) + "  no gene ID")
                try:
                    ProID = P.qualifiers["protein_id"][0]
                    aaSeq = P.qualifiers["translation"][0]
                except KeyError:
                    continue  #if there is no proteinID, this will not be checked for LCRs

                if GeneID in AllGeneIDs: #if same Gene, then check aaSeq length.
                    dup += 1
                    a = AllGeneIDs[GeneID][-1]  #this indexes the aaSeq
                    b = aaSeq
                    if len(a) >= len(b): #if previous aaSeq(withLCR)is longer, keep it and continue, otherwise check if it has LCR
                        ExtraGeneIDs[dupnotreplaced] = (GeneID, ProID, Acc) #if new seq has no LCR, it will not be replaced
                        dupnotreplaced += 1
                        continue     

                if P.location_operator == "join": #otherwise add in new information for that geneID
                    intron = True                 #or get info for that unique gene

                if intron == True:
                    for exon in P.location.parts:
                        ex = re.search(r"(\d+):[><]*(\d+)", str(exon)) #< or < cuz some exons are partial sequences
                        exstrand = re.search(r"\([+-]\)", str(exon))
                        exonPos[(int(ex.group(1))+1, int(ex.group(2)))] = exstrand.group() #1 indexed

                if (strand != 1 and strand != -1) and intron == False:    #exons on different strands
                    raise ValueError("There is no strand type and intron == False")
                if P.location_operator != "join" and P.location_operator != None:
                    raise ValueError("intron is neither 'join nor None")

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

                #Ensure DNAseq is 3X length of aaSeq
                if len(aaSeq)*3 != len(GeneSeq):
                    unequal_lengths.append((ProID, Acc, GeneSeq, aaSeq, str(len(GeneSeq)), str(len(aaSeq))))
                    warnings.warn("DNA seq and aaSeq do not have proportional lengths")
                    continue  #these will not be checked for LCRs

                #check for stop codons in ORF, these sequences will not be included
                codon_list = Split_Into_Codons(GeneSeq)
                if "tag" in codon_list or "tga" in codon_list or "taa" in codon_list:
                    warnings.warn("DNA sequence contains stop codon in ORF")
                    stopcod_OR_ambignts.append((ProID, Acc, codon_list, aaSeq)) #will only contain ambig nts for codon_other
                    hadstopcod += 1
                    continue  #do not run seg on it

                #get the sequences and variables to write to Seg
                if direction == "pro_dna" or direction == "pro_cod": #Pro_other
                    body = aaSeq
                    d = "prot"
                elif direction == "dna_pro" or direction == "dna_cod": #DNA_other
                    body = GeneSeq
                    d = "dna"
                elif direction == "cod_pro" or direction == "cod_dna":
                    CodonSeq = convert_to_ascii(ProID, Acc, codon_list) 
                    body = CodonSeq
                    d = "codon"
                    if CodonSeq == "NA":
                        stopcod_OR_ambignts.append((ProID, Acc, codon_list, aaSeq))
                        continue #cannot convert to ascii, cannot run through seg
                
                #random file generator to write to fasta file to feed seg
                tfName1 = tempfile.NamedTemporaryFile(delete = False).name
                header = ">" + ProID + "\n"
                with open(tfName1, "w") as fh:
                    fh.write(header)
                    fh.write(body)   #body is aaSeq, GeneSeq, or CodonSeq
              
                #run SEG on protein and pass all LCRs into new Fasta format file
                SEG = subprocess.check_output(["segA", tfName1, W, K1, K2, d, "-l", "-t", trim]) #default alphabet size = 20
                subprocess.call(["rm", tfName1])   #delete tempfile1
                tfName2 = tempfile.NamedTemporaryFile(delete = False).name  #make new tempfile2
                with open(tfName2, "w") as fh:
                    fh.write(SEG.decode())
              
                #find longest LCR 
                if AllLCR == True:   #if we need to take all the LCRs and not just 1
                    pass 
                #take only the longest LCR
                elif AllLCR == False:  
                    longestlcr_pos = LongestLCR(tfName2, KeepEnds)  #usually it will be a list of 1 tup
                    subprocess.call(["rm", tfName2])  #remove tempfile2

                    if longestlcr_pos == "noLCR":
                        NoLCR[nolcr_count] = [GeneID, ProID, Acc, "noLCR", Product, Description, aaSeq]
                        nolcr_count += 1
                        continue #this continue is very necessary
                                 #do not put it into the AllGeneIDs because this way

                    elif len(longestlcr_pos) > 1: #if LCRs are the same length, take the one with the lowest Entropy
                        #these variables are needed to give entropy calculator proper sequence + sequence type name
                        if direction == "pro_dna" or direction == "pro_cod":
                            seq = aaSeq
                            Type = "protein"
                        elif direction == "dna_pro" or direction == "dna_cod":
                            seq = GeneSeq  
                            Type = "dna"
                        elif direction == "cod_pro" or direction == "cod_dna":
                            seq = GeneSeq #entropy calculator uses DNAseq instead of codon seq
                            Type = "codon"
                        #now compare entropies of the LCR sequences from seg
                        EntCompare = {}
                        for tup in longestlcr_pos:
                            s = tup[0]-1 #where to index sequence to get the LCR
                            e = tup[1]
                            if direction == "cod_pro" or direction == "cod_dna":
                                s = tup[0]*3-3 #because positions are based on ascii chars, but need to index DNA
                                e = tup[1]*3   #for 0 index splicing
                            EntCompare[tup] = Entropy(ProID, Acc, seq[s:e], Type)
                        Min = min(list(EntCompare.values()))
                        for key, val in EntCompare.items():
                            if Min == val:
                                longestlcr_pos = [key]  #put tuple in list to keep the same as below
                        
                    #else: # only 1 tuple was returned for longestlcr_pos
                    #longestlcr_pos could be aa_pos or gene_pos (aa_pos if direction is pro_dna or pro_cod)
                    #gene_pos for dna_pro,  ascii_pos  (like aa_pos) for cod_pro
                    LCR_chrom_pos = MapToChrom(direction, strand, longestlcr_pos, GeneStart, GeneEnd, intron, exonPos)  

                    #get substring of DNA to feed to Entropy() to get Ed or Ec
                    if len(LCR_chrom_pos) == 1:  #only one dict entry, with tuple of positions means the LCR does not span and exon junction
                        spanJunc = "-"    #spanjuncvariables will be written to file
                        for pos, stra in LCR_chrom_pos.items(): #will only do 1 iteration
                            LCRDNASeq = full_chrom_seq[pos[0]-1:pos[1]]  
                            if stra == "(-)":
                                LCRDNASeq = str(Seq(LCRDNASeq).reverse_complement())  #used to calculate entropy DNA

                    #if the LCR spans an exon junction (LCR_chrom_pos[pos] is a tuple)
                    else: # len(LCR_chrom_pos) > 1:
                        spanJunc = "+"
                        LCRDNASeq = ""
                        for pos, stra in LCR_chrom_pos.items():
                            if stra == "(+)":                                   
                                LCRDNASeq += full_chrom_seq[pos[0]-1:pos[1]]
                            elif stra == "(-)":
                                s = full_chrom_seq[pos[0]-1:pos[1]]
                                LCRDNASeq += str(Seq(s).reverse_complement())
                    #get LCRAAseq
                    if direction == "pro_dna" or direction == "pro_cod"\
                        or direction == "cod_pro" or direction == "cod_dna":  #get LCR_aa_pos and LCRDNAseq
                        LCR_aa_pos = longestlcr_pos                          #for cod_pro, ascii positions == aa_pos
                    elif direction == "dna_pro" or direction == "dna_cod":
                        LCR_aa_pos = [(MapToProtein(longestlcr_pos))]  #longestlcr_pos is LCR_gene_pos #must be a list of tups to match

                    LCRAASeq = aaSeq[LCR_aa_pos[0][0]-1:LCR_aa_pos[0][1]]
                    #now determine entropies! based on direction, so unnecessary calcs aren't done
                    if direction == "pro_dna" or direction == "dna_pro":
                        Ep = Entropy(ProID, Acc, LCRAASeq, "protein")
                        Ed = Entropy(ProID, Acc, LCRDNASeq, "dna")
                        x = Ep     #x and y are variables to write to file
                        y = Ed
                    elif direction == "pro_cod" or direction == "cod_pro":
                        Ep = Entropy(ProID, Acc, LCRAASeq, "protein")
                        Ec = Entropy(ProID, Acc, LCRDNASeq, "codon")
                        x = Ep    #x and y are variables to write to file
                        y = Ec
                    elif direction == "cod_dna" or direction == "dna_cod":
                        Ec = Entropy(ProID, Acc, LCRDNASeq, "codon")  #use DNA seq to get codon ent (sequences cod_other will not contain ambig chars)
                        Ed = Entropy(ProID, Acc, LCRDNASeq, "dna")                                  #sequences other_codon may contain ambig chars
                        x = Ec
                        y = Ed

                    #get some variables for file writing!
                    S = []; E = []
                    for pos in LCR_chrom_pos.keys():
                        S.append(pos[0])
                        E.append(pos[1])
                    LCR_chrom_s = min(S)   #where the LCR starts and ends on the chromosome 
                    LCR_chrom_e = max(E)

                    LCR_aa_s = LCR_aa_pos[0][0]  #where the LCR starts and ends in the protein
                    LCR_aa_e = LCR_aa_pos[0][1]  #is equivalent to which codon the LCR starts and ends at

                    #add all the info to the dict                        

                    replacedORAdded += 1 #for counting 
                    AllGeneIDs[GeneID] = [ProID, Acc, strand, intron, spanJunc, GeneStart, \
                                         GeneEnd, LCR_aa_s, LCR_aa_e, x, LCR_chrom_s, LCR_chrom_e, y,\
                                         Product, Description, exonPos, LCR_chrom_pos, LCRDNASeq,\
                                         LCRAASeq, aaSeq]
                                         #must add aaSeq at [-1] position so in the next iteration, length can be compared

    #all the info is collected for 1 instance of each gene, so now sort through it and write those files.
    print("data collected, writing to files")
    if direction == "pro_dna":
        n1 = "#Protein to DNA\n"
        x1 = "EntropyProtein"
        y1 = "EntropyDNA"
    elif direction == "dna_pro":
        n1 = "#DNA to Protein\n"
        x1 = "EntropyProtein"
        y1 = "EntropyDNA"
    elif direction == "pro_cod":
        n1 = "#Protein to Codon"
        x1 = "EntropyProtein"
        y1 = "EntropyCodon"
    elif direction == "cod_pro":
        n1 = "#Codon to Protein"
        x1 = "EntropyProtein"
        y1 = "EntropyCodon"
    elif direction == "dna_cod":
        n1 = "#DNA to Codon"
        x1 = "EntropyCodon"
        y1 = "EntropyDNA"
    elif direction == "cod_dna":
        n1 = "#Codon to DNA"
        x1 = "EntropyCodon"
        y1 = "EntropyDNA"

    Note = n1 + "#Organism: " + Organism + "\n" + "#seg parameters used for W, K1, K2 were "\
            + W + ", " + K1 + ", " + K2 + ", respecively" + "\n"
    filename = Organism + W + K1 + K2 + "_"+ direction +"_" 

    #write the headers of the file
    HEAD1 = Note + "#This file contains only proteins with no LCRs.\
    To confirm the file has been parsed properly\n"
    with open(filename + "NoLCRs", "w") as fh:
        fh.write(HEAD1)
    HEAD2 = Note + "#This file contains only locations of LCRs which span 1 or more exon junctions.\
    Numbers are in pairs.\nProID\tAccession\tstrand\tLCRpositions\n"
    with open(filename + "multexon", "w") as fh:
        fh.write(HEAD2)
    HEAD3 = Note + "#dup_replaced: " + str(dup-dupnotreplaced) + "\treplacedORadded: "\
    + str(replacedORAdded-hadstopcod) + "\t"+"noLCRs: "+ str(nolcr_count) + "\ttotal: "+str(total)+"\n"\
    "ProID\t\taccession\tstrand\tintron\tJunc\tGeneS\
    \tGeneE\tLCRProS\tLCRProE\t"+x1+"\tLCRChromS\tLCRChromE\t"+y1+"\n" 
    with open(filename+ "Output", "w") as fh:
        fh.write(HEAD3)
    HEAD4 = Note + "ProID\tAcc\tProduct\tDescription\n"
    with open(filename + "Output_Des", "w") as fh:
        fh.write(HEAD4)
    HEAD5 = Note + "#This file contains the LCR DNA and AA sequences with ProID and Acc as a header\n"
    with open(filename + "Seq", "w") as fh:
        fh.write(HEAD5)
    HEAD6 = Note + "#This file contains all of the duplicate gene IDs that\
    were not used for LCR identification. Count them up to make sure\
    the file has been parsed properly\nGeneID\tProID\tAccession\n"
    with open(filename + "duplicateGeneIDs", "w") as fh:
        fh.write(HEAD6)
    HEAD7 = Note + "#This file contains DNA sequences which had stop codons in their ORF\
    \n#If direction is codon_other, this file also contains sequences with ambig nts which could not be converted to ascii.\n\
    #Uses ProID and Acc as header\n\n"
    with open(filename + "StopCodon", "w") as fh:
        fh.write(HEAD7)
    HEAD8 = Note + "#This file contains DNA + protein sequences in which the DNA seqeunce was not 3X the protein sequence.\n"
    with open(filename + "UnequalLengths", "w") as fh:
        fh.write(HEAD8)

    #add the info to the files
    INFO1 = ""
    for key in NoLCR:
        ProID = NoLCR[key][1]
        Acc = NoLCR[key][2]
        INFO1 += ProID+"\t"+Acc+"\tNoLCR\n"
    with open(filename + "NoLCRs", "a") as fh:
        fh.write(INFO1)

    INFO2 = ""; INFO3 = ""; INFO4 = ""; INFO5 = ""
    for key in AllGeneIDs:
        
        ProID = AllGeneIDs[key][0]
        Acc = AllGeneIDs[key][1]
        strand = AllGeneIDs[key][2]    #index all variables from dictionary
        intron = AllGeneIDs[key][3]
        spanJunc  = AllGeneIDs[key][4]
        GeneStart = AllGeneIDs[key][5]
        GeneEnd = AllGeneIDs[key][6]
        LCR_aa_s = AllGeneIDs[key][7]
        LCR_aa_e = AllGeneIDs[key][8]
        x = AllGeneIDs[key][9]             #EntropyProtein or #Entropy Codon
        LCR_chrom_s = AllGeneIDs[key][10]
        LCR_chrom_e = AllGeneIDs[key][11]
        y = AllGeneIDs[key][12]            #EntropyDNA or EntropyCodon (depending on direction)
        Product = AllGeneIDs[key][13]
        Description = AllGeneIDs[key][14]
        exonPos = AllGeneIDs[key][15]
        LCR_chrom_pos = AllGeneIDs[key][16]
        LCRDNASeq = AllGeneIDs[key][17]   
        LCRAASeq = AllGeneIDs[key][18]

        if spanJunc == "+":  #raise some warnings and write spanjunc file
            warnings.warn("LCR spans 1 or more exon junctions")
            pos = ""
            for tup, st in LCR_chrom_pos.items():
                pos += str(tup[0]) + st + "\t"
                pos += str(tup[1]) + st + "\t"
            plusminusLCR = [] #this is to raise a warning if LCR is on different strands
            for char in pos:
                if char == "+" or char == "-":
                    plusminusLCR.append(char)
            if "+" in plusminusLCR and "-" in plusminusLCR:
                warnings.warn(ProID + " " + Acc + " " + str(strand) + " LCR is on positive and negative strand")
            
            INFO2 += ProID + "\t" + Acc + "\t" + str(strand) + "\t" + pos + "\n" #write the spanjunc+ exon file

        INFO3 += ProID + "\t" + Acc + "\t" + str(strand) + "\t" + str(intron) + "\t" + spanJunc\
        + "\t" + str(GeneStart) + "\t" + str(GeneEnd) + "\t" + str(LCR_aa_s)  + "\t"\
        + str(LCR_aa_e) + "\t" + str(x) + "\t" + str(LCR_chrom_s) + "\t" + str(LCR_chrom_e)\
        + "\t" + str(y) + "\n"

        INFO4 += ProID + "\t" + Acc + "\t" + Product + "\t" + Description + "\n"

        INFO5 += ">"+ProID+"_"+Acc+"\n"+LCRDNASeq+"\n\n"+">"+ ProID+"_"+Acc\
        +"\n"+LCRAASeq+"\n\n\n"

    
    with open(filename + "NoLCRs", "a") as fh:
        fh.write(INFO1)
    with open(filename + "multexon", "a") as fh:
        fh.write(INFO2)
    with open(filename + "Output", "a") as fh:
        fh.write(INFO3)
    with open(filename + "Output_Des", "a") as fh:
        fh.write(INFO4)
    with open(filename + "Seq", "a") as fh:
        fh.write(INFO5)

    INFO6 = ""
    for val in ExtraGeneIDs.values():
        INFO6 += val[0] + "\t" + val[1] + "\t" + val[2] + "\n"
    with open(filename + "duplicateGeneIDs", "a") as fh:
        fh.write(INFO6)

    INFO7 = ""
    for tup in stopcod_OR_ambignts:
        INFO7 += ">"+tup[0]+"_"+tup[1]+"\n"+str(tup[2])+"\n\n"+\
                 ">"+tup[0]+"_"+tup[1]+"\n"+tup[3]+"\n\n"
    with open(filename + "StopCodon", "a") as fh:
        fh.write(INFO7)

    INFO8 = ""
    for tup in unequal_lengths:
        INFO8 += ">"+tup[0]+"_"+tup[1]+"_length"+str(tup[4])+"\n"+tup[2]+"\n\n"+\
                 ">"+tup[0]+"_"+tup[1]+"_length"+str(tup[5])+"\n"+tup[3]+"\n\n\n"
    with open(filename + "UnequalLengths", "a") as fh:
        fh.write(INFO8)


    print("total", total)
    print("replaced or added ", replacedORAdded-hadstopcod)
    print("duplicate ", dup)
    print("duplicate not replaced", dupnotreplaced)
    print("has stop codon", hadstopcod)

#write unequal Pos file

CombinedSysCall(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
#CombinedSysCall(sys.argv[1], sys.argv[2],  sys.argv[3], sys.argv[4], sys.argv[5])
#CombinedSysCall(input("gb file name: "), input("direction('DNA' or 'protein'): "), input("W: "), input("K1: "), input("K2: ")\
#,input("trim(d): "), input("AllLCR(d): "), input("KeepEnds(d): "))

#Pro_DNA(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

#directions:   pro_dna  dna_pro  pro_cod



#Pro_DNA("proj1/S.cer_gen_all", "S.cer_Output/100trim_Data/Pro_DNA/S.cer152225", "S.cer_Output/100trim_Data/Pro_DNA/S.cer152225_Des", "15", "2.2","2.5", input(), input())


#       fileR      fileWI               fileWD                W         K1      K2 
# trim(d)    AllLCR(d)

#"/scratch/Human/humanchrom_all"
# "HumfileW_ProtoDNA", "HumGeneDes_ProtoDNA"
#was originally Pro_DNA
