#!/usr/bin/python3
from sys import argv as a
from re import search
from math import log


def Split_Into_Codons(DNAseq):
    codon_list = []
    for i in range(0, len(DNAseq)-2, 3):
        codon = DNAseq[i] + DNAseq[i+1] + DNAseq[i+2]
        codon_list.append(codon)
    return codon_list


def Entropy(seqString, Type):

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
                #print(codon_list)
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


#this takes a dna string and searches for patterns (mono-nt repeat, di-nt repeat, tri-nt repeat, --> 9 nt repeat)
def searchfordnapattern(string, rep1, rep2, rep3, rep4, rep5, rep6, rep7, rep8, rep9):

    pat1 = r"(\w)\1{"+rep1+",}" 
    pat2 = r"(\w)(\w)(\1\2){"+rep2+",}"         #a 2-letter tandem repeat
    pat3 = r"(\w)(\w)(\w)(\1\2\3){"+rep3+",}"    #a 3-letter tandem repeat
    pat4 = r"(\w)(\w)(\w)(\w)(\1\2\3\4){"+rep4+",}"
    pat5 = r"(\w)(\w)(\w)(\w)(\w)(\1\2\3\4\5){"+rep5+",}"
    pat6 = r"(\w)(\w)(\w)(\w)(\w)(\w)(\1\2\3\4\5\6){"+rep6+",}"
    pat7 = r"(\w)(\w)(\w)(\w)(\w)(\w)(\w)(\1\2\3\4\5\6\7){"+rep7+",}"
    pat8 = r"(\w)(\w)(\w)(\w)(\w)(\w)(\w)(\w)(\1\2\3\4\5\6\7\8){"+rep8+",}"
    pat9 = r"(\w)(\w)(\w)(\w)(\w)(\w)(\w)(\w)(\w)(\1\2\3\4\5\6\7\8\9){"+rep9+",}"

    spat1 = search(pat1, string)
    spat2 = search(pat2, string)
    spat3 = search(pat3, string)
    spat4 = search(pat4, string)
    spat5 = search(pat5, string)
    spat6 = search(pat6, string)    
    spat7 = search(pat7, string)
    spat8 = search(pat8, string)
    spat9 = search(pat9, string)

    if not spat1 and not spat2 and not spat3 and not spat4\
    and not spat5 and not spat6 and not spat7 and not spat8 and not spat9:
        return False
    if spat1:
        spat = spat1.group(0)
        pattype = 1
    elif spat3:
        spat = spat3.group(0)
        pattype = 3   
    elif spat6:
        spat = spat6.group(0)
        pattype = 6   
    elif spat9:
        spat = spat9.group(0)
        pattype = 9
    elif spat2:
        spat = spat2.group(0)
        pattype = 2
    elif spat4:
        spat = spat4.group(0)
        pattype = 4
    elif spat5:
        spat = spat5.group(0)
        pattype = 5
    elif spat7:
        spat = spat7.group(0)
        pattype = 7
    elif spat8:
        spat = spat8.group(0)
        pattype = 8


    lowB = string.index(spat)+1  #1 indexed
    highB = lowB + len(spat)-1     #1 indexed
    return(pattype, spat, lowB, highB)           #maybe I should find a way to keep track of if there is multiple?

#this function searches for a pattern within a string
#assumption made that the pattern only contains one such string.
#if it has a single letter repeat, this will take precedence over di-repeat, which will take precedence over tri-repeat (etc.)
def searchforpropattern(string, rep1, rep2, rep3):

    pat1 = r"(\w)\1{"+rep1+",}"           #a single letter repeat
    pat2 = r"(\w)(\w)(\1\2){"+rep2+",}"         #a 2-letter tandem repeat
    pat3 = r"(\w)(\w)(\w)(\1\2\3){"+rep3+",}"    #a 3-letter tandem repeat
    
    spat1 = search(pat1, string)
    spat2 = search(pat2, string)
    spat3 = search(pat3, string)

    if not spat1 and not spat2 and not spat3:
        return False
    if spat1:
        spat = spat1.group(0)
        pattype = 1
      #  print(spat)
    elif spat2:
        spat = spat2.group(0)
        pattype = 2
      #  print(spat)
    elif spat3:
        spat = spat3.group(0)
        pattype = 3
      #  print(spat)
    lowB = string.index(spat)+1  #1 indexed
    highB = lowB + len(spat)-1     #1 indexed
   # print(lowB, highB)
    #print(string)
    return(pattype, spat, lowB, highB)

#trims ends of DNA to match protein
def trim_dna_ends(a, b):
    if (a-1)%3 == 0 and b%3 == 0: #both are fine (first and third)
        pass
    elif (a-1)%3 == 0 and (b+1)%3 == 0: #start in 1st pos, end in 2nd pos
        b = b-2
    elif (a-1)%3 == 0 and (b+2)%3 == 0: #start in 1st pos, end is in 1st pos
        b = b-1
    elif (a-2)%3 == 0 and b%3 == 0: #start in 2nd pos, end in 3rd pos
        a = a+2
    elif (a-2)%3 == 0 and (b+1)%3 == 0: #start in 2nd pos, end in 2nd pos
        a = a+2
        b = b-2
    elif (a-2)%3 == 0 and (b-1)%3 == 0: #start in 2nd pos, end in 1st pos
        a = a+2
        b = b-1
    elif a%3 == 0 and b%3 == 0: #start in 3rd pos, end in 3rd pos
        a = a+1
    elif a%3 == 0 and (b+1)%3 == 0: #start in 3rd pos, end in 2nd pos
        a = a+1
        b = b-2
    else:  #start in 3rd pos, end in 1st pos
        a = a+1
        b = b-1
    return (a,b)

#takes positions of a pattern within a protein sequence (as 1 indexed)
#returns the start and end positions of the corresponding DNA sequence
def maptodna(prostart, proend):
    dnastart = prostart*3 - 2
    dnaend = proend*3
    return (dnastart, dnaend) #1 indexed

def maptopro(dnastart, dnaend):
    d = {}
    dnastart = trim_dna_ends(dnastart,dnaend)[0]
    dnaend = trim_dna_ends(dnastart,dnaend)[1]
    prostart = ((dnastart-1)//3)+1
    proend = ((dnaend-1)//3)+1
     
    d["dna"]=(dnastart, dnaend) ; d["protein"]=(prostart, proend) #1 indexed
    return d

#this is the main function
#takes the file with sequences, searches the sequences for periodicity
#returns a file with the corresonding protein, codon, and dna entropy --> for either the whole sequence (True) or just the periodic part (False)
#output file: col1 will be the id, col2 will be repeat type (1,2, or 3), col3 will be dna ent, col4 will be codon ent, col5 will be protein ent
#output file will have header of the inputed seuqence filename, as well as the min number of tandem repeats which were searched for.

def periodicity(filename1, whichseq, wholeseq="True", rep1p="4", rep2p="2", rep3p="2", rep1d="11", rep2d="5", rep3d="4", rep4d="2", rep5d="2", rep6d="2", rep7d="2", rep8d="2", rep9d="2"):

    print("filename1, whichseq, wholeseq='True', rep1p='4', rep2p='2', rep3p='2', rep1d='11', rep2d='5', rep3d='4', rep4d='2', rep5d='2', rep6d='2', rep7d='2', rep8d='2', rep9d='2'")
    if whichseq != "protein" and whichseq != "dna":
        print("whichseq (2nd arg) must == 'protein' or 'dna'")
    if wholeseq != "True" and wholeseq != "False":
        print("wholeseq (3rd arg) must == True or False")
  #  if type(eval(rep1p)) != int or type(eval(rep2p)) != int or type(eval(rep3p)) != int or type(eval(rep1d)) or type(eval(rep2d)) != int or type(eval(rep3d)) != int or\
  #  type(eval(rep4d)) != int or type(eval(rep5d)) != int or type(eval(rep6d)) != int or type(eval(rep7d)) or type(eval(rep8d)) != int or type(eval(rep9d)) != int:
  #      print("rep must be an integer (value is 1 less than the minimum number of repeats you are looking for)")

    entinfo_per = {}
    entinfo_noper = {}
    with open(filename1, "r") as fh:
        f = fh.readlines()[4:]

        for line in range(len(f)):
            if line % 7 == 0:
                header = f[line][:-1]
                dnaseq = f[line+1][:-1]   #index to :-1 to get rid of the "\n" from the string
                proseq = f[line+4][:-1]
             #   print(proseq)
             #   print(dnaseq)
                #print(header)
                if whichseq == "protein":
                    pattern = searchforpropattern(proseq, rep1p, rep2p, rep3p)  #returns False or (pattern_type, pattern_seq, lowbound, highbound)
                elif whichseq == "dna":
                    pattern = searchfordnapattern(dnaseq, rep1d, rep2d, rep3d, rep4d, rep5d, rep6d, rep7d, rep8d, rep9d)

                if pattern == False:  #if pattern is false, should also take these sequences with no periodicity and do comparative entropy test.
                    dnaEnt = Entropy(dnaseq, "dna")
                    codEnt = Entropy(dnaseq, "codon")
                    proEnt = Entropy(proseq, "protein")
                    entinfo_noper[header] = [dnaEnt, codEnt, proEnt]
                else:   # there is a pattern, so get entropy 
                    pattype = pattern[0]
                    patseq = pattern[1]
                    lowB = pattern[2]
                    highB = pattern[3] 
                        
                    if wholeseq == "True":
                        dnaEnt = Entropy(dnaseq, "dna")
                        codEnt = Entropy(dnaseq, "codon")
                        proEnt = Entropy(proseq, "protein")

                        print(dnaseq)
                        print(patseq)
                        print(pattype, lowB, highB)
                        print(proseq, end="\n\n\n")

                    elif wholeseq == "False":
                        if whichseq == "protein":
                            dnapos = maptodna(lowB, highB)         #just some stuff to get the corresponding DNA seq 
                            dnastart = dnapos[0]-1  #0 indexed
                            dnaend = dnapos[1]

                            dnaEnt = Entropy(dnaseq[dnastart:dnaend], "dna")
                            codEnt = Entropy(dnaseq[dnastart:dnaend], "codon")
                            proEnt = Entropy(patseq, "protein")

                        elif whichseq == "dna":
                            dict_pos = maptopro(lowB, highB)
                            lowB = dict_pos["dna"][0]-1  #dna start  #0 indexed
                            highB = dict_pos["dna"][1] #dna end
                            prostart = dict_pos["protein"][0]-1   #0 indexed
                            proend = dict_pos["protein"][1]

                            dnaEnt = Entropy(dnaseq[lowB:highB], "dna")
                            codEnt = Entropy(dnaseq[lowB:highB], "codon")
                            proEnt = Entropy(proseq[prostart:proend], "protein")

                    entinfo_per[header] = [pattype, lowB, highB, dnaEnt, codEnt, proEnt]


    #now write this information to a file
    
    if wholeseq == "True":
        w = "wholeseq_"
        s = "The entropy values were taken from the WHOLE LCR sequence so long as it contained a repeat."
    elif wholeseq == "False":
        w = "partseq_"
        s = "The entropy values were taken from ONLY the portion of the LCR with a periodic repeat." 
    if whichseq == "protein":
        NOTE1 = "#this entropy information came from the sequence file: "+filename1+"\n#"\
        +s+"\n#Minimum repeat values for mono-, di-, and tri- amino acid repeats were "+str(int(rep1p)+1)+", "+str(int(rep2p)+1)+", and "+str(int(rep3p)+1)+", respectively.\n"
    elif whichseq == "dna":
        NOTE1 = "#this entropy information came from the sequence file: "+filename1+"\n#"\
        +s+"\n#Minimum repeat values for mono- to nine nucleotide repeats were "+str(int(rep1d)+1)+", "+str(int(rep2d)+1)+", "+str(int(rep3d)+1)+\
        str(int(rep4d)+1)+", "+str(int(rep5d)+1)+", "+str(int(rep6d)+1)+", "+str(int(rep7d)+1)+", "+str(int(rep8d)+1)+", and "+str(int(rep9d)+1)+", respectively.\n"

    header1 = "ID\t\t\tpattype\tlowB\thighB\tdnaEnt\tcodEnt\tproEnt\n"
   # print(to_file)
    to_file1 = NOTE1 + header1
    for k, v in entinfo_per.items():
        id = k
        pattype = str(v[0])
        lowB = str(v[1])
        highB = str(v[2])
        dnaEnt = str(v[3])
        codEnt = str(v[4])
        proEnt = str(v[5])
        to_file1 += id+"\t"+pattype+"\t"+lowB+"\t"+highB+"\t"+dnaEnt+"\t"+codEnt+"\t"+proEnt+"\n"

    NOTE2 = "#this entropy information came from the sequence file: "+filename1+"\n#"\
    "This file contains the entropy information for the sequences which did NOT have periodic repeats\n\
    Repeats were taken from the"+whichseq+"sequence"
    header2 = "ID\t\t\tdnaEnt\tcodEnt\tproEnt\n"
    to_file2 = NOTE2 + header2
    for k, v in entinfo_noper.items():
        id = k
        dnaEnt = str(v[0])
        codEnt = str(v[1])
        proEnt = str(v[2])
        to_file2 += id+"\t"+dnaEnt+"\t"+codEnt+"\t"+proEnt+"\n"
#    print(to_file)

    #generating filename becuase I hate creating filenames.
    namestart = filename1.rindex("/")+1
    nameend = filename1.rindex("_")
    if whichseq == "protein":
        rep = rep1p+rep2P+rep3p
    elif whichseq == "dna":
        rep = rep1d+rep2d+rep3d+rep4d+rep5d+rep6d+rep7d+rep8d+rep9d

    filename2 = "periodicity_" + whichseq + w + rep +"_" + filename1[namestart:nameend] + "EntOutput"
   # print(filename2)
    with open(filename2, "w") as fh:  #file containing periodic repeats
        fh.write(to_file1)
    filename3 = "noperiodicity_"+ rep +"_"+filename1[namestart:nameend]+"EntOutput"  #file not containing periodic repeats
    with open(filename3, "w") as fh:
        fh.write(to_file2)

###THE END####

#periodicity("~/proj1/S.cer_Output/DNA_Pro/ran0722/Seq/Saccharomyces_cerevisiae451.31.5_dna_pro_Seq", "dna")
periodicity(a[1], a[2])
#try:
#    periodicity(a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12])
#except IndexError:
#    try:
#        periodicity(a[1], a[2], a[3], a[4], a[5])
#    except IndexError:
#        try:
#            periodicity(a[0], a[1], a[2])
#        except IndexError:
#            try:
#                periodicity(a[0], a[1])
#            except IndexError:
#                print("put seq filename, whichseq ('protein' or 'dna') wholeseq (True or False), repeatnum1 --> 3 OR 9 (depending)\n\
#                Repeat numbers are 1 less than the minimum repeat search length.\n\
#                Put in order of mono-, di-, tri-, etc. repeat number.")

