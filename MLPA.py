from math import log                                            #import of modules
import operator
import os

zurücksetzen = open('blast_LHS_RHS.fasta.txt', 'w')             #reset of text-files
zurücksetzen6 = open("blast_sonden.txt", "w")                   #reset of text-files
zurücksetzen7 = open("blast_sonden2.txt", "w")                  #reset of text-files
zurücksetzen8 = open("blast_result_final.txt", "w")             #reset of text-files
zurücksetzen9 = open("auswahl.txt", "w")                        #reset of text-files

laengen = [20,21,22,23,24,25,26,27,28,29,30]                    #lenghts of HS used for generating HS
stuffer_lhs = "cgctactactattagtagaattgatgccaccttttcagctcgcg"    #sequence of stuffer sequence of LHS
stuffer_rhs = "accatttgcgaaatgtatctaatggtcaaactaaatctac"        #sequence of stuffer sequence of RHS
lhs_primer = "GGGTTCCCTAAGGGTTGGA"                              #primer sequence LHS
rhs_primer = "TCTAGATTGGATCTTGCTGGCAC"                          #primer sequence RHS
p300_laengen = [64,82,88,92,108,130,148,173,180,184,191,196,208,214,220,226,232,239,246,252,258,265,274]    #probe lengths of P300-Probemix


def reading_input():                                            #reading the sequences (in FASTA-format) of input.txt and saving them in a list
    with open("input.txt", "r") as input_file:
        namen = []
        sequenzen = []
        for zeile in input_file:
            zeile = zeile.rstrip()
            if zeile[0] == ">":                                 #if line starts with a '>' it is an indication that the line only contains sequence name/ exon name
                namen.append(zeile[1:])
            else:
                sequenzen.append(zeile.upper())                 #otherwise it is the sequence
        exone = []
        for i in range(0, len(namen)):
            temp = []
            temp.append(namen[i])
            temp.append(sequenzen[i])
            exone.append(temp)                                   #saving name and sequence of FASTA in a list
        for j in exone:
            name = "sonden_" + j[0]                              #text-file for each sequence in input.txt is created
            with open(name, "w+") as filename:
                file = filename
    return exone                                                 #output of method is the list containing sequence and name

def Tm_Value(hs):                                                #calculation of the Tm-value of HS
    nearest_neighbour = open("nearest_neighbour_values.txt", "r")   #text-file containing ΔH- and ΔS-values for base interactions
    nearest_neighbour_values = {}
    for zeile in nearest_neighbour:
        zeile = zeile.rstrip()
        liste = zeile.split("\t")
        basen = liste[0]
        delta_werte = []
        delta_werte.append(float(liste[1]))
        delta_werte.append(float(liste[2]))
        nearest_neighbour_values[basen] = delta_werte               #saving Δ-values of text-file in a list
    delta_H = 0.0
    delta_S = 0.0
    for i in range(0,len(hs)):
        neighbour = hs[i:i+2]
        if neighbour in nearest_neighbour_values:
            delta_H = nearest_neighbour_values[neighbour][0] + delta_H      #determine ΔH of neighbouring bases
            delta_S = nearest_neighbour_values[neighbour][1] + delta_S      #determine ΔS of neighbouring bases
    Tm = delta_H/(-0.0108+delta_S+0.00199*log(0.0000008/4))-273.15+log(0.035,16.6)-21   #formula for Tm
    return round(Tm)                                              #return of method is rounded Tm-value of HS


def primer_bs_neu(kandidat):                                      #calculating if primer binding region of hs falls under criteria
    kandidat_anfang =kandidat[0:3]                                #first 3 bases of 5' end from HS
    kandidat_ende = kandidat[len(kandidat)-3:len(kandidat)]       #first 3 bases of 3' end from HS
    gc_anfang = 0
    gc_ende = 0
    for i in kandidat_anfang:
        if i == "G" or i == "C":
            gc_anfang=gc_anfang+1                               #counting GC
    for j in kandidat_ende:
        if j == "G" or j == "C":
            gc_ende = gc_ende +1                                #counting GC
    if gc_anfang <=3 and gc_ende<=3:                            #if less than 3 GC at the 3 bases of 5' and 3' end
        return True                                             #return positive boolean

def ligation(kandidat,plaenge):                                 #determines if ligation region of LHS falls under criteria
    lhs = kandidat[0:plaenge]
    ligation = 0
    for j in lhs[len(lhs) - 5:len(lhs)]:                        #last 5 bases of LHS
        if j == "G" or j == "C":
            ligation = ligation + 1
    return ligation                                             #return amount of GC in last 5 bases of LHS at primer binding site

def gc_gehalt(kandidat,plaenge):                                #calculation of GC-content of LHS and RHS
    gc_lhs = 0
    gc_rhs = 0
    lhs = kandidat[0:plaenge]
    rhs = kandidat[plaenge:]
    gc_gesamt = []
    for k in range(0, len(lhs)):
        if lhs[k] == "G" or lhs[k] == "C":
            gc_lhs = gc_lhs + 1
        if rhs[k] == "G" or rhs[k] == "C":
            gc_rhs = gc_rhs + 1
    gc_gesamt.append(round(gc_lhs/len(lhs)*100))                #converting GC-content in precent and rounding
    gc_gesamt.append(round(gc_rhs/len(rhs)*100))
    return gc_gesamt                                            #return of a list containing GC-content of LHS and RHS

def gc_rek(hs,side,plaenge):                                    #GC-content calculation of LHS/ RHS with puffer length attached
    gc_liste = []
    puffer_laenge = len(hs)-plaenge
    maximum = 100
    ideal = 0
    if puffer_laenge == 0:
        return (hs, 0)
    for i in range(0,puffer_laenge):
        gc_hs = 0
        if side == "l":                                             #GC-content calculation of LHS
            for j in range(i, len(hs)):
                if hs[j] == "G" or hs[j] == "C":
                    gc_hs = gc_hs + 1
            gc_liste.append(gc_hs / len(hs[i:len(hs)]) * 100)
        if side == "r":                                             #GC-content calculation of RHS
            for j in range(0, len(hs)-i):
                if hs[j] == "G" or hs[j] == "C":
                    gc_hs = gc_hs + 1
            gc_liste.append(gc_hs / len(hs[0:len(hs)-i]) * 100)     #list with GC-content of HS plus every puffer length combination

    for i in sorted(gc_liste):                                      #sorts list by GC-content closest to be at 50%
        if abs(i-50) < maximum:
            maximum = abs(i-50)
            ideal = i
    if side=="l":
        maximum = hs[gc_liste.index(ideal):]
        return (maximum, round(ideal))                              #return of LHS plus puffer GC-content
    if side=="r":
        maximum = hs[:len(hs)-gc_liste.index(ideal)]
        return (maximum, round(ideal))                              #return of RHS plus puffer GC-content

def identifikation(pexon):                                          #identifies exon of calculated HS
    exone = reading_input()

    for i in exone:
        if pexon==i[1]:
            name =  i[0]
            file = "sonden_"+name+".txt"
            with open(file, "a") as sonden:
                ausgabe = [name,sonden]
                return ausgabe                                      #returns Exon number and text-file name of Exon

def local_blast(pdatei, pliste_sonden):                             #compares text-files of every exon containig generated HS to each other and deletes HS with an alignment of over 12 bases
    with open(pdatei, "a") as datei:
        liste_sonden = pliste_sonden                                #list containing text-files of HS for every exon including P300-probemix sequences
        liste_used = {}
        file_iterator = 0

        for i in liste_sonden:  # For every file                    #HS of one exon only compared to HS of other exons not to same exon
            file_iterator = file_iterator + 1
            with open(i.name, "r") as sonde:  # sonde = current file
                for j in sonde:  # for every line j in file sonde
                    zeile = j.rstrip().split("\t")  # split line in array
                    hs = zeile[1]                                   #HS-seqeunce of exon X
                    hs_daten = zeile[0]
                    isgud = True

                    for n in range(0, len(hs)):
                        if n + 11 > len(hs): break
                        check = hs[n:n + 12]                        #goes through HS and takes the next 12 bases for each step
                        if check in liste_used.keys():
                            if liste_used[check] != file_iterator:  #saves 12 bases of HS in dictionary if it is not in yet
                                isgud = False
                                break
                        else:
                            liste_used[check] = file_iterator

                    if not i.name == "p300_geordnet.txt":
                        if isgud:
                            # writable
                            datei.write(hs_daten + "\t" + hs + "\n")    #if no alignment with current HS it is saved in text-file
                        else:
                            datei.write("*" + hs_daten + "\t" + hs + " mit " + str(liste_used[check]) + "\n")   #if alignment present it is saved in text-file and marked with '*' for further deletion

def sortieren(pliste_daten):                                            #text-file containing all HS sorts LHS/RHS by GC-content closest to 50%
    liste_dateien = pliste_daten

    for i in liste_dateien:
        liste_sonden = []
        sorting = {}
        with open(i.name, "r") as sonde:
            for j in sonde:
                zeile_test = j.rstrip()
                zeile = zeile_test.split("\t")
                liste_sonden.append(zeile)                               #reading text-file and saving HS in list
            for j in range(0, len(liste_sonden),2):
                gc_diff = abs(50-int(liste_sonden[j][2]))+abs(50-int(liste_sonden[j+1][2]))         #calculation of how close GC-contetn of current LHS/RHS is to 50%
                tm_diff = abs(142-(int(liste_sonden[j][1])+int(liste_sonden[j+1][1])))              #calculation of how close Tm-contetn of current LHS/RHS is to 70°C
                gesamt = gc_diff + tm_diff
                sorting[j] = gc_diff                                                                #HS and GC-difference is saved in dictionary
            sorted_dict = sorted(sorting.items(), key=operator.itemgetter(1))                       #dictioanry is sorted by GC-differnece

            with open(i.name, "w") as zurücksetzen:
                for k in sorted_dict:
                    if k[0]-1 % 2: #ungerade
                        zurücksetzen.write(str(liste_sonden[k[0]][0]) + "\t" + str(liste_sonden[k[0]][3]) + "\n")   #seqeunce plus information about it is written in text-file sorted by GC-contetn closest to 50%
                        zurücksetzen.write(str(liste_sonden[k[0]+1][0]) + "\t" + str(liste_sonden[k[0]+1][3]) + "\n")
                    else: #gerade
                        zurücksetzen.write(str(liste_sonden[k[0]-1][0])+ "\t" + str(liste_sonden[k[0]-1][3]) + "\n")
                        zurücksetzen.write(str(liste_sonden[k[0]][0])+ "\t" + str(liste_sonden[k[0]][3]) + "\n")

def berechnung(pexon, plaenge):                                     #main method that goes through whole sequence of each exon and generates HS by using other methods of this script
    laenge = plaenge * 2                                            #HS generating is done with every defined HS-length
    anzahl_kandidaten = 0
    zexon = identifikation(pexon)[0]
    sonden_datei = identifikation(pexon)[1]

    for i in range(0, len(pexon)):
        if not pexon[i] == "A" and len(pexon[i:len(pexon)]) >= laenge:          #if HS starts with 'A' HS is discarded
            kandidat = pexon[i:i+laenge]
            if i <= 10:
                puffer_lhs = pexon[0:i]                                         #puffer of LHS is determined
            else:
                puffer_lhs = pexon[i-10:i]
            if len(pexon)-(i+laenge) <= 10:
                puffer_rhs = pexon[i+laenge:len(pexon)]                         #puffer of RHS is determined
            else:
                puffer_rhs = pexon[i+laenge:i+laenge+10]

            if primer_bs_neu(kandidat):
                if ligation(kandidat, plaenge) <= 2:
                    lhs = kandidat[0:plaenge]                                    #splitting HS into LHS
                    rhs = kandidat[plaenge:]                                     #and RHS
                    anzahl_kandidaten = anzahl_kandidaten +1
                    gc_gehalt_lhs = gc_gehalt(kandidat, plaenge)[0]
                    gc_gehalt_rhs = gc_gehalt(kandidat, plaenge)[1]
                    datei = open('blast_LHS_RHS.fasta.txt', 'a')

                    if gc_gehalt_rhs >=35:
                        if gc_gehalt_lhs >=35:
                            lhs_anders = gc_rek(puffer_lhs+lhs,"l",plaenge)
                            rhs_anders = gc_rek(rhs+puffer_rhs,"r",plaenge)
                            lhs_temp = lhs_anders[0]+"AAAAA"
                            rhs_temp = "AAAAA"+rhs_anders[0]
                            if lhs_anders[1]>gc_gehalt_lhs and lhs_anders[1] < 55:
                               if not lhs_anders[0][0] == "A":
                                    if primer_bs_neu(lhs_temp):
                                        lhs = lhs_anders[0]
                                        gc_gehalt_lhs = lhs_anders[1]
                            if rhs_anders[1] > gc_gehalt_rhs and rhs_anders[1] < 55:
                                if primer_bs_neu(rhs_temp):
                                    rhs = rhs_anders[0]
                                    gc_gehalt_rhs = rhs_anders[1]

                            datei.write(">" +"LHS(Nr:"+ str(anzahl_kandidaten)+ ", " +zexon + ", Länge: " + str(len(lhs)) + ")\n" + lhs + "\n" + ">" + "RHS(Nr:" + str(anzahl_kandidaten) + ", " +zexon + ", Länge: " + str(len(rhs)) + ")\n" + rhs + "\n" + "\n")                  #writing generated LHS/RHS sequences and information about them in text-file
                            #sonden_datei.write(">LHS(Nr:"+ str(anzahl_kandidaten)+ ", " +zexon + ", Länge: " + str(len(lhs)) + " ,Tm:" + str(Tm_Value(lhs)) + " ,Originallänge: " + str(plaenge) + ")\t" + lhs + "\n" + ">RHS(Nr:" + str(anzahl_kandidaten) + ", " +zexon + ", Länge: " + str(len(rhs)) + " ,Tm: " + str(Tm_Value(rhs)) + " ,Originallänge: " + str(plaenge) + ")\t" +rhs + "\n")
                            #sonden_datei.write(">LHS(Nr:" + str(anzahl_kandidaten) + ", " + zexon + ", Länge: " + str(len(lhs))  + " ,Originallänge: " + str(plaenge) + ")\t" + str(Tm_Value(lhs)) + "\t" + str(gc_gehalt_lhs) +  "\t" + lhs + "\n" + ">RHS(Nr:" + str(anzahl_kandidaten) + ", " + zexon + ", Länge: " + str(len(rhs))  + " ,Originallänge: " + str(plaenge) + ")\t" + str(Tm_Value(rhs)) + "\t" + str(gc_gehalt_rhs) + "\t" + rhs + "\n")
                            with open(sonden_datei.name, "a") as sonden_datei2:                         #creating a text-file for each exon and writing HS generated for them in file
                                sonden_datei2.write(">LHS("+ zexon + ", LHS-Länge: " + str(
                                    len(lhs)) + ", GC-Gehalt: " + str(gc_gehalt_lhs) + "%, " +"Tm: " + str(Tm_Value(lhs)) + "°C)\t" + str(
                                    Tm_Value(lhs)) + "\t" + str(gc_gehalt_lhs) + "\t" + lhs + "\n" + ">RHS("  + zexon + ", RHS-Länge: " + str(
                                    len(rhs)) + ", GC-Gehalt: " + str(gc_gehalt_rhs) + "%, " + "Tm: " + str(Tm_Value(rhs)) +"°C)\t" + str(
                                    Tm_Value(rhs)) + "\t" + str(gc_gehalt_rhs) + "\t" + rhs + "\n")
                                #print(zexon + " " + str(plaenge) + "\nLHS: "  + lhs + " GC-Gehalt: " + str(gc_gehalt_lhs) + "%"  + " Nr.: " + str(anzahl_kandidaten) + "\n" + "RHS: " + rhs +" GC-Gehalt: " + str(gc_gehalt_rhs)  + "% Nr.: " + str(anzahl_kandidaten) + "\n")

def entfernen():                                                          #goes through text-file containing HS of every exon and saves them in new text-file excluding marked HS with too high alignment
    liste_sonden = []
    liste_duplikate = []
    with open("blast_sonden2.txt", "r") as datei:

        behalten = True
        letztes = ""
        index = -1
        for i in datei:
            zeile = i.rstrip().split("\t")
            if zeile[0][0]==">":
                if behalten:
                    letztes = zeile
                    liste_sonden.append(letztes)                            #HS not marked for deletion are saved in list
                    index = index+1
                behalten = True
            if behalten:
                if zeile[0][0]=="*":                                        #checks if HS is marked for deletion
                    if zeile[0][2]=="L":                                    #since LHS and RHS belong together it is checked which one is marked for deletion to delete the corresponding HS
                        behalten = False
                    elif zeile[0][2]=="R":
                        liste_duplikate.append(index)
                        behalten = True
    for i in liste_duplikate:
        liste_sonden.pop(i)
    with open("blast_result_final.txt", "a") as finale:
        vorherige = ""
        for i in liste_sonden:
            exon_name = ""
            for x in range(5,len(i[0])):
                if not i[0][x] == ",":
                    exon_name = exon_name + i[0][x]
                else:
                    break

            if not exon_name == vorherige:
                vorherige = exon_name
                for j in range(0,60):
                    finale.write("%")
                finale.write(" " + vorherige + " ")
                for j in range(0,60):
                    finale.write("%")
                finale.write("\n")
            if i[0][1] == "L":
                finale.write(i[0] + "\t" + i[1] +"\n")                                  #HS with not too high alignment are saved in new text-file and are subdivided by exon to witch they belong to
            elif i[0][1] == "R":
                finale.write(i[0] + "\t" + i[1] + "\n\n")

liste_exone = reading_input()

for i in liste_exone:
    sequenz = i[1]
    for j in laengen:
        berechnung(sequenz, j)                                                  #execute main method berechnung() with each exon sequence

file_names = []
only_names = []
for i in liste_exone:
    name = "sonden_" + i[0] + ".txt"                                            #creating text-file for each exon
    name2 = "sonden_" + i[0]
    with open(name, "r") as name3:
        file_names.append(name3)
        only_names.append(name2)

sortieren(file_names)
local_blast("blast_sonden.txt",file_names)
local_blast("blast_sonden2.txt",[open("p300_geordnet.txt", "r"), open("blast_sonden.txt", "r")])
entfernen()

for i in only_names:                        #removal of text-files only used for calculations
    os.remove(i)
for i in file_names:
    os.remove(i.name)

print(stuffer_lhs.upper())
print(stuffer_rhs.upper())