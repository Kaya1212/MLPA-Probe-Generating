from math import  floor
import operator

stuffer_lhs = "CGCTACTACTATTAGTAGAATTGATGCCACCTTTTCAGCTCGCG"                        #sequence of stuffer sequence of LHS
stuffer_rhs = "ACCATTTGCGAAATGTATCTAATGGTCAAACTAAATCTAC"                            #sequence of stuffer sequence of RHS
lhs_primer = "GGGTTCCCTAAGGGTTGGA"                                                  #primer sequence LHS
rhs_primer = "TCTAGATTGGATCTTGCTGGCAC"                                              #primer sequence RHS
p300_laengen = [64,82,88,92,108,130,148,173,180,184,191,196,208,214,220,226,232,239,246,252,258,265,274]        #probe lengths of P300-Probemix
kandidaten_liste = []

kandidaten = []
zurücksetzten = open("MLPA_SONDEN.txt", "w")
schreiben = open("MLPA_SONDEN.txt", "a")

def datei_lesen():                                                                      #reading input of sequences
    with open("auswahl.txt", "r") as auswahl:
        for i in auswahl:
            zeile = i.rstrip().split("\t")
            kandidaten_liste.append(zeile)                                               #saving sequences in list

def listen_erstellen():                                                                   #attaches primer-binding seqeunce to LHS/RHS and saves them in a list sorted by LPO+RPO length
    for i in range(0, len(kandidaten_liste), 2):
        lpo = lhs_primer + kandidaten_liste[i][1]                                         #attaching left primer-binding sequence to LHS
        rpo = kandidaten_liste[i+1][1] + rhs_primer                                       #attaching right primer-binding sequence to RHS
        exon_name = ""
        for x in range(5,len(kandidaten_liste[i][0])):
            if not kandidaten_liste[i][0][x] == ",":
                exon_name = exon_name + kandidaten_liste[i][0][x]                         #saves exon-name of LHS/RHS
            else:
                break
        info = exon_name
        lpo_länge = len(lpo)
        rpo_länge = len(rpo)
        temp = []
        temp.append(info)                                                                 #saves exon-name in list
        temp.append(lpo_länge+rpo_länge)                                                  #saves LPO+RPO length in list
        temp.append(lpo)                                                                  #saves LPO sequence in list
        temp.append(rpo)                                                                  #saves RPO sequence in list
        kandidaten.append(temp)
        #print(info +"\t" + str(lpo_länge) + "\t" + lpo + "\t" + rpo + "\t" + str(rpo_länge))
    print(kandidaten_liste)
    kandidaten.sort(key=operator.itemgetter(1))                                           #list gets sorted by LPO+RPO length
    print(kandidaten)

def laengen_berechnen():                                                                #calculates free lengths for probes to guarantee maximum distance between probes
    verfügbare_langen = []
    abstand= 0
    for k in range(4,13):

        aktuelle_laengen=[]
        for i in range(0, len(p300_laengen)):
            if not i == len(p300_laengen)-1:
                differenz = p300_laengen[i+1]-p300_laengen[i]          #calcualtes how much bases are between neighbouring probes of P300-probemix
                if differenz>=k*2:                                     #checks if difference is great enough for usage
                    space = floor((differenz-k)/k)
                    zwischenlänge = p300_laengen[i]
                    for j in range(0,space):
                        zwischenlänge = zwischenlänge +k
                        aktuelle_laengen.append(zwischenlänge)          #calculates all possible free lenghts for probes
                    platz = 0
                    if len(aktuelle_laengen) >= len(kandidaten):        #checks if there are enough possible lengths for probes
                        for l in aktuelle_laengen:
                            if l>= kandidaten[0][1] and l<=140:         #checks if probe is longer than possible lenght
                                platz=platz+1
                    if platz >= len(kandidaten):
                        verfügbare_langen = aktuelle_laengen            #lenghts that can be used for probes
                        abstand = k                                     #maximum distance between probes
    print("Abstand zwischen Sonden: " + str(abstand))
    print("verfügbare Längen: " +str(verfügbare_langen))
    return verfügbare_langen



def formatierung(kandidaten_info,stuffer_lpo,stuffer_rpo):              #edits the final output text-file

    with open("MLPA_SONDEN.txt", "a") as sonden_file:
        for k in range(0, 60):
            sonden_file.write("%")
        sonden_file.write(" " + kandidaten_info[0] + " ")
        for k in range(0, 60):
            sonden_file.write("%")
        sonden_file.write("\n")

        sonden_file.write(">LPO:\n"+lhs_primer+stuffer_lpo+kandidaten_info[2][len(lhs_primer):]+"\n" )
        for i in range(0,len(lhs_primer)):
            sonden_file.write("-")
        for i in range(0,len(stuffer_lpo)):
            sonden_file.write("+")
        for i in range(0, len(kandidaten_info[2])-len(lhs_primer)):
            sonden_file.write("*")
        sonden_file.write("\n\n")

        sonden_file.write(">RPO:\n" + kandidaten_info[3][:len(kandidaten_info[3])-len(rhs_primer)]+stuffer_rpo+rhs_primer+"\n")
        for i in range(0,len(kandidaten_info[3])-len(rhs_primer)):
            sonden_file.write("*")
        for i in range(0,len(stuffer_rpo)):
            sonden_file.write("+")
        for i in range(0,len(rhs_primer)):
            sonden_file.write("-")
        sonden_file.write("\n"+"Gesamtlänge: " + str(kandidaten_info[1]+len(stuffer_lpo)+len(stuffer_rpo)))
        sonden_file.write("\n\n")

def stuffer_wählen():                                               #calculates ideal stuffer sequence lenght for each probe half
    verfügbare_langen = laengen_berechnen()
    for i in kandidaten:
        kandidat_laenge = i[1]
        for j in verfügbare_langen:
            if kandidat_laenge==j:                                  #if probe is as long as the free lenght calculated in laengen_berechnen(), the probe becomes an input for formatierung()
                verfügbare_langen.remove(j)
                formatierung(i,"")

            elif j>kandidat_laenge:
                stuffer_laenge = j-kandidat_laenge
                verfügbare_langen.remove(j)
                lpo_laenge = len(i[2])
                rpo_leange = len(i[3])
                differenz = lpo_laenge-rpo_leange
                lpo_stuffer = ""
                rpo_stuffer = ""
                rpo_temp_benutzen = False
                lpo_temp_benutzen = False

                if differenz<0:#RPO größer                          #if RPO is longer, LPO will have more stuffer sequence to be as long as RPO
                    lpo_stuffer = stuffer_lhs[0:abs(differenz)]
                    stuffer_laenge = stuffer_laenge-abs(differenz)
                    lpo_temp_benutzen = True
                    lpo_stuffer_temp = lpo_stuffer
                    if stuffer_laenge<1:
                        if abs(differenz) > abs(stuffer_laenge):
                            lpo_temp_benutzen = True
                            lpo_stuffer = stuffer_lhs[0:abs(stuffer_laenge)]
                            rpo_stuffer = ""
                            lpo_stuffer_temp = lpo_stuffer
                            formatierung(i, lpo_stuffer, rpo_stuffer)
                        else:
                            lpo_temp_benutzen = True
                            lpo_stuffer = stuffer_lhs[0:abs(differenz)-1]
                            rpo_stuffer = ""
                            lpo_stuffer_temp = lpo_stuffer
                            formatierung(i, lpo_stuffer, rpo_stuffer)

                if differenz>0:#LPO größer                          #if LPO is longer, RPO will have more stuffer sequence to be as long as LPO
                    rpo_stuffer =stuffer_rhs[len(stuffer_rhs)-abs(differenz):]
                    stuffer_laenge = stuffer_laenge-abs(differenz)
                    rpo_temp_benutzen = True
                    rpo_stuffer_temp = rpo_stuffer
                    if stuffer_laenge<1:
                        if abs(differenz) > abs(stuffer_laenge):
                            lpo_stuffer = ""
                            rpo_stuffer = stuffer_rhs[len(stuffer_rhs) - abs(stuffer_laenge):]
                            rpo_stuffer_temp = rpo_stuffer
                            formatierung(i, lpo_stuffer, rpo_stuffer)
                        else:
                            lpo_stuffer = ""
                            rpo_stuffer = stuffer_rhs[len(stuffer_rhs)-abs(differenz)+1:]
                            rpo_stuffer_temp = rpo_stuffer
                            formatierung(i, lpo_stuffer, rpo_stuffer)

                if differenz==0:
                    hälfte = floor(stuffer_laenge / 2)
                    überstand = stuffer_laenge % 2
                    lpo_stuffer = stuffer_lhs[0:hälfte+überstand]
                    rpo_stuffer = stuffer_rhs[len(stuffer_rhs)-hälfte:]
                    formatierung(i, lpo_stuffer, rpo_stuffer)
                    break

                if stuffer_laenge==1:
                    lpo_stuffer =lpo_stuffer+stuffer_lhs[len(lpo_stuffer)]
                    formatierung(i, lpo_stuffer, rpo_stuffer)
                    break

                if stuffer_laenge>1:
                    hälfte = floor(stuffer_laenge/2)
                    überstand = stuffer_laenge % 2

                    if lpo_temp_benutzen :
                        lpo_stuffer = lpo_stuffer_temp+stuffer_lhs[len(lpo_stuffer_temp):len(lpo_stuffer_temp)+hälfte+überstand]
                        rpo_stuffer = stuffer_rhs[
                                      len(stuffer_rhs) - hälfte:len(stuffer_rhs) - len(rpo_stuffer)] + rpo_stuffer
                        formatierung(i,lpo_stuffer,rpo_stuffer)

                    elif rpo_temp_benutzen:
                        rpo_stuffer = stuffer_rhs[len(stuffer_rhs) - hälfte:len(stuffer_rhs) - len(rpo_stuffer_temp)] + rpo_stuffer_temp
                        lpo_stuffer = lpo_stuffer + stuffer_lhs[len(lpo_stuffer):len(lpo_stuffer) + hälfte + überstand + abs(differenz)]
                        formatierung(i, lpo_stuffer, rpo_stuffer)

                    else:
                        lpo_stuffer = lpo_stuffer + stuffer_lhs[len(lpo_stuffer):len(lpo_stuffer) + hälfte + überstand]
                        rpo_stuffer = stuffer_rhs[len(stuffer_rhs) - hälfte:len(stuffer_rhs) - len(rpo_stuffer)] + rpo_stuffer
                        formatierung(i, lpo_stuffer, rpo_stuffer)


                break


datei_lesen()
listen_erstellen()
stuffer_wählen()
