#! /usr/bin/python3
import sys
import os
import pandas as pd
from datetime import datetime
from datetime import timedelta

sys.path.insert(0,os.getcwd())

from REUS_input import *

"""
Documentation: 
This scripts goal is to benchmanrk the REUS run, 
It serves as 3 main purposes:
-it can count roundtrips (in casethey occur),
-it is able to follow the path of a conformation during a REUS run 
-It is able to benchmark time and exchange number during a run.

Files produced:
RElaps.dat - where it counts the number of roundtrips each US window
does during our REUS run. 
Roundtrips.dat - Outputs total roundtrips of a REUS simulation for each 
US window

RElog.dat -

REexc.dat

REnumberexch.dat - 


List index: 
RE_key - tuple with each US window as a list entry (-10,-9,-8,...)
RE_list - list where the 
RE_exc
RE_n - list where the number of exchanges per segment are counted 

"""

#### Get the Replicate Exchange log file ####

def calc_laps (data):
    REs     = data[0][1:len(data[0])]
    RElaps  = [0     for i in range(len(REs))]
    REctrl  = [0     for i in range(len(REs))]
    REinit  = [False for i in range(len(REs))]
    # Get min and max ##
    f=[float(x) for x in REs]
    REmin= float(REs[f.index(min(f))])
    REmax= float(REs[f.index(max(f))])

    with open('./Benchmark/RElaps.dat', 'w') as f:
        f.write("%s\t" %'0.0')
        for RE in RElaps:
            f.write("%8s" % float(RE))
        f.write('\n')

    for l in data:
        for i in range(len(pHs)):
            if REinit[i] == False:
                if float(l[i+1]) == REmin:
                    REinit[i] = l[0]
                    REctrl[i] = -1
        else:
            if REctrl[i] == -1 and float(l[i+1]) == REmax:
                # if counter started on min (-1) and reach max add halftrip
                RElaps[i] = RElaps[i] + 0.5
                ## and update the controler to the max counter (1)
                REctrl[i] == 1
            if REctrl[i] == 1 and float(l[i+1]) == REmin:
                # if counter started on max (1) and reached min add halftrip
                RElaps[i] = RElaps[i] + 0.5
                ## and update the controler to the min counter (-1)
                REctrl[i] == -1
        
        with open('./Benchmark/RElaps.dat', 'a') as f:
            f.write("%s\t" % l[0])
            for RE in RElaps:
                f.write("%8s" % float(RE))
            f.write('\n')
            
    roundns = [0     for i in range(len(REs))]
    for i in range(len(roundns)):
        roundns[i] = '{:4.2f}'.format(1000*RElaps[i]/(float(data[len(data)-1][0])-float(REinit[i])))
    

    with open('./Benchmark/Roundtrips.dat','w') as r:
        for rl in roundns:
            r.write("%8s" %rl)

    return roundns
############################################

def make_log(RE_l):
    RE_key = tuple(RE_l)
    RE_list = list(RE_key)
    RE_exc = list(RE_key)
    RE_n = [0]
    #RE_win = [0     for i in range(len(RE_l))] #number of exchanges per window
    
    ## Making exchange combination dictionary ##
    RE_win = {}
    for r in range(0, len(RE_key)-1, 1):
        USpair = str(RE_key[r] + "_" + RE_key[r+1] )
        RE_win[USpair] = 0

    for files in ('./Benchmark/RElog.dat', './Benchmark/REexc.dat', './Benchmark/REnumberexch.dat'):
        with open(files, 'w') as log:
            log.write("%s\t" %'0.0')
    
    glob = sorted(os.listdir('GLOB/'))
    for dir in glob:
        if ('GLOBLOG' in dir):
            with open('GLOB/{0}'.format(dir)) as file:
                for line in file:
                    if "REUS step at t" in line:
                        l = line.split()
                        for fi, lis in [['./Benchmark/RElog.dat', RE_list],['./Benchmark/REexc.dat', RE_exc],['./Benchmark/REnumberexch.dat', RE_n]]:
                            with open(fi, 'a') as log :      
                                for RE in lis:
                                    log.write("%8s" % float(RE))
                                log.write('\n')
                                log.write("%s\t" % l[5])
                        RE_n = [0]
                    elif "accepted" in line:
                        l = line.split()
                        ## REgister the exchange on a pair
                        RE_win[l[0]] +=1

                        RE1, RE2 = l[0].split('_')
                        RE1i, RE2i = RE_key.index(RE1), RE_key.index(RE2) # index of 03 and 04 on original
                        RE_list[RE2i], RE_list[RE1i] = RE_list[RE1i], RE_list[RE2i] # value on the list on position 03 (+5) is exchanged 

                        ##### Insert code for the file with the window of each conformation
                        """
                        example original [01, 02, 03, 04, 05] - list [02, 01, 04, 03, 05]
                            exchange of the 04 window with 05
                            In terms of conformation it is an exchange of 05 conf with 03
                            REa1 and REa2 are the number of the conformations exchanged ()
                            REai1 and REai2 are the indexes 
                        """
                        REa1, REa2 = RE_list[RE1i] , RE_list[RE2i] #conformations exchanged (example 05 exchanged with 03)
                        ### Get the column of the exc confo on the exc list
                        REai1, REai2 = RE_key.index(REa1), RE_key.index(REa2)
                        ## change their entries on the RE_exc to the window they are in(got from the original)
                        RE_exc[REai1],RE_exc[REai2] = RE_key[RE1i], RE_key[RE2i]

                        ##### Sum the number of exchanged on each of the windows and overall ####
                        #RE_win[RE1i] += 1
                        #RE_win[RE2i] += 1

                        RE_n[0] += 1


    for fi, lis in [['./Benchmark/RElog.dat', RE_list],['./Benchmark/REexc.dat', RE_exc],['./Benchmark/REnumberexch.dat', RE_n]]:
        with open(fi, 'a') as log :
            for RE in lis:
                log.write("%8s" % float(RE))

    data = []
    ## put the data on a list 
    with open('./Benchmark/REexc.dat', 'r') as log :
        for line in log:
            l = line.split()  
            data.append(l)

    return data, RE_win

def report_benchmark(RE_l, rate):

    ### Hours per block 
    ## Get starting date ##
    totblocktime=(ncycles*(effsteps*dt)/1000) #in ns
    glob = sorted(os.listdir('GLOB/'))
    avgtimes=[]
    ntot = 0

    for dir in glob:
        if ('GLOBLOG' in dir):
            ntot += 1
            with open('GLOB/{0}'.format(dir)) as file:
                timesCpHMD=[]
                timesREUS=[]
                for line in file:                
                    if "CpHMD step" in line:
                        l = line.split()
                        timesCpHMD.append(datetime.strptime(str(l[9]+" "+l[10]), "%Y/%m/%d %H:%M:%S"))
                    if "REUS step at" in line:
                        timesREUS.append(datetime.strptime(str(l[9]+" "+l[10]), "%Y/%m/%d %H:%M:%S"))
                    
            avgtimes.append(timesREUS[-1] - timesCpHMD[0])

    ns = sum(avgtimes, start=timedelta())/(len(avgtimes)*totblocktime)

    ### Make average exchanges per ns ###
    with open('./Benchmark/REnumberexch.dat') as inf:
        texc = 0
        for line in inf:
            line = line.split()
            texc += float(line[1])
    print(rate)
    avgexch=texc/(len(avgtimes)*totblocktime)
    #####################################
    with open('./Benchmark/Benchmark.dat', 'w') as log :
        log.write("%s\n\n" %'Reporting REUS benchmark data:')
        log.write("%s\t%s\n\n" %('Nanoseconds per day:', ns))
        log.write("%s\n\n" %'Reporting REUS rates')
        log.write("%s\t%5.1f\n\n" %('Exchanges per ns:', avgexch))
        log.write("%s\n" %'Reaction Coordinate exchange per ns :')
        for i in rate.keys():
            log.write("%4s\t%d\t%.2f\n" %(i,(rate[i]),(rate[i]/ntot)))
        #for i in range(len(rate)):
        #    log.write("%4s\t%d\t%.2f\n" %(RE_l[i],(rate[i]),(rate[i]/ntot)))
        



if __name__ == "__main__":
    ### make benchmark directory ###
    os.system('mkdir -p ./Benchmark')

    US_list=list(USs.keys())
    data,rate = make_log(US_list)

    roun = calc_laps(data)

    report_benchmark(US_list, rate)

    #print(roun)

