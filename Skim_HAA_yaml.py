# !/usr/bin/env python
#########################
#Author: Sam Higginbotham
'''

* File Name : Skim_HAA_yaml.py

* Purpose : For skimming and combining ntuples. CAREFUL! May cause cuts not expected

* Creation Date : june-30-2020

* Last Modified :

'''
#########################
import sys
import os
import ROOT
import numpy as np
#import argsparse
from array import array
import itertools
import operator
import csv
import datetime
import time
import threading 
import yaml
import glob
#user definitions 
from utils.Parametrization import *

#Function that combines neighboring bins into a list 
def returnBinList(flatbins):
    comboBins=[]
    for sl in range(0,len(flatbins)):
        comboBins.append([])
        for bin in range(0,len(flatbins[sl])-1):
            comboBins[sl].append([flatbins[sl][bin],flatbins[sl][bin+1]])
    return comboBins

def returnPerms(setlist):

    combinations=list(itertools.product(*setlist))

    return combinations


if __name__ == "__main__":
    
    begin_time = datetime.datetime.now()



    import io
    import os

    import argparse

    parser = argparse.ArgumentParser(description="This file generates root files containing Histograms ... files in utils contain selections and settings")
    #parser.add_arguement('--CategoryFiles',nargs="+",help="Select the files containing the categories for the datacards")
    parser.add_argument("-o",  "--outname", default="",  help="postfix string")
    parser.add_argument("-fi",  "--ffin", default="",  help="fake factor files")
    parser.add_argument("-fo",  "--ffout", default="",  help="fake factor files to output")
    parser.add_argument("-c",  "--categories", default="categories.yaml",  help="categories yaml file")
    parser.add_argument("-csv",  "--csvfile", default="MCsamples_2016_v6_yaml.csv",  help="categories yaml file")
    parser.add_argument("-p",  "--processes", default="processes_special.yaml",  help="processes yaml file")
    parser.add_argument("-dm",  "--datameasure", default=False,action='store_true',  help="Use DataDriven Method measure part")
    parser.add_argument("-dd",  "--datadriven", default=False,action='store_true',  help="Use DataDriven Method")
    parser.add_argument("-dmZH",  "--datameasureZH", default=False,action='store_true',  help="Use DataDriven Method measure part")
    parser.add_argument("-ddZH",  "--datadrivenZH", default=False,action='store_true',  help="Use DataDriven Method")
    parser.add_argument("-ff",  "--makeFakeHistos", default=False,action='store_true',  help="Just make fake rate histos")
    parser.add_argument("-v",  "--verbose", default=False,action='store_true',  help="print per event")
    parser.add_argument("-t",  "--test", default=False,action='store_true',  help="only do 1 event to test code")
    parser.add_argument("-m",  "--merge", default=False,action='store_true',  help="Merge the files into processes?")
    parser.add_argument("-mt",  "--mt", default=False,action='store_true',  help="Use Multithreading")
    parser.add_argument("-pt",  "--maxprint", default=False,action='store_true',  help="Print Info on cats and processes")
    args = parser.parse_args()


    #gather functions for computing variables in the event loop
    from utils.functions import functs

    #Structure for plotting variables 
    from utils.Parametrization import Category
    from utils.Parametrization import Process

    #Gather the analysis datasets and info 
    sampleDict = {}
 

    for line in open(args.csvfile,'r').readlines() :
            #[nickname]        = [category,xsec,numberOfEvents,finishedEvents,idk?,DASDataset]
            if len(line.split(',')[2].split("*"))>1:
                tempval=1.0
                for val in line.split(',')[2].split("*"):
                    tempval = tempval * float(val)
                row = line.split(',')
                sampleDict[row[0]] = [row[1],tempval,row[3],row[4],row[5],row[6]]
            else:
                row = line.split(',')
                sampleDict[row[0]] = [row[1],float(row[2]),row[3],row[4],row[5],row[6]]

    #importing analysis categories and conventions
    with io.open(args.categories,'r') as catsyam:
        categories = yaml.load(catsyam)

    #loading fake factor and data driven methods
    with io.open(args.processes,'r') as prosyam:
        processes_special = yaml.load(prosyam)

    allcats={}

    for category in categories:
        #print category
        #print categories[category]['name']
        tempcat = Category() 
        tempcat.name=categories[category]['name']
        tempcat.cuts=categories[category]['cuts']
        tempcat.newvariables=categories[category]['newvariables']
        tempcat.vars=categories[category]['vars']
        allcats[tempcat.name]=tempcat

    
    #loading standard processes
    HAA_processes={}
    for sample in sampleDict.keys():
        #print processes[process]
        temppro = Process() 
        temppro.nickname=sample
        temppro.file=sample+"_2016.root"
        temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3]}
        temppro.cuts={sampleDict[sample][0]:""}
        if "ggTo2mu2tau" in sample:
            temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)}
        if "W" in sample and "Jets" in sample: 
            temppro.file=sample+"_2016.root"
            temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3]}
            temppro.cuts={"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]}
        if "TT" in sample and not "TTTT" in sample and not "TTHH" in sample: 
            temppro.file=sample+"_2016.root"
            temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3]}
            temppro.cuts={"TT":"","TTT":[["gen_match_4","==",5]],"TTL":[["gen_match_4",">=",5]],"TTJ":[["gen_match_4",">",5]]}
        HAA_processes[temppro.nickname]=temppro

    #loading special processes ... fake factor and data
    for process in processes_special:
        temppro = Process() 
        temppro.nickname=processes_special[process]['nickname']
        temppro.cuts=processes_special[process]['cuts']
        temppro.weights=processes_special[process]['weights']
        temppro.file=processes_special[process]['file']
        HAA_processes[temppro.nickname]=temppro
        
    if args.maxprint:
        print "categories "
        for cat in allcats:
            print cat.name
            print cat.cuts
        

        print "processes "
        for pro in HAA_processes.keys():
            print HAA_processes[pro].nickname
            print HAA_processes[pro].cuts
            print HAA_processes[pro].weights
            print HAA_processes[pro].file
    
             
    
    #Gather the Analysis Files
    os.system("rm -r /afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/skimtest/")
    os.system("cp -r /afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/ /afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/skimtest/")
    dir = "/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/skimtest/"
    #filelist = HAA
    filelist = {}
    processes = []

    for pro in HAA_processes.values():
        for subpro in pro.cuts.keys():
            processes.append(subpro)
    processes=list(dict.fromkeys(processes))
    print "Will plot the following physics processes ",processes
    print "added SF weights"

    #for testing only... one cat a15
    for proObj in HAA_processes.keys():
       filelist[proObj]=HAA_processes[proObj].file

    print filelist


    treename = "Events"

    print "initializing histograms"
    
    try:
        os.mkdir(dir+"/skim_"+str(args.outname))
    except:
        print "directory exists"

    filedict = {}
    skimDict = {}
    outFiles = {}
    treeList = {}
    outTrees = {}
    
    numvar=0

    #combining cuts
    print allcats["mmmt_inclusive"].cuts["preselection"]
    cutstring = ""
    for cut in allcats["mmmt_inclusive"].cuts["preselection"]:
        #print cut,cut[0],cut[0][0]
        if cut[0][0]=="EQT" and cut[2]=="mult":
            cutstring = cutstring + cut[1][0]+"*"+cut[1][1]+cut[3]+str(cut[4])
        elif cut[1]=="absl":
            cutstring = cutstring + "-"+str(cut[2])+"<"+cut[0]+"&&"+cut[0]+"<"+str(cut[2])
        else:
            cutstring = cutstring + cut[0]+cut[1]+str(cut[2])
            
        if cut!=allcats["mmmt_inclusive"].cuts["preselection"][-1]:
            cutstring = cutstring + "&&"

    print "The preselections ",cutstring
    #making skims
    for nickname in filelist.keys():

        print "working on file ",filelist[nickname]
        skimDict[nickname] = [ROOT.TFile.Open(dir+filelist[nickname],"update")]
        skimDict[nickname].append(skimDict[nickname][0].Get(treename))
        
        #newtree=tree.CloneTree()
        print "tree entries ",skimDict[nickname][1].GetEntries()

        
        
        skimDict[nickname].append(skimDict[nickname][1].CopyTree(cutstring))
        print "new tree entries ",skimDict[nickname][2].GetEntries()
        skimDict[nickname][0].cd()
        skimDict[nickname][2].Write(skimDict[nickname][1].GetName(),ROOT.TObject.kOverwrite)
        

        #processes will be a dictionary of key: process name value: cuts
        #for process in processes:
        #    treeList[process]=ROOT.TList()

        #for process in processes:
        skimDict[nickname][0].Close()
        #del skimDict[nickname][1]
        #del skimDict[nickname][0]
        

    #Now merging TTree/Files         
    if args.merge:
        print "merging by TList"
        for process in processes:
            mergefiles="" 
            for nickname in filelist.keys():
                for pro in HAA_processes[nickname].cuts.keys():
                    if process==pro:
                        mergefiles = mergefiles +" "+ filelist[nickname] 
            print "hadd "+dir+"skim_"+str(args.outname)+"/"+process+".root "+mergefiles
            #os.system("hadd "+dir+"/skim_"+str(args.outname)+process+".root "+mergefiles)

        
        

    print "computation time"        
    print(datetime.datetime.now() - begin_time)

