# !/usr/bin/env python
#########################
#Author: Sam Higginbotham
'''

* File Name : MakeDistributions.py

* Purpose : For visualizing datasets. Will create ultra skimmed root files or histograms for visualization and for limit setting later

* Usage:
python MakeDistributions.py -c cat_mmtt_2016.yaml  -csv MCsamples_2016_v6_yaml.csv  -i /path/to/input/root/files/ -p processes_special_mmtt.yaml -o main_outputdirectory -fo fake_factor_measure_outputdirectory -ch mmtt -s

* Creation Date : june-20-2020

* Last Modified : feb-16-2021

'''
#########################
import sys
import os
import ROOT
import uproot
import pandas
import root_numpy
import numpy as np
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



#Setting up operators for cut string iterator
ops = { "==": operator.eq, "!=": operator.eq, ">": operator.gt, "<": operator.lt, ">=": operator.ge, "<=": operator.le, "band": operator.and_,"bor":operator.or_}

class myThread (threading.Thread):
   def __init__(self, threadID, name, counter):
      threading.Thread.__init__(self)
      self.threadID = threadID
      self.name = name
      self.counter = counter
   def run(self):
      print("Starting " + self.name)
      # Get lock to synchronize threads
      threadLock.acquire()
      print_time(self.name, self.counter, 3)
      # Free lock to release next thread
      threadLock.release()

def print_time(threadName, delay, counter):
   while counter:
      time.sleep(delay)
      print("%s: %s" % (threadName, time.ctime(time.time())))
      counter -= 1

def subArrayLogic(evt,array):
    boo=ops[array[1]](returnValue(evt,array[0]),array[2])
    return boo


def returnArray(masterArray,variable):
    #if variable in ["njets","jpt_1","jeta_1","jpt_2","jeta_2","bpt_1","bpt_2","nbtag","beta_1","beta_2"]:
    #    val = masterArray[variable][:,0]
    #else:
    #print "working on branch ",variable
    #if len(masterArray[variable].shape) == 1:   # flat array important for jagged arrays of input data
    #    val = masterArray[variable]
    if "[" in variable:
        basevar = variable.split("[")[0]
        index = int(variable.split("[")[1].split("]")[0])
        val = masterArray[basevar][:,index]
    else:
        val = masterArray[variable]
    return val

#create mask for all the events that do or don't pass cuts
def cutOnArray(masterArray,cuts):

    mask=np.full(len(masterArray["evt"]),True)

    for cut in cuts:
        #Every USE CASE MUST MULTIPLY BY MASK!
        #print cut
        if type(cut[0][0])==type([0]):
            mask = cutOnArray(masterArray,cut)

        if cut[0][0]=="EQT":
            for i,var in enumerate(cut[1]):
                if cut[2]=="mult":
                    if(i==0):
                        tempmask=np.full(len(masterArray["evt"]),1.0)
                    tempmask = tempmask * returnArray(masterArray,var)

                if cut[2]=="div":
                    if(i==0):
                        tempmask=np.full(len(masterArray["evt"]),1.0)
                    tempmask = tempmask / returnArray(masterArray,var)
                if cut[2]=="add":
                    if(i==0):
                        tempmask=np.full(len(masterArray["evt"]),0.0)
                    tempmask = tempmask + returnArray(masterArray,var)
                if cut[2]=="sub":
                    if(i==0):
                        tempmask=np.full(len(masterArray["evt"]),0.0)
                    tempmask = tempmask - returnArray(masterArray,var)

            tempmask = ops[cut[3]](tempmask,cut[4])

            mask *= tempmask.astype(bool)
            continue

        if cut[0][0]=="OR":
            tempmasks=[]
            tempmask=np.full(len(masterArray["evt"]),False)
            for oe in range(1,len(cut)):
                oneor = ops[cut[oe][1]](returnArray(masterArray,cut[oe][0]),cut[oe][2])
                tempmask = ops["bor"](tempmask,oneor).astype(bool)
            mask *= tempmask.astype(bool)
            continue

        else:
            statement = cut[1]
            if statement=="absg": # lessthan OR greater than
                tempmask = ops["<"](returnArray(masterArray,cut[0]),(-1 * cut[2]))
                tempmask2= ops[">"](returnArray(masterArray,cut[0]),cut[2])
                tempmask = ops["bor"](tempmask,tempmask2)
            if statement=="absl": # greater than AND less tahn
                tempmask = ops[">"](returnArray(masterArray,cut[0]),(-1 * cut[2]))
                tempmask *= ops["<"](returnArray(masterArray,cut[0]),cut[2])
            else:
                tempmask = ops[statement](returnArray(masterArray,cut[0]),cut[2])
            mask *= tempmask.astype(bool)


    return mask


#for fake factor return the bin value from event content
def ptFun(pt,numpyArr):
    newArr = np.full(len(numpyArr),1.0)
    newArr = np.vectorize(pt.FindBin)(numpyArr)
    newArr = np.vectorize(pt.GetBinContent)(newArr)

    return newArr

def getEventWeightDicitonary():
    #import copyreg, copy, pickle # for picking the bound methods due to the multiplrocess tool

    from TauPOG.TauIDSFs.TauIDSFTool import TauIDSFTool
    from TauPOG.TauIDSFs.TauIDSFTool import TauESTool
    from TauPOG.TauIDSFs.TauIDSFTool import TauFESTool
    tauIDSF = TauIDSFTool('2016Legacy','DeepTau2017v2p1VSjet','Medium').getSFvsPT

    #def _pickle_method(method):
    #    func_name = method.im_func.__name__
    #    obj = method.im_self
    #    cls = method.im_class
    #return _unpickle_method, (func_name, obj, cls)
#
    #def _unpickle_method(func_name, obj, cls):
    #    for cls in cls.mro():
    #    try:
    #    func = cls.__dict__[func_name]
    #    except KeyError:
    #    pass
    #    else:
    #    break
    #return func.__get__(obj, cls)
    #copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)


    EventWeights={
        #"name":[[if statements],[weight to apply]]
        "3_mt_lt0p4":   [[["cat","==",6],["decayMode_3","==",0],["eta_3","<",0.4]],[0.80]],
        "3_mj_lt0p4":   [[["cat","==",6],["decayMode_3","==",0],["eta_3","<",0.4]],[1.21]],
        "3_mt_0p4to0p8":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",0.4],["eta_3","<",0.8]],[0.81]],
        "3_mj_0p4to0p8":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",0.4],["eta_3","<",0.8]],[1.11]],
        "3_mt_0p8to1p2":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",0.8],["eta_3","<",1.2]],[0.79]],
        "3_mj_0p8to1p2":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",0.8],["eta_3","<",1.2]],[1.2]],
        "3_mt_1p2to1p7":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",1.2],["eta_3","<",1.7]],[0.68]],
        "3_mj_1p2to1p7":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",1.2],["eta_3","<",1.7]],[1.16]],
        "3_mt_1p7to2p3":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",1.7],["eta_3","<",2.3]],[0.68]],
        "3_mj_1p7to2p3":[[["cat","==",6],["decayMode_3","==",0],["eta_3",">",1.7],["eta_3","<",2.3]],[2.25]],

        "4_mt_lt0p4":   [[["cat","==",6],["decayMode_4","==",0],["eta_4","<",0.4]],[0.80]],
        "4_mj_lt0p4":   [[["cat","==",6],["decayMode_4","==",0],["eta_4","<",0.4]],[1.21]],
        "4_mt_0p4to0p8":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",0.4],["eta_4","<",0.8]],[0.81]],
        "4_mj_0p4to0p8":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",0.4],["eta_4","<",0.8]],[1.11]],
        "4_mt_0p8to1p2":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",0.8],["eta_4","<",1.2]],[0.79]],
        "4_mj_0p8to1p2":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",0.8],["eta_4","<",1.2]],[1.2]],
        "4_mt_1p2to1p7":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",1.2],["eta_4","<",1.7]],[0.68]],
        "4_mj_1p2to1p7":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",1.2],["eta_4","<",1.7]],[1.16]],
        "4_mt_1p7to2.3":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",1.7],["eta_4","<",2.3]],[0.68]],
        "4_mj_1p7to2.3":[[["cat","==",6],["decayMode_4","==",0],["eta_4",">",1.7],["eta_4","<",2.3]],[2.25]],

        "3_et_lt1p479_DM0":[[["cat","==",5],["decayMode_3","==",0],["eta_3","<",1.479]],[0.80]],
        "3_ej_lt1p479_DM0":[[["cat","==",5],["decayMode_3","==",0],["eta_3","<",1.479]],[1.18]],
        "3_et_gt1p479_DM0":[[["cat","==",5],["decayMode_3","==",0],["eta_3",">",1.479]],[0.72]],
        "3_ej_gt1p479_DM0":[[["cat","==",5],["decayMode_3","==",0],["eta_3",">",1.479]],[0.93]],
        "3_et_lt1p479_DM1":[[["cat","==",5],["decayMode_3","==",1],["eta_3","<",1.479]],[1.14]],
        "3_ej_lt1p479_DM1":[[["cat","==",5],["decayMode_3","==",1],["eta_3","<",1.479]],[1.18]],
        "3_et_gt1p479_DM1":[[["cat","==",5],["decayMode_3","==",1],["eta_3",">",1.479]],[0.64]],
        "3_ej_gt1p479_DM1":[[["cat","==",5],["decayMode_3","==",1],["eta_3",">",1.479]],[1.07]],

        "4_et_lt1p479_DM0":[[["cat","==",5],["decayMode_4","==",0],["eta_4","<",1.479]],[0.80]],
        "4_ej_lt1p479_DM0":[[["cat","==",5],["decayMode_4","==",0],["eta_4","<",1.479]],[1.18]],
        "4_et_gt1p479_DM0":[[["cat","==",5],["decayMode_4","==",0],["eta_4",">",1.479]],[0.72]],
        "4_ej_gt1p479_DM0":[[["cat","==",5],["decayMode_4","==",0],["eta_4",">",1.479]],[0.93]],
        "4_et_lt1p479_DM1":[[["cat","==",5],["decayMode_4","==",1],["eta_4","<",1.479]],[1.14]],
        "4_ej_lt1p479_DM1":[[["cat","==",5],["decayMode_4","==",1],["eta_4","<",1.479]],[1.18]],
        "4_et_gt1p479_DM1":[[["cat","==",5],["decayMode_4","==",1],["eta_4",">",1.479]],[0.64]],
        "4_ej_gt1p479_DM1":[[["cat","==",5],["decayMode_4","==",1],["eta_4",">",1.479]],[1.07]],

        "8_3_et_lt1p479_DM0":[[["cat","==",8],["decayMode_3","==",0],["eta_3","<",1.479]],[0.80]],
        "8_3_ej_lt1p479_DM0":[[["cat","==",8],["decayMode_3","==",0],["eta_3","<",1.479]],[1.18]],
        "8_3_et_gt1p479_DM0":[[["cat","==",8],["decayMode_3","==",0],["eta_3",">",1.479]],[0.72]],
        "8_3_ej_gt1p479_DM0":[[["cat","==",8],["decayMode_3","==",0],["eta_3",">",1.479]],[0.93]],
        "8_3_et_lt1p479_DM1":[[["cat","==",8],["decayMode_3","==",1],["eta_3","<",1.479]],[1.14]],
        "8_3_ej_lt1p479_DM1":[[["cat","==",8],["decayMode_3","==",1],["eta_3","<",1.479]],[1.18]],
        "8_3_et_gt1p479_DM1":[[["cat","==",8],["decayMode_3","==",1],["eta_3",">",1.479]],[0.64]],
        "8_3_ej_gt1p479_DM1":[[["cat","==",8],["decayMode_3","==",1],["eta_3",">",1.479]],[1.07]],

        "8_4_mt_lt0p4":   [[["cat","==",8],["decayMode_4","==",0],["eta_4","<",0.4]],[0.80]],
        "8_4_mj_lt0p4":   [[["cat","==",8],["decayMode_4","==",0],["eta_4","<",0.4]],[1.21]],
        "8_4_mt_0p4to0p8":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",0.4],["eta_4","<",0.8]],[0.81]],
        "8_4_mj_0p4to0p8":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",0.4],["eta_4","<",0.8]],[1.11]],
        "8_4_mt_0p8to1p2":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",0.8],["eta_4","<",1.2]],[0.79]],
        "8_4_mj_0p8to1p2":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",0.8],["eta_4","<",1.2]],[1.2]],
        "8_4_mt_1p2to1p7":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",1.2],["eta_4","<",1.7]],[0.68]],
        "8_4_mj_1p2to1p7":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",1.2],["eta_4","<",1.7]],[1.16]],
        "8_4_mt_1p7to2.3":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",1.7],["eta_4","<",2.3]],[0.68]],
        "8_4_mj_1p7to2.3":[[["cat","==",8],["decayMode_4","==",0],["eta_4",">",1.7],["eta_4","<",2.3]],[2.25]],

        #bound method ... last input list is func parameters! so this uses the tauIDSF function
        "3_tauSF_8":[[["cat","==",8],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
        "4_tauSF_8":[[["cat","==",8],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]],
        "3_tauSF_7":[[["cat","==",7],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
        "4_tauSF_7":[[["cat","==",7],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]],
        "3_tauSF_6":[[["cat","==",6],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
        "4_tauSF_6":[[["cat","==",6],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]],
        "3_tauSF_5":[[["cat","==",5],["gen_match_3","==",5]],[np.vectorize(tauIDSF),["pt_3","gen_match_3"]]],
        "4_tauSF_5":[[["cat","==",5],["gen_match_4","==",5]],[np.vectorize(tauIDSF),["pt_4","gen_match_4"]]]

        #recoil corrections? not affecting current fit variable
        #"recoilcorr2":[[["cat","==",6],["njets","==",2]],[recoilCorrector,["met_x","met_y","gen_match_3"]]],
    }

    return EventWeights

def initialize():
    import os
    import datetime
    from copy import copy
    import yaml
    begin_time = datetime.datetime.now()

    #Structure for mapping between root files and processes
    from utils.ProcessCuts import HAA_ProCuts

    #Structure for mapping between processes and weights
    from utils.Weights import CommonWeights
    from utils.Weights import HAAWeights


    #sys.path.append('../SFs/')
    #import utils.SFs.ScaleFactor as SF

    import io
    import os
    from shutil import copyfile


    #gather functions for computing variables in the event loop
    from utils.functions import functs
    from ROOT import gInterpreter

    #Structure for plotting variables
    from utils.Parametrization import Category#, DylanCategory
    from utils.Parametrization import Process

    #Gather the analysis datasets and info
    sampleDict = {}
    #csvfile = "MCsamples_2016_v7_ZZ.csv"
    csvfile = "MCsamples_2016_v7.csv"
    categories = "cat_mmtt_2016.yaml"
    processes = "processes_special_mmtt.yaml"

    dir = "/eos/home-s/shigginb/HAA_ntuples/2016_v7/"


    weightHistoDict = {}
    for nickname in filelist.keys():
        try:
            frooin = ROOT.TFile.Open(dir+filelist[nickname],"read")
            weightHistoDict[nickname]=frooin.Get("hWeights")
        except:
            print nickname," hWeights prob doesn't exist?"
    for line in open(csvfile,'r').readlines() :
            # structure for csv to sampleDict conversion
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
    with io.open(categories,'r') as catsyam:
        categories = yaml.load(catsyam)

    #loading fake factor and data driven methods
    with io.open(processes,'r') as prosyam:
        processes_special = yaml.load(prosyam)

    allcats={}

    for category in categories:
        #print(category)
        #print(categories[category]['name'])
        if categories[category]['name']=="mmtt_inclusive":
            tempcat = Category()
            tempcat.name=categories[category]['name']
            tempcat.cuts=categories[category]['cuts']
            tempcat.newvariables=categories[category]['newvariables']
            tempcat.vars=categories[category]['vars']
            tempcat.systematics=categories[category]['systematics']
            allcats[tempcat.name]=tempcat

    #allcats = {name: DylanCategory(cat) for name, cat in categories.items()}


    #for cat in allcats.keys():
        #print(allcats[cat].name)
        #print(allcats[cat].cuts)
        #print(allcats[cat].systematics)

    #loading standard processes
    HAA_processes={}

    #MC for background

    #HAA_process = {name: Process(dir, name, info) for name, info in sampleDict.items()}
    weightHistoDict = {}


    for sample in sampleDict.keys():
        temppro = Process()
        temppro.nickname=sample
        temppro.file=dir+sample+"_2016.root"
        frooin = ROOT.TFile.Open(temppro.file,"read")
        #frooin.ls()
        try:
            htemp=frooin.Get("hWeights")
            print "entries in hWeights ",htemp.GetEntries()
            weightHistoDict[sample]=copy(htemp)
        except:
            print "SKIPPING FILE NOT WORKING ",sample
            continue
        #frooin.Close()
        temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3]}
        temppro.cuts={sampleDict[sample][0]:""}
        if "ggTo2mu2tau" in sample:
            temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} # worked before
            #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.0001*5.0)} # SM Higgs xsec x BR Haa x 5 for DataMC control plots
            #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.001* 5.0)} # SM Higgs xsec x BR Haa  for signal extraction data MC control plots(AN 17029)
        HAA_processes[temppro.nickname]=temppro


    jetWeightMultiplicity = {}
    try:
        DYJetsFile = ROOT.TFile.Open(HAA_processes["DYJetsToLLext1"].file,"read")
        jetWeightMultiplicity["DYJetsToLLext1"]=DYJetsFile.Get("hWeights").GetSumOfWeights()
        DY1JetsFile = ROOT.TFile.Open(HAA_processes["DY1JetsToLL"].file,"read")
        jetWeightMultiplicity["DY1JetsToLL"]=DY1JetsFile.Get("hWeights").GetSumOfWeights()
        DY2JetsFile = ROOT.TFile.Open(HAA_processes["DY2JetsToLL"].file,"read")
        jetWeightMultiplicity["DY2JetsToLL"]=DY2JetsFile.Get("hWeights").GetSumOfWeights()
        DY3JetsFile = ROOT.TFile.Open(HAA_processes["DY3JetsToLL"].file,"read")
        jetWeightMultiplicity["DY3JetsToLL"]=DY3JetsFile.Get("hWeights").GetSumOfWeights()
        DY4JetsFile = ROOT.TFile.Open(HAA_processes["DY4JetsToLL"].file,"read")
        jetWeightMultiplicity["DY4JetsToLL"]=DY4JetsFile.Get("hWeights").GetSumOfWeights()
    #
        WJetsFile = ROOT.TFile.Open(HAA_processes["WJetsToLNu"].file,"read")
        jetWeightMultiplicity["WJetsToLNu"]=WJetsFile.Get("hWeights").GetSumOfWeights()
        W1JetsFile = ROOT.TFile.Open(HAA_processes["W1JetsToLNu"].file,"read")
        jetWeightMultiplicity["W1JetsToLNu"]=W1JetsFile.Get("hWeights").GetSumOfWeights()
        W2JetsFile = ROOT.TFile.Open(HAA_processes["W2JetsToLNu"].file,"read")
        jetWeightMultiplicity["W2JetsToLNu"]=W2JetsFile.Get("hWeights").GetSumOfWeights()
        W3JetsFile = ROOT.TFile.Open(HAA_processes["W3JetsToLNu"].file,"read")
        jetWeightMultiplicity["W3JetsToLNu"]=W3JetsFile.Get("hWeights").GetSumOfWeights()
        W4JetsFile = ROOT.TFile.Open(HAA_processes["W4JetsToLNu"].file,"read")
        jetWeightMultiplicity["W4JetsToLNu"]=W4JetsFile.Get("hWeights").GetSumOfWeights()
    except:
        print "check jetWeightMultiplicity ... the files may not be loaded"


    Bkg = ["DY","W","TT","ST","EWK"]
    irBkg = ["ZZ","ZHToTauTau","vbf","WHTT"]
    TrialphaBkg = ["ttZ","ttW","WWZ","WZZ","ZZZ","WWW_4F","HZJ"]
    rareBkg = ["Other","rare","WZ"]
    finalDistributions = {}
    finalDistributions["Bkg"]=Bkg
    finalDistributions["data_obs"]=["data_obs"]
    finalDistributions["a15"]=["a15"]
    finalDistributions["a20"]=["a20"]
    finalDistributions["a25"]=["a25"]
    finalDistributions["a30"]=["a30"]
    finalDistributions["a35"]=["a35"]
    finalDistributions["a40"]=["a40"]
    finalDistributions["a45"]=["a45"]
    finalDistributions["a50"]=["a50"]
    finalDistributions["a55"]=["a55"]
    finalDistributions["a60"]=["a60"]
    finalDistributions["irBkg"]=irBkg
    finalDistributions["TrialphaBkg"]=TrialphaBkg
    finalDistributions["rareBkg"]=rareBkg

    allSkims = {}
    finalSkims ={}


    #EventWeights = getEventWeightDicitonary()

    return allcats, HAA_processes,finalDistributions,weightHistoDict,jetWeightMultiplicity





def info(title):
    print(title)
    print('module name:', __name__)
    print('parent process:', os.getppid())
    print('process id:', os.getpid())

def fstar(args):
    print("arggggg matey!")
    return f(*args)

def slimskimstar(args):
    return slimskim(*args)

def f(process, categories):
    info('function f')
    print('hello', process.nickname)
    for cat in categories.keys():
        print("category ",cat)
    return 1

def makeCutsOnTreeArray(process, masterArray,allcats,weightHistoDict,systematic):

    from utils.functions import functs
    from utils.Weights import CommonWeights
    from ROOT import gInterpreter
    commonweight = CommonWeights["lumi"][0]
    skimArrayPerCat = {}
    print "working on process obj ",process.nickname
    for cat in allcats.keys():
        #print(process)
        #print(process.cuts)
        #print(masterArray.keys())
        plottedVars = []

        if process.nickname=="data_obs":
            masterArray['finalweight']=np.full(len(masterArray['evt']),1.0)

            for variableHandle in allcats[cat].vars.keys():
                variable = allcats[cat].vars[variableHandle][0]
                if "[" in variable:
                    #print "adding jagged variable to array ",variable
                    basevar = variable.split("[")[0]
                    index = int(variable.split("[")[1].split("]")[0])
                    val = masterArray[basevar][:,index]
                    #masterArray[basevar+"_"+str(index)]=val
                    #plottedVars.append(basevar+"_"+str(index))
                    masterArray[variableHandle+"_"+str(index)]=val
                    plottedVars.append(variableHandle+"_"+str(index))
                else:
                    plottedVars.append(variableHandle)


            for var in allcats[cat].newvariables.keys():
                newVarVals[var]=0.0

            cuts=[]
            for cuttype in allcats[cat].cuts.keys():
                for cut in allcats[cat].cuts[cuttype]:
                    cuts.append(cut)


            for var in newVarVals.keys():
                arguments = allcats[cat].newvariables[var][2]
                tempvals=[]
                for ag in arguments:
                    tempvals.append(returnArray(masterArray,ag))
                masterArray[var] = functs[allcats[cat].newvariables[var][0]](*tempvals)

            mask = cutOnArray(masterArray,cuts)
            masterArray["mask"]=mask
            masterArray["finalweight"] *= mask.astype(int)
            weightfinal = 1.0   #don't weight the data!!

            skipEvents = np.where(mask==0)[0]
            skimArray={}
            print("before skim", len(masterArray["finalweight"]))
            for key in masterArray.keys():
                try:
                    skimArray[key] = masterArray[key][mask]
                except:
                    #print "length problem? length of key in master ",len(masterArray[key])," length of mask ",len(mask)
                    #print "skipping branch ",key
                    continue
            print("after skim", len(skimArray["mll"]), process.file)
            if len(skimArray["mll"])==0:
                continue

            for key in skimArray.keys():
                if key not in plottedVars and key != "finalweight":
                    del skimArray[key]
            print "working on category ",cat
            skimArrayPerCat[systematic+":"+cat+":"+process.nickname+":"+process.cuts.keys()[0]] = skimArray




        if(process.nickname=="FF" and datadrivenPackage["bool"]):
            masterArray['finalweight']=np.full(len(masterArray['evt']),1.0)

            for variableHandle in allcats[cat].vars.keys():
                variable = allcats[cat].vars[variableHandle][0]
                if "[" in variable:
                    #print "adding jagged variable to array ",variable
                    basevar = variable.split("[")[0]
                    index = int(variable.split("[")[1].split("]")[0])
                    val = masterArray[basevar][:,index]
                    #masterArray[basevar+"_"+str(index)]=val
                    #plottedVars.append(basevar+"_"+str(index))
                    masterArray[variableHandle+"_"+str(index)]=val
                    plottedVars.append(variableHandle+"_"+str(index))
                else:
                    plottedVars.append(variableHandle)


            for var in allcats[cat].newvariables.keys():
                newVarVals[var]=0.0

            cuts=[]
            for cuttype in allcats[cat].cuts.keys():
                for cut in allcats[cat].cuts[cuttype]:
                    cuts.append(cut)


            for var in newVarVals.keys():
                arguments = allcats[cat].newvariables[var][2]
                tempvals=[]
                for ag in arguments:
                    tempvals.append(returnArray(masterArray,ag))
                masterArray[var] = functs[allcats[cat].newvariables[var][0]](*tempvals)

            tempmask=np.full(len(masterArray["evt"]),1.0)

            #events that pass the FF_1 criteria
            tempmask_1 = cutOnArray(masterArray,HAA_processes["FF"].cuts["FF_1"])
            #print tempmask_1[:1000]
            tempmask_2 = cutOnArray(masterArray,HAA_processes["FF"].cuts["FF_2"])
            #print tempmask_2[:1000]
            tempmask_12 = cutOnArray(masterArray,HAA_processes["FF"].cuts["FF_12"])
            #print tempmask_12[:1000]


            #FF_1
            #causal catching ... the pt may be outside the shape from the histogram... if so we need the constant fit value for extrapolation
            fitmask_1 = cutOnArray(masterArray,[["pt_3","<",datadrivenPackage["fakerate1"].GetBinLowEdge(datadrivenPackage["fakerate1"].GetNbinsX())],["pt_3",">",datadrivenPackage["fakerate1"].GetBinLowEdge(2)]])
            fitmask_1 = fitmask_1.astype(int)
            ptarr_1 = masterArray["pt_3"]

            ffweight_1 = ptFun(datadrivenPackage["fakerate1"],ptarr_1)
            ffweight_1 = ffweight_1/(1.0000000001 - ffweight_1)



            #FF_2
            fitmask_2 = cutOnArray(masterArray,[["pt_4","<",datadrivenPackage["fakerate2"].GetBinLowEdge(datadrivenPackage["fakerate2"].GetNbinsX())],["pt_4",">",datadrivenPackage["fakerate2"].GetBinLowEdge(2)]])
            fitmask_2 = fitmask_2.astype(int)
            ptarr_2 = masterArray["pt_4"]

            ffweight_2 = ptFun(datadrivenPackage["fakerate2"],ptarr_2)
            ffweight_2 = ffweight_2/(1.0000000001 - ffweight_2)


            ffweight = ffweight_1 * ffweight_2

            #replace 0s with constant fit value
            ffweight_1 *= fitmask_1
            ffweight_1[np.where(ffweight_1==0)] =  datadrivenPackage["fitrate1"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate1"].GetParameter(0))
            ffweight_2 *= fitmask_2
            ffweight_2[np.where(ffweight_2==0)] =  datadrivenPackage["fitrate2"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate2"].GetParameter(0))

            #FF_12
            #replace 0s with constant fit value
            ffweight[np.where(ffweight==0)] =  (datadrivenPackage["fitrate1"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate1"].GetParameter(0)))*(datadrivenPackage["fitrate2"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate2"].GetParameter(0)))

            fitmask_1 *= fitmask_2
            ffweight *= fitmask_1
            ffweight *= tempmask_12

            ffweight_1[~tempmask_1] = 0.0
            #print "ffweight_1 ",ffweight_1[:1000]
            ffweight_2[~tempmask_2] = 0.0
            finalWeight = ffweight_1 + ffweight_2
            intersection = np.all((tempmask_1,tempmask_2), axis=0)

            finalWeight -= ffweight


            #print "check on fake factor final weight ",finalWeight[:1000]


            masterArray["finalweight"] *= finalWeight
            print "summed final weight ",np.sum(finalWeight)

            keepEvents = ~np.where(finalWeight==0.0)[0]

            skimArray={}
            for key in masterArray.keys():
                skimArray[key] = masterArray[key][keepEvents]

            print("after skim", len(skimArray["mll"]), process.file)
            if len(skimArray["mll"])==0:
                continue

            for key in skimArray.keys():
                if key not in plottedVars and key != "finalweight":
                    del skimArray[key]
            #return skimArray
            skimArrayPerCat[systematic+":"+cat+":"+process.nickname+":"+process.cuts.keys()[0]] = skimArray




        if process.nickname not in ["data_obs","FF","FF_1","FF_2","FF_12"]:
            EventWeights = getEventWeightDicitonary()

            newVarVals={}
            masterArray['finalweight']=np.full(len(masterArray['evt']),1.0)

            for variableHandle in allcats[cat].vars.keys():
                variable = allcats[cat].vars[variableHandle][0]
                if "[" in variable:
                    #print "adding jagged variable to array ",variable
                    basevar = variable.split("[")[0]
                    index = int(variable.split("[")[1].split("]")[0])
                    val = masterArray[basevar][:,index]
                    #masterArray[basevar+"_"+str(index)]=val
                    #plottedVars.append(basevar+"_"+str(index))
                    masterArray[variableHandle+"_"+str(index)]=val
                    plottedVars.append(variableHandle+"_"+str(index))
                else:
                    plottedVars.append(variableHandle)


            for var in allcats[cat].newvariables.keys():
                newVarVals[var]=0.0

            cuts=[]
            for cuttype in allcats[cat].cuts.keys():
                for cut in allcats[cat].cuts[cuttype]:
                    cuts.append(cut)


            for var in newVarVals.keys():
                arguments = allcats[cat].newvariables[var][2]
                tempvals=[]
                for ag in arguments:
                    tempvals.append(returnArray(masterArray,ag))
                masterArray[var] = functs[allcats[cat].newvariables[var][0]](*tempvals)
                plottedVars.append(var)

            mask = cutOnArray(masterArray,cuts)
            masterArray["mask"]=mask
            masterArray["finalweight"] *= mask.astype(int)
            weightfinal = 1.0   #don't weight the data!!

            skipEvents = np.where(mask==0)[0]
            skimArray={}
            if len(masterArray["mll"])==0:
                continue

            weightDict = process.weights
            weightfinal = commonweight
            nickname = process.nickname
            for scalefactor in weightDict.keys():
                if scalefactor == "kfactor":
                    weightfinal =  weightfinal * (1 / float(weightDict[scalefactor]))
                elif scalefactor in ["PU"]:
                    masterArray["finalweight"] *= (returnArray(masterArray,weightDict[scalefactor]))
                elif scalefactor =="theoryXsec":
                    weightfinal =  weightfinal * float(weightDict[scalefactor])

            #print "finalweight after PU and kFactor ",masterArray["finalweight"][:100]

            if nickname =="DY1JetsToLL":
                norm1 = jetWeightMultiplicity["DYJetsToLLext1"]/HAA_processes["DYJetsToLLext1"].weights["xsec"]
                norm2 = jetWeightMultiplicity["DY1JetsToLL"]/HAA_processes["DY1JetsToLL"].weights["xsec"]
                weightfinal = weightfinal * 1/(norm1+norm2)
            if nickname =="DY2JetsToLL":
                norm1 = jetWeightMultiplicity["DYJetsToLLext1"]/HAA_processes["DYJetsToLLext1"].weights["xsec"]
                norm2 = jetWeightMultiplicity["DY2JetsToLL"]/HAA_processes["DY2JetsToLL"].weights["xsec"]
                weightfinal = weightfinal * 1/(norm1+norm2)
            if nickname =="DY3JetsToLL":
                norm1 = jetWeightMultiplicity["DYJetsToLLext1"]/HAA_processes["DYJetsToLLext1"].weights["xsec"]
                norm2 = jetWeightMultiplicity["DY3JetsToLL"]/HAA_processes["DY3JetsToLL"].weights["xsec"]
                weightfinal = weightfinal * 1/(norm1+norm2)
            if nickname =="DY4JetsToLL":
                norm1 = jetWeightMultiplicity["DYJetsToLLext1"]/HAA_processes["DYJetsToLLext1"].weights["xsec"]
                norm2 = jetWeightMultiplicity["DY4JetsToLL"]/HAA_processes["DY4JetsToLL"].weights["xsec"]
                weightfinal = weightfinal * 1/(norm1+norm2)
            if nickname =="W1JetsToLNu":
                norm1 = jetWeightMultiplicity["WJetsToLNu"]/HAA_processes["WJetsToLNu"].weights["xsec"]
                norm2 = jetWeightMultiplicity["W1JetsToLNu"]/HAA_processes["W1JetsToLNu"].weights["xsec"]
                weightfinal = weightfinal * 1/(norm1+norm2)
            elif nickname =="W2JetsToLNu":
                norm1 = jetWeightMultiplicity["WJetsToLNu"]/HAA_processes["WJetsToLNu"].weights["xsec"]
                norm2 = jetWeightMultiplicity["W2JetsToLNu"]/HAA_processes["W2JetsToLNu"].weights["xsec"]
                weightfinal = weightfinal * 1/(norm1+norm2)
            elif nickname =="W3JetsToLNu":
                norm1 = jetWeightMultiplicity["WJetsToLNu"]/HAA_processes["WJetsToLNu"].weights["xsec"]
                norm2 = jetWeightMultiplicity["W3JetsToLNu"]/HAA_processes["W3JetsToLNu"].weights["xsec"]
                weightfinal = weightfinal * 1/(norm1+norm2)
            else:
                tempDenom=1.0
                sumOfWeights = 1.0
                for nic in weightHistoDict.keys():
                    if nickname==nic:
                        sumOfWeights = weightHistoDict[nickname].GetSumOfWeights()
                    if nickname==nic.strip("_ext")[0]:
                        sumOfWeights*=tempDenom*weightHistoDict[nic].GetSumOfWeights()
                weightfinal = weightfinal * HAA_processes[nickname].weights["xsec"]/ sumOfWeights
                print "xsec/SoW ",HAA_processes[nickname].weights["xsec"]/ sumOfWeights

            #multiply by scalar weight
            masterArray["finalweight"] *= weightfinal
            print "finalweight before per event scaling ",masterArray["finalweight"][:100]


            #eventWeightDict = process.eventWeights
            eventWeightDict = EventWeights
            if eventWeightDict:
                for scalefactor in eventWeightDict.keys():

                    cutlist = eventWeightDict[scalefactor][0]
                    weightMask=cutOnArray(masterArray,cutlist)
                    weightMask=weightMask.astype(float)

                    if type(eventWeightDict[scalefactor][1][0])==float:
                        weightMask *= eventWeightDict[scalefactor][1][0]

                    if hasattr(eventWeightDict[scalefactor][1][0],'__call__'):
                        arguments = eventWeightDict[scalefactor][1][1]
                        tempvals=[]

                        for ag in arguments:
                           tempvals.append(returnArray(masterArray,ag))
                        weightMask*= eventWeightDict[scalefactor][1][0](*tempvals)

                    if scalefactor!="fake" and scalefactor!="fake1" and scalefactor!="fake2":
                        weightMask[np.where(weightMask==0.0)]=1.0
                    else:
                        print "subtracting fakes "
                        print weightMask[:100]

                    masterArray["finalweight"] *= weightMask
            print "finalweight after per event scaling ",masterArray["finalweight"][:100]

            print("before skim", len(masterArray["finalweight"]))
            for key,value in masterArray.iteritems():
                if ( key in plottedVars or key == "finalweight" ) \
                    and (len(mask)==len(value)):
                    skimArray[key] = value[mask]
            #for key in masterArray.keys():
                #try:
                    #skimArray[key] = masterArray[key][mask]
                #    skimArray[key] = value[mask]
                #except:
                    #print "length problem? length of key in master ",len(masterArray[key])," length of mask ",len(mask)
                    #print "skipping branch ",key
                #    continue
            print("after skim", len(skimArray["mll"]), process.file)
            #for key in skimArray.keys():
            #    if key not in plottedVars and key != "finalweight":
            #        del skimArray[key]
            print "working on category ",cat
            skimArrayPerCat[systematic+":"+cat+":"+process.nickname+":"+process.cuts.keys()[0]] = skimArray
    return skimArrayPerCat





def slimskim(process,allcats,weightHistoDict):
    import multiprocessing as mp

    skimArrayPerSysCats={}
    fin = uproot.open(process.file)
    try:
        tree = fin["Events"]
    except:
        return skimArrayPerSysCats






    skimArrayPerSysCats.update(makeCutsOnTreeArray(process,tree.arrays(),allcats,weightHistoDict,"Nominal"))
    print skimArrayPerSysCats.keys()
    try:
        print "final weight after nominal for Nominal:mmtt_inclusive:ZZTo4L:ZZ",skimArrayPerSysCats["Nominal:mmtt_inclusive:ZZTo4L:ZZ"]["finalweight"][:100]
    except:
        print "Nominal:mmtt_inclusive:ZZTo4L:ZZ doesn't exist for sample"
    #treatment of systematics

    systematic = "scale_m_etalt1p2Up"
    work_dict = dict()
    try:
        tree = fin[systematic]
    except:
        return skimArrayPerSysCats
    syst_names = set(fin[systematic].keys())
    nom_names = set(fin["Events"].keys()) - syst_names
    work_dict.update(fin[systematic].arrays(list(syst_names)))
    work_dict.update(fin["Events"].arrays(list(nom_names)))

    skimArrayPerSysCats.update(makeCutsOnTreeArray(process,work_dict,allcats,weightHistoDict,systematic))
    print skimArrayPerSysCats.keys()


    #pool  = mp.Pool(2)
    ##nums=pool.map(fstar,payloads)
#
    ##skims = slimskimstar(payloads[0])
    #sytematics = pool.map(slimskimstar,payloads)
#
    #pool.close()
    #pool.join()
    ##print(nums)





    return skimArrayPerSysCats












def createOutputSystematics(skimmedArraysSet,finalDistributions):
    import root_numpy

    finalSkims={}
    cats = []
    systematics = []
    #print skimmedArrays

    #print skimmedArrays.keys()

    for skimmedArrays in skimmedArraysSet:
        for sysCatNicPro in skimmedArrays.keys():
            print "key ",sysCatNicPro
            sys = sysCatNicPro.split(":")[0]
            systematics.append(sys)
            channel = sysCatNicPro.split(":")[1].split("_")[0]
            print "channel ",channel
            print "systematic ",sys
            cat = sysCatNicPro.split(":")[1]
            cats.append(cat)
            if sys in finalSkims:
                if not cat in finalSkims[sys]:
                    finalSkims[sys][cat]={}
            else:
                finalSkims[sys]={}
                if not cat in finalSkims[sys]:
                    finalSkims[sys][cat]={}
            process= sysCatNicPro.split(":")[-1]
            print "cat ",cat
            print "process ",process

            for catDist in finalDistributions.keys():
                for processOut in finalDistributions[catDist]:
                    #if (catDist not in finalSkims[sys][cat]):
                    #print " is catDist ",catDist,"  inside keys ? ",finalSkims[sys][cat].keys()
                    if (processOut==process) and (catDist not in finalSkims[sys][cat].keys()):
                        #finalSkims[sys][cat][catDist] = processSkims[process]
                        #print "first output to finalskims ", sysCatNicPro,"  for process ",process," finalDist cat ",catDist
                        finalSkims[sys][cat][catDist] = skimmedArrays[sysCatNicPro]
                    #elif (catDist in finalSkims[sys][cat]):
                    elif (processOut==process) and (catDist in finalSkims[sys][cat].keys()):
                        #print "adding to finalskims ", sysCatNicPro,"  for process ",process," finalDist cat ",catDist
                        for branch in finalSkims[sys][cat][catDist].keys():
                            #finalSkims[sys][cat][catDist][branch]=np.concatenate((finalSkims[sys][cat][catDist][branch],processSkims[process][branch]))
                            finalSkims[sys][cat][catDist][branch]=np.concatenate((finalSkims[sys][cat][catDist][branch],skimmedArrays[sysCatNicPro][branch]))
                    else:
                        continue

    #print finalSkims.keys()
    #print finalSkims["Nominal"].keys()

    for cat in cats:
        localtest="sys"
        skimFile = ROOT.TFile("skimmed_"+localtest+"_"+cat+".root","recreate")
        skimFile.cd()
        for sys in finalSkims.keys():
            dataTypes =[[],[]]
            #print finalSkims[sys][cat].values()
            random_sample = finalSkims[sys][cat].values()[0]
            for branch in random_sample.keys():
                dataTypes[0].append(branch)
                dataTypes[1].append(random_sample[branch].dtype)
            for catDist in finalSkims[sys][cat].keys():
                #print "on the final dist ",catDist
                data = np.zeros(len(finalSkims[sys][cat][catDist][branch]),dtype={'names':dataTypes[0],'formats':dataTypes[1]})
                for branch in data.dtype.names:
                    #print "working on branch ",branch
                    #print "branch datatype ",type(finalSkims[sys][cat][catDist][branch])
                    if len(finalSkims[sys][cat][catDist][branch].shape) == 1:   # flat array important for jagged arrays of input data
                        data[branch] = finalSkims[sys][cat][catDist][branch]
                    else:
                        data[branch] = finalSkims[sys][cat][catDist][branch][:,0]
                treeOut = root_numpy.array2tree(data, name=sys+"_"+catDist)
                treeOut.Write()
        skimFile.Close()

    return 1

def slimskimorg(process,allcats,weightHistoDict):

    from utils.functions import functs
    from utils.Weights import CommonWeights
    from ROOT import gInterpreter
    commonweight = CommonWeights["lumi"][0]
    fin = uproot.open(process.file)
    tree = fin["Events"]
    #print(fin)
    #print(tree)

    #for cat, catInfo in allcats.iteritems():
    #    for variableHandle in catInfo.vars.keys():
    for cat in allcats.keys():
        masterArray = tree.arrays()
        #print(process)
        #print(process.cuts)
        #print(masterArray.keys())
        plottedVars = []

        EventWeights = getEventWeightDicitonary()

        newVarVals={}
        skimArrayPerCat = {}
        masterArray['finalweight']=np.full(len(masterArray['evt']),1.0)

        for variableHandle in allcats[cat].vars.keys():
            variable = allcats[cat].vars[variableHandle][0]
            if "[" in variable:
                #print "adding jagged variable to array ",variable
                basevar = variable.split("[")[0]
                index = int(variable.split("[")[1].split("]")[0])
                val = masterArray[basevar][:,index]
                #masterArray[basevar+"_"+str(index)]=val
                #plottedVars.append(basevar+"_"+str(index))
                masterArray[variableHandle+"_"+str(index)]=val
                plottedVars.append(variableHandle+"_"+str(index))
            else:
                plottedVars.append(variableHandle)


        for var in allcats[cat].newvariables.keys():
            newVarVals[var]=0.0

        cuts=[]
        for cuttype in allcats[cat].cuts.keys():
            for cut in allcats[cat].cuts[cuttype]:
                cuts.append(cut)


        for var in newVarVals.keys():
            arguments = allcats[cat].newvariables[var][2]
            tempvals=[]
            for ag in arguments:
                tempvals.append(returnArray(masterArray,ag))
            masterArray[var] = functs[allcats[cat].newvariables[var][0]](*tempvals)

        mask = cutOnArray(masterArray,cuts)
        masterArray["mask"]=mask
        masterArray["finalweight"] *= mask.astype(int)
        weightfinal = 1.0   #don't weight the data!!
        print "finalweight right after mask ",masterArray["finalweight"][:100]

        weightDict = process.weights
        #print " weight dicitonary ? ",weightDict
        print "sum of weights dicitonary ? ",weightHistoDict
        weightfinal = commonweight
        nickname = process.nickname
        print "working on the process nickname ",nickname
        print "sum of weights dicitonary keys ? ",weightHistoDict.keys()
        for scalefactor in weightDict.keys():
            if scalefactor == "kfactor":
                weightfinal =  weightfinal * (1 / float(weightDict[scalefactor]))
            elif scalefactor in ["PU"]:
                masterArray["finalweight"] *= (returnArray(masterArray,weightDict[scalefactor]))
            elif scalefactor =="theoryXsec":
                weightfinal =  weightfinal * float(weightDict[scalefactor])

        print "finalweight after PU and kFactor ",masterArray["finalweight"][:100]

        if nickname =="DY1JetsToLL":
            norm1 = jetWeightMultiplicity["DYJetsToLLext1"]/HAA_processes["DYJetsToLLext1"].weights["xsec"]
            norm2 = jetWeightMultiplicity["DY1JetsToLL"]/HAA_processes["DY1JetsToLL"].weights["xsec"]
            weightfinal = weightfinal * 1/(norm1+norm2)
        if nickname =="DY2JetsToLL":
            norm1 = jetWeightMultiplicity["DYJetsToLLext1"]/HAA_processes["DYJetsToLLext1"].weights["xsec"]
            norm2 = jetWeightMultiplicity["DY2JetsToLL"]/HAA_processes["DY2JetsToLL"].weights["xsec"]
            weightfinal = weightfinal * 1/(norm1+norm2)
        if nickname =="DY3JetsToLL":
            norm1 = jetWeightMultiplicity["DYJetsToLLext1"]/HAA_processes["DYJetsToLLext1"].weights["xsec"]
            norm2 = jetWeightMultiplicity["DY3JetsToLL"]/HAA_processes["DY3JetsToLL"].weights["xsec"]
            weightfinal = weightfinal * 1/(norm1+norm2)
        if nickname =="DY4JetsToLL":
            norm1 = jetWeightMultiplicity["DYJetsToLLext1"]/HAA_processes["DYJetsToLLext1"].weights["xsec"]
            norm2 = jetWeightMultiplicity["DY4JetsToLL"]/HAA_processes["DY4JetsToLL"].weights["xsec"]
            weightfinal = weightfinal * 1/(norm1+norm2)
        if nickname =="W1JetsToLNu":
            norm1 = jetWeightMultiplicity["WJetsToLNu"]/HAA_processes["WJetsToLNu"].weights["xsec"]
            norm2 = jetWeightMultiplicity["W1JetsToLNu"]/HAA_processes["W1JetsToLNu"].weights["xsec"]
            weightfinal = weightfinal * 1/(norm1+norm2)
        elif nickname =="W2JetsToLNu":
            norm1 = jetWeightMultiplicity["WJetsToLNu"]/HAA_processes["WJetsToLNu"].weights["xsec"]
            norm2 = jetWeightMultiplicity["W2JetsToLNu"]/HAA_processes["W2JetsToLNu"].weights["xsec"]
            weightfinal = weightfinal * 1/(norm1+norm2)
        elif nickname =="W3JetsToLNu":
            norm1 = jetWeightMultiplicity["WJetsToLNu"]/HAA_processes["WJetsToLNu"].weights["xsec"]
            norm2 = jetWeightMultiplicity["W3JetsToLNu"]/HAA_processes["W3JetsToLNu"].weights["xsec"]
            weightfinal = weightfinal * 1/(norm1+norm2)
        else:
            tempDenom=1.0
            sumOfWeights = 1.0
            for nic in weightHistoDict.keys():
                if nickname==nic:
                    sumOfWeights = weightHistoDict[nickname].GetSumOfWeights()
                if nickname==nic.strip("_ext")[0]:
                    sumOfWeights*=tempDenom*weightHistoDict[nic].GetSumOfWeights()
            weightfinal = weightfinal * HAA_processes[nickname].weights["xsec"]/ sumOfWeights
            print "xsec/SoW ",HAA_processes[nickname].weights["xsec"]/ sumOfWeights

        #multiply by scalar weight
        masterArray["finalweight"] *= weightfinal
        print "finalweight before per event scaling ",masterArray["finalweight"][:100]


        eventWeightDict = process.eventWeights
        if eventWeightDict:
            for scalefactor in eventWeightDict.keys():

                cutlist = eventWeightDict[scalefactor][0]
                weightMask=cutOnArray(masterArray,cutlist)
                weightMask=weightMask.astype(float)

                if type(eventWeightDict[scalefactor][1][0])==float:
                    weightMask *= eventWeightDict[scalefactor][1][0]

                if hasattr(eventWeightDict[scalefactor][1][0],'__call__'):
                    arguments = eventWeightDict[scalefactor][1][1]
                    tempvals=[]

                    for ag in arguments:
                       tempvals.append(returnArray(masterArray,ag))
                    weightMask*= eventWeightDict[scalefactor][1][0](*tempvals)

                if scalefactor!="fake" and scalefactor!="fake1" and scalefactor!="fake2":
                    weightMask[np.where(weightMask==0.0)]=1.0
                else:
                    print "subtracting fakes "
                    print weightMask[:100]

                masterArray["finalweight"] *= weightMask
        print "finalweight after per event scaling ",masterArray["finalweight"][:100]
        skipEvents = np.where(mask==0)[0]
        skimArray={}
        print("before skim", len(masterArray["finalweight"]))
        for key in masterArray.keys():
            skimArray[key] = masterArray[key][mask]
        print("after skim", len(skimArray["mll"]))

        for key in skimArray.keys():
            if key not in plottedVars and key != "finalweight":
                del skimArray[key]
        skimArrayPerCat[cat+":"+process.nickname+":"+process.cuts.keys()[0]] = skimArray


    return skimArrayPerCat

def createOutput(skimmedArrays,finalDistributions):
    import root_numpy

    finalSkims={}
    cats = []
    for skimArrayPerCat in skimmedArrays:
        #print skimArrayPerCat.keys()
        for catNicPro in skimArrayPerCat.keys():
            print "key ",catNicPro
            channel = catNicPro.split(":")[0].split("_")[0]
            print "channel ",channel
            cat = catNicPro.split(":")[0]
            cats.append(cat)
            try:
                finalSkims[cat]={}
            except:
                continue
            process= catNicPro.split(":")[-1]
            print "cat ",cat
            print "process ",process

            for catDist in finalDistributions.keys():
                for processOut in finalDistributions[catDist]:
                    #if (catDist not in finalSkims[cat]):
                    if (processOut==process) and (catDist not in finalSkims[cat]):
                        #finalSkims[cat][catDist] = processSkims[process]
                        print "first output to finalskims ", catNicPro,"  for process ",process," finalDist cat ",catDist
                        finalSkims[cat][catDist] = skimArrayPerCat[catNicPro]
                        continue
                    #elif (catDist in finalSkims[cat]):
                    elif (processOut==process) and (catDist in finalSkims[cat]):
                        print "adding to finalskims ", catNicPro,"  for process ",process," finalDist cat ",catDist
                        for branch in finalSkims[cat][catDist].keys():
                            #finalSkims[cat][catDist][branch]=np.concatenate((finalSkims[cat][catDist][branch],processSkims[process][branch]))
                            finalSkims[cat][catDist][branch]=np.concatenate((finalSkims[cat][catDist][branch],skimArrayPerCat[catNicPro][branch]))
                    else:
                        continue


    for cat in cats:
        localtest="ZZ"
        skimFile = ROOT.TFile("skimmed_"+localtest+"_"+cat+".root","recreate")
        skimFile.cd()

        dataTypes =[[],[]]
        #print finalSkims[cat].values()
        random_sample = finalSkims[cat].values()[0]
        for branch in random_sample.keys():
            dataTypes[0].append(branch)
            dataTypes[1].append(random_sample[branch].dtype)
        for catDist in finalSkims[cat].keys():
            #print "on the final dist ",catDist
            data = np.zeros(len(finalSkims[cat][catDist][branch]),dtype={'names':dataTypes[0],'formats':dataTypes[1]})
            for branch in data.dtype.names:
                #print "working on branch ",branch
                #print "branch datatype ",type(finalSkims[cat][catDist][branch])
                if len(finalSkims[cat][catDist][branch].shape) == 1:   # flat array important for jagged arrays of input data
                    data[branch] = finalSkims[cat][catDist][branch]
                else:
                    data[branch] = finalSkims[cat][catDist][branch][:,0]
            treeOut = root_numpy.array2tree(data, name=catDist)
            treeOut.Write()
        skimFile.Close()

    return 1


if __name__ == "__main__":
    import datetime
    begin_time = datetime.datetime.now()
    from multiprocessing import *
    import multiprocessing as mp
    import logging
    import os
    allcats={}
    HAA_processes={}
    finalDistributions={}
    filelist = {}
    weightHistoDict={}
    jetWeightMultiplicity={}
    EventWeights={}


    allcats,HAA_processes,finalDistributions,weightHistoDict,jetWeightMultiplicity = initialize()

    #print(allcats)
    info('main line')
    #print(f)
    nums=[]
    skims=[]
    #skims={}
    payloadsdict={}
    payloads=[]

    for nickname, process in HAA_processes.items():
        payloadsdict[nickname]=[process,allcats]
        payloads.append((process,allcats,weightHistoDict))

    #print(payloads)


    #slimskimstar(payloads[0])
    pool  = mp.Pool(12)
    #nums=pool.map(fstar,payloads)

    #skims = slimskimstar(payloads[0])
    skims = pool.map(slimskimstar,payloads)
    #with mp.Pool(5) as pool:
        #p = mp.Process(target=f, args=('bob',))
        #p.start()
        #p.join()
        #print(pool.map(f,'bob'))
        #nums=pool.map(f,iterargs)
        #nums=pool.map(fstar,payloads)
        #print(payloads)

        #skims = pool.map(slimskimstar,payloads)

        #skims[nickname]=pool.map(slimskim,)
        #pool.map(f,'bob')

    pool.close()
    pool.join()
    #print(nums)
    #print skims
    #for item in skims:
    #    print(item.keys())

    print("combining the output")
    #if createOutput(skims,finalDistributions) : print "successful output"
    #else: print "ouch ..."
    if createOutputSystematics(skims,finalDistributions) : print "successful output"
    else: print "ouch ..."



    print("computation time")
    print(datetime.datetime.now() - begin_time)