#########################
#Author: Sam Higginbotham
'''

* File Name : Processes.py

* Purpose : Stores dictionaries for analysis signal and background files 
            THIS DRIVES ALL HISTOGRAMS STORED! AND IN PLOTS FOR THE STACK

* Creation Date : 04-02-2020

* Last Modified :

'''
#########################
import ROOT
import csv
import sys
from Parametrization import Process
from functions import *
from TauPOG.TauIDSFs.TauIDSFTool import TauIDSFTool
from TauPOG.TauIDSFs.TauIDSFTool import TauESTool
from TauPOG.TauIDSFs.TauIDSFTool import TauFESTool
#sys.path.append('../SFs/')
import SFs.ScaleFactor as SF

#This is a list of Process objects
HAA_processes={}
HAA_signals={}
HAA_processes_test={}

#Gathering extra information

#Gather the analysis datasets and info 
sampleDict = {}

csvfile="/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/AnalysisVisualization/MCsamples_2016_v6.csv"
#csvfile="/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/AnalysisVisualization/MCsamples_2016_v6_sam.csv"
#with open("/afs/cern.ch/work/s/shigginb/cmssw/HAA/HAA_Run2_May2020_10_2_9/src/AnalysisVisualization/MCsamples_2016_v6.csv")  as csvfile:
#    reader = csv.reader(csvfile, delimiter=',')
#    for row in reader:
#        #[nickname]        = [category,xsec,numberOfEvents,finishedEvents,idk?,DASDataset]
#        if len(row.split(',')[2].split("*"))>1:
#            tempval=1.0
#            for val in row.split(',')[2].split("*"):
#                tempval = tempval * float(val)
#            sampleDict[row[0]] = [row[1],tempval,row[3],row[4],row[5],row[6]]
#        else:
#            sampleDict[row[0]] = [row[1],float(row[2]),row[3],row[4],row[5],row[6]]

for line in open(csvfile,'r').readlines() :
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

#Complex set of event weights
# tau id scale factor from object
tauIDSF = TauIDSFTool('2016Legacy','DeepTau2017v2p1VSjet','Medium').getSFvsPT
# muon scale factor as function of pt
sf_MuonId = SF.SFs()
sf_MuonId.ScaleFactor("/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/AnalysisVisualization/ScaleFactors/LeptonEffs/Muon/Muon_Run2016_IdIso_0p2.root")
EventWeights={
    #"name":[[if statements],[weight to apply]]
    "3_mt_lt0p4":[[["decayMode_3","==",0],["eta_3","<",0.4]],[0.80]],
    "3_mj_lt0p4":[[["decayMode_3","==",0],["eta_3","<",0.4]],[1.21]],
    "3_mt_0p4to0p8":[[["decayMode_3","==",0],["eta_3",">",0.4],["eta_3","<",0.8]],[0.81]],
    "3_mj_0p4to0p8":[[["decayMode_3","==",0],["eta_3",">",0.4],["eta_3","<",0.8]],[1.11]],
    "3_mt_0p8to1p2":[[["decayMode_3","==",0],["eta_3",">",0.8],["eta_3","<",1.2]],[0.79]],
    "3_mj_0p8to1p2":[[["decayMode_3","==",0],["eta_3",">",0.8],["eta_3","<",1.2]],[1.2]],
    "3_mt_1p2to1p7":[[["decayMode_3","==",0],["eta_3",">",1.2],["eta_3","<",1.7]],[0.68]],
    "3_mj_1p2to1p7":[[["decayMode_3","==",0],["eta_3",">",1.2],["eta_3","<",1.7]],[1.16]],
    "3_mt_1p7to2p3":[[["decayMode_3","==",0],["eta_3",">",1.7],["eta_3","<",2.3]],[0.68]],
    "3_mj_1p7to2p3":[[["decayMode_3","==",0],["eta_3",">",1.7],["eta_3","<",2.3]],[2.25]],

    "4_mt_lt0p4":[[["decayMode_4","==",0],["eta_4","<",0.4]],[0.80]],
    "4_mj_lt0p4":[[["decayMode_4","==",0],["eta_4","<",0.4]],[1.21]],
    "4_mt_0p4to0p8":[[["decayMode_4","==",0],["eta_4",">",0.4],["eta_4","<",0.8]],[0.81]],
    "4_mj_0p4to0p8":[[["decayMode_4","==",0],["eta_4",">",0.4],["eta_4","<",0.8]],[1.11]],
    "4_mt_0p8to1p2":[[["decayMode_4","==",0],["eta_4",">",0.8],["eta_4","<",1.2]],[0.79]],
    "4_mj_0p8to1p2":[[["decayMode_4","==",0],["eta_4",">",0.8],["eta_4","<",1.2]],[1.2]],
    "4_mt_1p2to1p7":[[["decayMode_4","==",0],["eta_4",">",1.2],["eta_4","<",1.7]],[0.68]],
    "4_mj_1p2to1p7":[[["decayMode_4","==",0],["eta_4",">",1.2],["eta_4","<",1.7]],[1.16]],
    "4_mt_1p7to2.3":[[["decayMode_4","==",0],["eta_4",">",1.7],["eta_4","<",2.3]],[0.68]],
    "4_mj_1p7to2.3":[[["decayMode_4","==",0],["eta_4",">",1.7],["eta_4","<",2.3]],[2.25]],

    #bound method ... last input list is func parameters
    "3_tauSF":[[["cat","==",7],["gen_match_3","==",5]],[tauIDSF,["pt_3","gen_match_3"]]],
    "4_tauSF":[[["cat","==",7],["gen_match_4","==",5]],[tauIDSF,["pt_4","gen_match_4"]]],

    "3_tauSF":[[["cat","==",6],["gen_match_3","==",5]],[tauIDSF,["pt_3","gen_match_3"]]],
    "4_tauSF":[[["cat","==",6],["gen_match_4","==",5]],[tauIDSF,["pt_4","gen_match_4"]]]

}

data = Process()
data.nickname = "data"
data.plotname = "data"
data.weights = {"xsec":1.0} #higgs xsec [pb] * 2hdm type Branching ratio
data.cuts = {"data_obs":""}
data.file = "data.root"
data.classification = "data_obs"
data.color = [ROOT.kBlack]    #do i need a list here??? to embed it
HAA_processes[data.nickname]=data

#This is a very special case ... datadriven approach where we have 3 different regions 
#FF_1 leg 1 fails invert muon criteria 
#FF_2 leg 2 fails invert tau criteria 
#FF_12 both leg fail
from Categories import HAA_Inc_mmmt
from Categories import ffworkingpoint
FF = Process()
FF.nickname = "FF"
FF.plotname = "FF"
FF.weights = {"xsec":1.0} 
#FF.cuts = {"FF":"","FF_1":[HAA_Inc_mmmt.cuts["preselection"][0],[["OR"],["iso_3",">",0.15],["mediumId_3","<",1],["isGlobal_3","<",1],["isTracker_3","<",1]],["idDeepTau2017v2p1VSjet_4",">=",ffworkingpoint]],"FF_2":[HAA_Inc_mmmt.cuts["preselection"][0],["mediumId_3",">=",1],["iso_3","<=",0.15],[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]],["idDeepTau2017v2p1VSjet_4","<",ffworkingpoint]],"FF_12":[HAA_Inc_mmmt.cuts["preselection"][0]]}
FF.cuts = {"FF":"","FF_1":[HAA_Inc_mmmt.cuts["preselection"][0],HAA_Inc_mmmt.cuts["trigger"][0],["idDeepTau2017v2p1VSjet_4",">=",ffworkingpoint],[["EQT"],["q_3","q_4"],"mult","<",0]],"FF_2":[HAA_Inc_mmmt.cuts["preselection"][0],HAA_Inc_mmmt.cuts["trigger"][0],["mediumId_3",">=",1],["iso_3","<=",0.15],[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]],[["EQT"],["q_3","q_4"],"mult","<",0]],"FF_12":[HAA_Inc_mmmt.cuts["preselection"][0],HAA_Inc_mmmt.cuts["trigger"][0],[["EQT"],["q_3","q_4"],"mult","<",0]]}
FF.file = "data.root"
FF.classification = "FF"
#FF.color = [ROOT.kBlack]    #do i need a list here??? to embed it
HAA_processes[FF.nickname]=FF


a15 = Process()
a15.nickname = "ggTo2mu2tau_15"
a15.plotname = "a15"
a15.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a15.eventWeights = EventWeights
a15.cuts = {"a15":""}
a15.file = "ggTo2mu2tau_15_2016.root"
a15.classification = sampleDict[a15.nickname][3]
a15.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a15.nickname]=a15
HAA_signals[a15.nickname]=a15

a20 = Process()
a20.nickname = "ggTo2mu2tau_20"
a20.plotname = "a20"
a20.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a20.eventWeights = EventWeights
a20.cuts = {"a20":""}
a20.file = "ggTo2mu2tau_20_2016.root"
a20.classification = sampleDict[a20.nickname][3]
a20.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a20.nickname]=a20
HAA_signals[a20.nickname]=a20

a25 = Process()
a25.nickname = "ggTo2mu2tau_25"
a25.plotname = "a25"
a25.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a25.eventWeights = EventWeights
a25.cuts = {"a25":""}
a25.file = "ggTo2mu2tau_25_2016.root"
a25.classification = sampleDict[a25.nickname][3]
a25.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a25.nickname]=a25
HAA_signals[a25.nickname]=a25

a30 = Process()
a30.nickname = "ggTo2mu2tau_30"
a30.plotname = "a30"
a30.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a30.eventWeights = EventWeights
a30.cuts = {"a30":""}
a30.file = "ggTo2mu2tau_30_2016.root"
a30.classification = sampleDict[a30.nickname][3]
a30.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a30.nickname]=a30
HAA_signals[a30.nickname]=a30

a35 = Process()
a35.nickname = "ggTo2mu2tau_35"
a35.plotname = "a35"
a35.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a35.eventWeights = EventWeights
a35.cuts = {"a35":""}
a35.file = "ggTo2mu2tau_35_2016.root"
a35.classification = sampleDict[a35.nickname][3]
a35.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a35.nickname]=a35
HAA_signals[a35.nickname]=a35

a40 = Process()
a40.nickname = "ggTo2mu2tau_40"
a40.plotname = "a40"
a40.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a40.eventWeights = EventWeights
a40.cuts = {"a40":""}
a40.file = "ggTo2mu2tau_40_2016.root"
a40.classification = sampleDict[a40.nickname][3]
a40.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a40.nickname]=a40
HAA_signals[a40.nickname]=a40

a45 = Process()
a45.nickname = "ggTo2mu2tau_45"
a45.plotname = "a45"
a45.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a45.eventWeights = EventWeights
a45.cuts = {"a45":""}
a45.file = "ggTo2mu2tau_45_2016.root"
a45.classification = sampleDict[a45.nickname][3]
a45.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a45.nickname]=a45
HAA_signals[a45.nickname]=a45

a50 = Process()
a50.nickname = "ggTo2mu2tau_50"
a50.plotname = "a50"
a50.cuts = {"a50":""}
a50.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a50.eventWeights = EventWeights
a50.file = "ggTo2mu2tau_50_2016.root"
a50.classification = sampleDict[a50.nickname][3]
a50.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a50.nickname]=a50
HAA_signals[a50.nickname]=a50

a55 = Process()
a55.nickname = "ggTo2mu2tau_55"
a55.plotname = "a55"
a55.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a55.eventWeights = EventWeights
a55.cuts = {"a55":""}
a55.file = "ggTo2mu2tau_55_2016.root"
a55.classification = sampleDict[a55.nickname][3]
a55.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a55.nickname]=a55
HAA_signals[a55.nickname]=a55

a60 = Process()
a60.nickname = "ggTo2mu2tau_60"
a60.plotname = "a60"
a60.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a60.eventWeights = EventWeights
a60.cuts = {"a60":""}
a60.file = "ggTo2mu2tau_60_2016.root"
a60.classification = sampleDict[a60.nickname][0]
a60.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a60.nickname]=a60
HAA_signals[a60.nickname]=a60

GluGluToContinToZZTo2mu2tau = Process()
GluGluToContinToZZTo2mu2tau.nickname = "GluGluTo2mu2tau"
GluGluToContinToZZTo2mu2tau.plotname = "GluGluTo2mu2tau"
GluGluToContinToZZTo2mu2tau.weights = {"xsec":sampleDict[GluGluToContinToZZTo2mu2tau.nickname][1],"nevents":sampleDict[GluGluToContinToZZTo2mu2tau.nickname][3]} 
GluGluToContinToZZTo2mu2tau.cuts = {"GluGluTo2mu2tau":""}
GluGluToContinToZZTo2mu2tau.file = GluGluToContinToZZTo2mu2tau.nickname+"_2016.root"
GluGluToContinToZZTo2mu2tau.classification = sampleDict[GluGluToContinToZZTo2mu2tau.nickname][3]
GluGluToContinToZZTo2mu2tau.color = ["#65E114"]    #do i need a list here??? to embed it
HAA_processes[GluGluToContinToZZTo2mu2tau.nickname]=GluGluToContinToZZTo2mu2tau

GluGluToContinToZZTo2mu2nu = Process()
GluGluToContinToZZTo2mu2nu.nickname = "GluGluTo2mu2nu"
GluGluToContinToZZTo2mu2nu.plotname = "GluGluTo2mu2nu"
GluGluToContinToZZTo2mu2nu.weights = {"xsec":sampleDict[GluGluToContinToZZTo2mu2nu.nickname][1],"nevents":sampleDict[GluGluToContinToZZTo2mu2nu.nickname][3]} 
GluGluToContinToZZTo2mu2nu.cuts = {"GluGluTo2mu2nu":""}
GluGluToContinToZZTo2mu2nu.file = GluGluToContinToZZTo2mu2nu.nickname+"_2016.root"
GluGluToContinToZZTo2mu2nu.classification = sampleDict[GluGluToContinToZZTo2mu2nu.nickname][3]
GluGluToContinToZZTo2mu2nu.color = ["#65E114"]    #do i need a list here??? to embed it
HAA_processes[GluGluToContinToZZTo2mu2nu.nickname]=GluGluToContinToZZTo2mu2nu

GluGluTo4mu = Process()
GluGluTo4mu.nickname = "GluGluTo4mu"
GluGluTo4mu.plotname = "GluGluTo4mu"
GluGluTo4mu.weights = {"xsec":sampleDict[GluGluTo4mu.nickname][1],"nevents":sampleDict[GluGluTo4mu.nickname][3]} 
GluGluTo4mu.cuts = {"GluGluTo4mu":""}
GluGluTo4mu.file = GluGluTo4mu.nickname+"_2016.root"
GluGluTo4mu.classification = sampleDict[GluGluTo4mu.nickname][3]
GluGluTo4mu.color = ["#65E114"]    #do i need a list here??? to embed it
HAA_processes[GluGluTo4mu.nickname]=GluGluTo4mu

GluGluToHToTauTau = Process()
GluGluToHToTauTau.nickname = "GluGluToHToTauTau"
GluGluToHToTauTau.plotname = "GluGluToHToTauTau"
GluGluToHToTauTau.weights = {"xsec":sampleDict[GluGluToHToTauTau.nickname][1],"nevents":sampleDict[GluGluToHToTauTau.nickname][3]} 
GluGluToHToTauTau.cuts = {"GluGluToHToTauTau":""}
GluGluToHToTauTau.file = GluGluToHToTauTau.nickname+"_2016.root"
GluGluToHToTauTau.classification = sampleDict[GluGluToHToTauTau.nickname][3]
GluGluToHToTauTau.color = ["#65E114"]    #do i need a list here??? to embed it
HAA_processes[GluGluToHToTauTau.nickname]=GluGluToHToTauTau

ttZ = Process()
ttZ.nickname = "ttZ"
ttZ.plotname = "ttZ"
ttZ.weights = {"xsec":sampleDict[ttZ.nickname][1],"nevents":sampleDict[ttZ.nickname][3]} 
ttZ.cuts = {"ttZ":""}
ttZ.file = ttZ.nickname+"_2016.root"
ttZ.classification = sampleDict[ttZ.nickname][3]
ttZ.color = ["#65E114"]    #do i need a list here??? to embed it
HAA_processes[ttZ.nickname]=ttZ

ttW = Process()
ttW.nickname = "ttW"
ttW.plotname = "ttW"
ttW.weights = {"xsec":sampleDict[ttW.nickname][1],"nevents":sampleDict[ttW.nickname][3]} 
ttW.cuts = {"ttW":""}
ttW.file = ttW.nickname+"_2016.root"
ttW.classification = sampleDict[ttW.nickname][3]
ttW.color = ["#65E114"]    #do i need a list here??? to embed it
HAA_processes[ttW.nickname]=ttW

WZZ = Process()
WZZ.nickname = "WZZ"
WZZ.plotname = "WZZ"
WZZ.weights = {"xsec":sampleDict[WZZ.nickname][1],"nevents":sampleDict[WZZ.nickname][3]} 
WZZ.cuts = {"WZZ":""}
WZZ.file = WZZ.nickname+"_2016.root"
WZZ.classification = sampleDict[WZZ.nickname][3]
WZZ.color = ["#65E114"]    #do i need a list here??? to embed it
HAA_processes[WZZ.nickname]=WZZ

WWZ = Process()
WWZ.nickname = "WWZ"
WWZ.plotname = "WWZ"
WWZ.weights = {"xsec":sampleDict[WWZ.nickname][1],"nevents":sampleDict[WWZ.nickname][3]} 
WWZ.cuts = {"WWZ":""}
WWZ.file = WWZ.nickname+"_2016.root"
WWZ.classification = sampleDict[WWZ.nickname][3]
WWZ.color = ["#65E114"]    #do i need a list here??? to embed it
HAA_processes[WWZ.nickname]=WWZ

HZJ = Process()
HZJ.nickname = "HZJ"
HZJ.plotname = "HZJ"
HZJ.weights = {"xsec":sampleDict[HZJ.nickname][1],"nevents":sampleDict[HZJ.nickname][3]} 
HZJ.cuts = {"HZJ":""}
HZJ.file = HZJ.nickname+"_2016.root"
HZJ.classification = sampleDict[HZJ.nickname][3]
HZJ.color = ["#65E114"]    #do i need a list here??? to embed it
HAA_processes[HZJ.nickname]=HZJ

ZZZ = Process()
ZZZ.nickname = "ZZZ"
ZZZ.plotname = "ZZZ"
ZZZ.weights = {"xsec":sampleDict[ZZZ.nickname][1],"nevents":sampleDict[ZZZ.nickname][3]} 
ZZZ.cuts = {"ZZZ":""}
ZZZ.file = ZZZ.nickname+"_2016.root"
ZZZ.classification = sampleDict[ZZZ.nickname][3]
ZZZ.color = ["#65E114"]    #do i need a list here??? to embed it
HAA_processes[ZZZ.nickname]=ZZZ

WWW_4F = Process()
WWW_4F.nickname = "WWW_4F"
WWW_4F.plotname = "WWW_4F"
WWW_4F.weights = {"xsec":sampleDict[WWW_4F.nickname][1],"nevents":sampleDict[WWW_4F.nickname][3]} 
WWW_4F.cuts = {"WWW_4F":""}
WWW_4F.file = WWW_4F.nickname+"_2016.root"
WWW_4F.classification = sampleDict[WWW_4F.nickname][3]
WWW_4F.color = ["#65E114"]    #do i need a list here??? to embed it
HAA_processes[WWW_4F.nickname]=WWW_4F

ZHToTauTau = Process()
ZHToTauTau.nickname = "ZHToTauTau"
ZHToTauTau.plotname = "ZHToTauTau"
ZHToTauTau.weights = {"xsec":sampleDict[ZHToTauTau.nickname][1],"nevents":sampleDict[ZHToTauTau.nickname][3]} 
ZHToTauTau.cuts = {"ZHToTauTau":""}
ZHToTauTau.file = ZHToTauTau.nickname+"_2016.root"
ZHToTauTau.classification = sampleDict[ZHToTauTau.nickname][3]
ZHToTauTau.color = ["#13E2FE"]    #do i need a list here??? to embed it
HAA_processes[ZHToTauTau.nickname]=ZHToTauTau

ZZ = Process()
ZZ.nickname = "ZZ"
ZZ.plotname = "ZZ"
ZZ.weights = {"xsec":sampleDict[ZZ.nickname][1],"nevents":sampleDict[ZZ.nickname][3]} 
ZZ.cuts = {"ZZ":""}
ZZ.file = ZZ.nickname+"_2016.root"
ZZ.classification = sampleDict[ZZ.nickname][3]
ZZ.color = ["#13E2FE"]    #do i need a list here??? to embed it
HAA_processes[ZZ.nickname]=ZZ

ZZTo2L2Q = Process()
ZZTo2L2Q.nickname = "ZZTo2L2Q"
ZZTo2L2Q.plotname = "ZZTo2L2Q"
ZZTo2L2Q.weights = {"xsec":sampleDict[ZZTo2L2Q.nickname][1],"nevents":sampleDict[ZZTo2L2Q.nickname][3]} 
ZZTo2L2Q.cuts = {"ZZTo2L2Q":""}
ZZTo2L2Q.file = ZZTo2L2Q.nickname+"_2016.root"
ZZTo2L2Q.classification = sampleDict[ZZTo2L2Q.nickname][3]
ZZTo2L2Q.color = ["#13E2FE"]    #do i need a list here??? to embed it
HAA_processes[ZZTo2L2Q.nickname]=ZZTo2L2Q

ZZTo4Lmcnloext1 = Process()
ZZTo4Lmcnloext1.nickname = "ZZTo4Lmcnloext1"
ZZTo4Lmcnloext1.plotname = "ZZTo4Lmcnloext1"
ZZTo4Lmcnloext1.weights = {"xsec":sampleDict[ZZTo4Lmcnloext1.nickname][1],"nevents":sampleDict[ZZTo4Lmcnloext1.nickname][3]} 
ZZTo4Lmcnloext1.cuts = {"ZZTo4Lmcnloext1":""}
ZZTo4Lmcnloext1.file = ZZTo4Lmcnloext1.nickname+"_2016.root"
ZZTo4Lmcnloext1.classification = sampleDict[ZZTo4Lmcnloext1.nickname][3]
ZZTo4Lmcnloext1.color = ["#13E2FE"]    #do i need a list here??? to embed it
HAA_processes[ZZTo4Lmcnloext1.nickname]=ZZTo4Lmcnloext1

ZZTo4L = Process()
ZZTo4L.nickname = "ZZTo4L"
ZZTo4L.plotname = "ZZTo4L"
ZZTo4L.weights = {"xsec":sampleDict[ZZTo4L.nickname][1],"nevents":sampleDict[ZZTo4L.nickname][3]} 
ZZTo4L.cuts = {"ZZTo4L":""}
ZZTo4L.file = ZZTo4L.nickname+"_2016.root"
ZZTo4L.classification = sampleDict[ZZTo4L.nickname][3]
ZZTo4L.color = ["#13E2FE"]    #do i need a list here??? to embed it
HAA_processes[ZZTo4L.nickname]=ZZTo4L

ZZTo4L_ext1 = Process()
ZZTo4L_ext1.nickname = "ZZTo4L_ext1"
ZZTo4L_ext1.plotname = "ZZTo4L_ext1"
ZZTo4L_ext1.weights = {"xsec":sampleDict[ZZTo4L_ext1.nickname][1],"nevents":sampleDict[ZZTo4L_ext1.nickname][3]} 
ZZTo4L_ext1.cuts = {"ZZTo4L_ext1":""}
ZZTo4L_ext1.file = ZZTo4L_ext1.nickname+"_2016.root"
ZZTo4L_ext1.classification = sampleDict[ZZTo4L_ext1.nickname][3]
ZZTo4L_ext1.color = ["#13E2FE"]    #do i need a list here??? to embed it
HAA_processes[ZZTo4L_ext1.nickname]=ZZTo4L_ext1

ZZext = Process()
ZZext.nickname = "ZZ_ext1"
ZZext.plotname = "ZZ"
ZZext.weights = {"xsec":sampleDict[ZZext.nickname][1],"nevents":sampleDict[ZZext.nickname][3]} 
ZZext.cuts = {"ZZ":""}
ZZext.file = ZZext.nickname+"_2016.root"
ZZext.classification = sampleDict[ZZext.nickname][3]
ZZext.color = ["#13E2FE"]    #do i need a list here??? to embed it
HAA_processes[ZZext.nickname]=ZZext

DY = Process()
DY.nickname = "DYJetsToLLext1"
DY.plotname = "ZTT"
DY.weights = {"xsec":sampleDict[DY.nickname][1],"nevents":sampleDict[DY.nickname][3]} #higgs xsec [pb] * 2hdm type Branching ratio
DY.cuts = {"Z":"","ZTT":[["gen_match_4","==",5]],"ZL":[["gen_match_4",">=",5]],"ZJ":[["gen_match_4",">",5]]} #higgs xsec [pb] * 2hdm type Branching ratio
DY.file = DY.nickname+"_2016.root"
DY.classification = sampleDict[DY.nickname][3]
DY.color = ["#CF8AC8"]    #do i need a list here??? to embed it
HAA_processes[DY.nickname]=DY

DYext = Process()
DYext.nickname = "DYJetsToLLext2"
DYext.plotname = "ZTT"
DYext.weights = {"xsec":sampleDict[DYext.nickname][1],"nevents":sampleDict[DYext.nickname][3]} #higgs xsec [pb] * 2hdm type Branching ratio
DYext.cuts = {"Z":"","ZTT":[["gen_match_4","==",5]],"ZL":[["gen_match_4",">=",5]],"ZJ":[["gen_match_4",">",5]]} #higgs xsec [pb] * 2hdm type Branching ratio
DYext.file = DYext.nickname+"_2016.root"
DYext.classification = sampleDict[DYext.nickname][3]
DYext.color = ["#CF8AC8"]    #do i need a list here??? to embed it
HAA_processes[DYext.nickname]=DYext

DY1 = Process()
DY1.nickname = "DY1JetsToLL"
DY1.plotname = "ZTT"
DY1.weights = {"xsec":sampleDict[DY1.nickname][1],"nevents":sampleDict[DY1.nickname][3]} #higgs xsec [pb] * 2hdm type Branching ratio
DY1.cuts = {"Z":"","ZTT":[["gen_match_4","==",5]],"ZL":[["gen_match_4",">=",5]],"ZJ":[["gen_match_4",">",5]]} #higgs xsec [pb] * 2hdm type Branching ratio
DY1.file = DY1.nickname+"_2016.root"
DY1.classification = sampleDict[DY1.nickname][3]
DY1.color = ["#CF8AC8"]    #do i need a list here??? to embed it
HAA_processes[DY1.nickname]=DY1

DY2 = Process()
DY2.nickname = "DY2JetsToLL"
DY2.plotname = "ZTT"
DY2.weights = {"xsec":sampleDict[DY2.nickname][1],"nevents":sampleDict[DY2.nickname][3]} #higgs xsec [pb] * 2hdm type Branching ratio
DY2.cuts = {"Z":"","ZTT":[["gen_match_4","==",5]],"ZL":[["gen_match_4",">=",5]],"ZJ":[["gen_match_4",">",5]]} #higgs xsec [pb] * 2hdm type Branching ratio
DY2.file = DY2.nickname+"_2016.root"
DY2.classification = sampleDict[DY2.nickname][3]
DY2.color = ["#CF8AC8"]    #do i need a list here??? to embed it
HAA_processes[DY2.nickname]=DY2

DY3 = Process()
DY3.nickname = "DY3JetsToLL"
DY3.plotname = "ZTT"
DY3.weights = {"xsec":sampleDict[DY3.nickname][1],"nevents":sampleDict[DY3.nickname][3]} #higgs xsec [pb] * 2hdm type Branching ratio
DY3.cuts = {"Z":"","ZTT":[["gen_match_4","==",5]],"ZL":[["gen_match_4",">=",5]],"ZJ":[["gen_match_4",">",5]]} #higgs xsec [pb] * 2hdm type Branching ratio
DY3.file = DY3.nickname+"_2016.root"
DY3.classification = sampleDict[DY3.nickname][3]
DY3.color = ["#CF8AC8"]    #do i need a list here??? to embed it
HAA_processes[DY3.nickname]=DY3

DY4 = Process()
DY4.nickname = "DY4JetsToLL"
DY4.plotname = "ZTT"
DY4.weights = {"xsec":sampleDict[DY4.nickname][1],"nevents":sampleDict[DY4.nickname][3]} #higgs xsec [pb] * 2hdm type Branching ratio
DY4.cuts = {"Z":"","ZTT":[["gen_match_4","==",5]],"ZL":[["gen_match_4",">=",5]],"ZJ":[["gen_match_4",">",5]]} #higgs xsec [pb] * 2hdm type Branching ratio
DY4.file = DY4.nickname+"_2016.root"
DY4.classification = sampleDict[DY4.nickname][3]
DY4.color = ["#CF8AC8"]    #do i need a list here??? to embed it
HAA_processes[DY4.nickname]=DY4

DYM = Process()
DYM.nickname = "DYJetsToLLM10to50"
DYM.plotname = "ZTT"
DYM.weights = {"xsec":sampleDict[DYM.nickname][1],"nevents":sampleDict[DYM.nickname][3]} #higgs xsec [pb] * 2hdm type Branching ratio
DYM.cuts = {"Z":"","ZTT":[["gen_match_4","==",5]],"ZL":[["gen_match_4",">=",5]],"ZJ":[["gen_match_4",">",5]]} #higgs xsec [pb] * 2hdm type Branching ratio
DYM.file = DYM.nickname+"_2016.root"
DYM.classification = sampleDict[DYM.nickname][3]
DYM.color = ["#CF8AC8"]    #do i need a list here??? to embed it
HAA_processes[DYM.nickname]=DYM

WJetsToLNu = Process()
WJetsToLNu.nickname = "WJetsToLNu"
WJetsToLNu.plotname = "W"
WJetsToLNu.weights = {"xsec":sampleDict[WJetsToLNu.nickname][1],"nevents":sampleDict[WJetsToLNu.nickname][3]}
WJetsToLNu.cuts = {"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]} 
WJetsToLNu.file = WJetsToLNu.nickname+"_2016.root"
WJetsToLNu.classification = sampleDict[WJetsToLNu.nickname][3]
WJetsToLNu.color = ["#CF8AC8"]    
HAA_processes[WJetsToLNu.nickname]=WJetsToLNu

WJetsToLNu_ext2 = Process()
WJetsToLNu_ext2.nickname = "WJetsToLNu_ext2"
WJetsToLNu_ext2.plotname = "W"
WJetsToLNu_ext2.weights = {"xsec":sampleDict[WJetsToLNu_ext2.nickname][1],"nevents":sampleDict[WJetsToLNu_ext2.nickname][3]}
WJetsToLNu_ext2.cuts = {"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]} 
WJetsToLNu_ext2.file = WJetsToLNu_ext2.nickname+"_2016.root"
WJetsToLNu_ext2.classification = sampleDict[WJetsToLNu_ext2.nickname][3]
WJetsToLNu_ext2.color = ["#CF8AC8"]    
HAA_processes[WJetsToLNu_ext2.nickname]=WJetsToLNu_ext2

WJetsext = Process()
WJetsext.nickname = "WJetsToLNuext"
WJetsext.plotname = "W"
WJetsext.weights = {"xsec":sampleDict[WJetsext.nickname][1],"nevents":sampleDict[WJetsext.nickname][3]}
WJetsext.cuts = {"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]} 
WJetsext.file = WJetsext.nickname+"_2016.root"
WJetsext.classification = sampleDict[WJetsext.nickname][3]
WJetsext.color = ["#CF8AC8"]    
HAA_processes[WJetsext.nickname]=WJetsext

W1Jets = Process()
W1Jets.nickname = "W1JetsToLNu"
W1Jets.plotname = "W"
W1Jets.weights = {"xsec":sampleDict[W1Jets.nickname][1],"nevents":sampleDict[W1Jets.nickname][3]}
W1Jets.cuts = {"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]} 
W1Jets.file = W1Jets.nickname+"_2016.root"
W1Jets.classification = sampleDict[W1Jets.nickname][3]
W1Jets.color = ["#CF8AC8"]    
HAA_processes[W1Jets.nickname]=W1Jets

W2Jets = Process()
W2Jets.nickname = "W2JetsToLNu"
W2Jets.plotname = "W"
W2Jets.weights = {"xsec":sampleDict[W2Jets.nickname][1],"nevents":sampleDict[W2Jets.nickname][3]}
W2Jets.cuts = {"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]} 
W2Jets.file = W2Jets.nickname+"_2016.root"
W2Jets.classification = sampleDict[W2Jets.nickname][3]
W2Jets.color = ["#CF8AC8"]    
HAA_processes[W2Jets.nickname]=W2Jets

W2Jetsext = Process()
W2Jetsext.nickname = "W2JetsToLNuext1"
W2Jetsext.plotname = "W"
W2Jetsext.weights = {"xsec":sampleDict[W2Jetsext.nickname][1],"nevents":sampleDict[W2Jetsext.nickname][3]}
W2Jetsext.cuts = {"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]} 
W2Jetsext.file = W2Jetsext.nickname+"_2016.root"
W2Jetsext.classification = sampleDict[W2Jetsext.nickname][3]
W2Jetsext.color = ["#CF8AC8"]    
HAA_processes[W2Jetsext.nickname]=W2Jetsext

W3Jets = Process()
W3Jets.nickname = "W3JetsToLNu"
W3Jets.plotname = "W"
W3Jets.weights = {"xsec":sampleDict[W3Jets.nickname][1],"nevents":sampleDict[W3Jets.nickname][3]}
W3Jets.cuts = {"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]} 
W3Jets.file = W3Jets.nickname+"_2016.root"
W3Jets.classification = sampleDict[W3Jets.nickname][3]
W3Jets.color = ["#CF8AC8"]    
HAA_processes[W3Jets.nickname]=W3Jets

W3Jetsext = Process()
W3Jetsext.nickname = "W3JetsToLNu_ext1"
W3Jetsext.plotname = "W"
W3Jetsext.weights = {"xsec":sampleDict[W3Jetsext.nickname][1],"nevents":sampleDict[W3Jetsext.nickname][3]}
W3Jetsext.cuts = {"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]} 
W3Jetsext.file = W3Jetsext.nickname+"_2016.root"
W3Jetsext.classification = sampleDict[W3Jetsext.nickname][3]
W3Jetsext.color = ["#CF8AC8"]    
HAA_processes[W3Jetsext.nickname]=W3Jetsext

W4Jets = Process()
W4Jets.nickname = "W4JetsToLNu"
W4Jets.plotname = "W"
W4Jets.weights = {"xsec":sampleDict[W4Jets.nickname][1],"nevents":sampleDict[W4Jets.nickname][3]}
W4Jets.cuts = {"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]} 
W4Jets.file = W4Jets.nickname+"_2016.root"
W4Jets.classification = sampleDict[W4Jets.nickname][3]
W4Jets.color = ["#CF8AC8"]    
HAA_processes[W4Jets.nickname]=W4Jets

W4Jetsext1 = Process()
W4Jetsext1.nickname = "W4JetsToLNu_ext1"
W4Jetsext1.plotname = "W"
W4Jetsext1.weights = {"xsec":sampleDict[W4Jetsext1.nickname][1],"nevents":sampleDict[W4Jetsext1.nickname][3]}
W4Jetsext1.cuts = {"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]} 
W4Jetsext1.file = W4Jetsext1.nickname+"_2016.root"
W4Jetsext1.classification = sampleDict[W4Jetsext1.nickname][3]
W4Jetsext1.color = ["#CF8AC8"]    
HAA_processes[W4Jetsext1.nickname]=W4Jetsext1

W4Jetsext2 = Process()
W4Jetsext2.nickname = "W4JetsToLNu_ext2"
W4Jetsext2.plotname = "W"
W4Jetsext2.weights = {"xsec":sampleDict[W4Jetsext2.nickname][1],"nevents":sampleDict[W4Jetsext2.nickname][3]}
W4Jetsext2.cuts = {"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]} 
W4Jetsext2.file = W4Jetsext2.nickname+"_2016.root"
W4Jetsext2.classification = sampleDict[W4Jetsext2.nickname][3]
W4Jetsext2.color = ["#CF8AC8"]    
HAA_processes[W4Jetsext2.nickname]=W4Jetsext2

W4JetsToLNuext2 = Process()
W4JetsToLNuext2.nickname = "W4JetsToLNuext2"
W4JetsToLNuext2.plotname = "W"
W4JetsToLNuext2.weights = {"xsec":sampleDict[W4JetsToLNuext2.nickname][1],"nevents":sampleDict[W4JetsToLNuext2.nickname][3]}
W4JetsToLNuext2.cuts = {"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]} 
W4JetsToLNuext2.file = W4JetsToLNuext2.nickname+"_2016.root"
W4JetsToLNuext2.classification = sampleDict[W4JetsToLNuext2.nickname][3]
W4JetsToLNuext2.color = ["#CF8AC8"]    
HAA_processes[W4JetsToLNuext2.nickname]=W4JetsToLNuext2

TTJets_1LT = Process()
TTJets_1LT.nickname = "TTJets_1LT"
TTJets_1LT.plotname = "TT"
TTJets_1LT.weights = {"xsec":sampleDict[TTJets_1LT.nickname][1],"nevents":sampleDict[TTJets_1LT.nickname][3]}
TTJets_1LT.cuts = {"TT":"","TTT":[["gen_match_4","==",5]],"TTL":[["gen_match_4",">=",5]],"TTJ":[["gen_match_4",">",5]]} 
TTJets_1LT.file = TTJets_1LT.nickname+"_2016.root"
TTJets_1LT.classification = sampleDict[TTJets_1LT.nickname][3]
TTJets_1LT.color = ["#CF8AC8"]    
HAA_processes[TTJets_1LT.nickname]=TTJets_1LT

TTJets_1LT_ext1 = Process()
TTJets_1LT_ext1.nickname = "TTJets_1LT_ext1"
TTJets_1LT_ext1.plotname = "TT"
TTJets_1LT_ext1.weights = {"xsec":sampleDict[TTJets_1LT_ext1.nickname][1],"nevents":sampleDict[TTJets_1LT_ext1.nickname][3]}
TTJets_1LT_ext1.cuts = {"TT":"","TTT":[["gen_match_4","==",5]],"TTL":[["gen_match_4",">=",5]],"TTJ":[["gen_match_4",">",5]]} 
TTJets_1LT_ext1.file = TTJets_1LT_ext1.nickname+"_2016.root"
TTJets_1LT_ext1.classification = sampleDict[TTJets_1LT_ext1.nickname][3]
TTJets_1LT_ext1.color = ["#CF8AC8"]    
HAA_processes[TTJets_1LT_ext1.nickname]=TTJets_1LT_ext1

TTJets_1LTbar = Process()
TTJets_1LTbar.nickname = "TTJets_1LTbar"
TTJets_1LTbar.plotname = "T"
TTJets_1LTbar.weights = {"xsec":sampleDict[TTJets_1LTbar.nickname][1],"nevents":sampleDict[TTJets_1LTbar.nickname][3]}
TTJets_1LTbar.cuts = {"TT":"","TTT":[["gen_match_4","==",5]],"TTL":[["gen_match_4",">=",5]],"TTJ":[["gen_match_4",">",5]]} 
TTJets_1LTbar.file = TTJets_1LTbar.nickname+"_2016.root"
TTJets_1LTbar.classification = sampleDict[TTJets_1LTbar.nickname][3]
TTJets_1LTbar.color = ["#CF8AC8"]    
HAA_processes[TTJets_1LTbar.nickname]=TTJets_1LTbar

TTJets_1LTbar_ext1 = Process()
TTJets_1LTbar_ext1.nickname = "TTJets_1LTbar_ext1"
TTJets_1LTbar_ext1.plotname = "T"
TTJets_1LTbar_ext1.weights = {"xsec":sampleDict[TTJets_1LTbar_ext1.nickname][1],"nevents":sampleDict[TTJets_1LTbar_ext1.nickname][3]}
TTJets_1LTbar_ext1.cuts = {"TT":"","TTT":[["gen_match_4","==",5]],"TTL":[["gen_match_4",">=",5]],"TTJ":[["gen_match_4",">",5]]} 
TTJets_1LTbar_ext1.file = TTJets_1LTbar_ext1.nickname+"_2016.root"
TTJets_1LTbar_ext1.classification = sampleDict[TTJets_1LTbar_ext1.nickname][3]
TTJets_1LTbar_ext1.color = ["#CF8AC8"]    
HAA_processes[TTJets_1LTbar_ext1.nickname]=TTJets_1LTbar_ext1

TTJets_DiLept = Process()
TTJets_DiLept.nickname = "TTJets_DiLept"
TTJets_DiLept.plotname = "T"
TTJets_DiLept.weights = {"xsec":sampleDict[TTJets_DiLept.nickname][1],"nevents":sampleDict[TTJets_DiLept.nickname][3]}
TTJets_DiLept.cuts = {"TT":"","TTT":[["gen_match_4","==",5]],"TTL":[["gen_match_4",">=",5]],"TTJ":[["gen_match_4",">",5]]} 
TTJets_DiLept.file = TTJets_DiLept.nickname+"_2016.root"
TTJets_DiLept.classification = sampleDict[TTJets_DiLept.nickname][3]
TTJets_DiLept.color = ["#CF8AC8"]    
HAA_processes[TTJets_DiLept.nickname]=TTJets_DiLept

TTJets_DiLept_ext1 = Process()
TTJets_DiLept_ext1.nickname = "TTJets_DiLept_ext1"
TTJets_DiLept_ext1.plotname = "T"
TTJets_DiLept_ext1.weights = {"xsec":sampleDict[TTJets_DiLept_ext1.nickname][1],"nevents":sampleDict[TTJets_DiLept_ext1.nickname][3]}
TTJets_DiLept_ext1.cuts = {"TT":"","TTT":[["gen_match_4","==",5]],"TTL":[["gen_match_4",">=",5]],"TTJ":[["gen_match_4",">",5]]} 
TTJets_DiLept_ext1.file = TTJets_DiLept_ext1.nickname+"_2016.root"
TTJets_DiLept_ext1.classification = sampleDict[TTJets_DiLept_ext1.nickname][3]
TTJets_DiLept_ext1.color = ["#CF8AC8"]    
HAA_processes[TTJets_DiLept_ext1.nickname]=TTJets_DiLept_ext1

ST_s = Process()
ST_s.nickname = "ST_s"
ST_s.plotname = "ST"
ST_s.weights = {"xsec":sampleDict[ST_s.nickname][1],"nevents":sampleDict[ST_s.nickname][3]}
ST_s.cuts = {"ST":"","STT":[["gen_match_4","==",5]],"STL":[["gen_match_4",">=",5]],"STJ":[["gen_match_4",">",5]]} 
ST_s.file = ST_s.nickname+"_2016.root"
ST_s.classification = sampleDict[ST_s.nickname][3]
ST_s.color = ["#CF8AC8"]    
HAA_processes[ST_s.nickname]=ST_s

ST_tW_antitop = Process()
ST_tW_antitop.nickname = "ST_tW_antitop_ext1"
ST_tW_antitop.plotname = "ST"
ST_tW_antitop.weights = {"xsec":sampleDict[ST_tW_antitop.nickname][1],"nevents":sampleDict[ST_tW_antitop.nickname][3]}
ST_tW_antitop.cuts = {"ST":"","STT":[["gen_match_4","==",5]],"STL":[["gen_match_4",">=",5]],"STJ":[["gen_match_4",">",5]]} 
ST_tW_antitop.file = ST_tW_antitop.nickname+"_2016.root"
ST_tW_antitop.classification = sampleDict[ST_tW_antitop.nickname][3]
ST_tW_antitop.color = ["#CF8AC8"]    
HAA_processes[ST_tW_antitop.nickname]=ST_tW_antitop

ST_tW_top = Process()
ST_tW_top.nickname = "ST_tW_top_ext1"
ST_tW_top.plotname = "ST"
ST_tW_top.weights = {"xsec":sampleDict[ST_tW_top.nickname][1],"nevents":sampleDict[ST_tW_top.nickname][3]}
ST_tW_top.cuts = {"ST":"","STT":[["gen_match_4","==",5]],"STL":[["gen_match_4",">=",5]],"STJ":[["gen_match_4",">",5]]} 
ST_tW_top.file = ST_tW_top.nickname+"_2016.root"
ST_tW_top.classification = sampleDict[ST_tW_top.nickname][3]
ST_tW_top.color = ["#CF8AC8"]    
HAA_processes[ST_tW_top.nickname]=ST_tW_top

ST_t_antitop = Process()
ST_t_antitop.nickname = "ST_t_antitop"
ST_t_antitop.plotname = "ST"
ST_t_antitop.weights = {"xsec":sampleDict[ST_t_antitop.nickname][1],"nevents":sampleDict[ST_t_antitop.nickname][3]}
ST_t_antitop.cuts = {"ST":"","STT":[["gen_match_4","==",5]],"STL":[["gen_match_4",">=",5]],"STJ":[["gen_match_4",">",5]]} 
ST_t_antitop.file = ST_t_antitop.nickname+"_2016.root"
ST_t_antitop.classification = sampleDict[ST_t_antitop.nickname][3]
ST_t_antitop.color = ["#CF8AC8"]    
HAA_processes[ST_t_antitop.nickname]=ST_t_antitop

ST_t_top = Process()
ST_t_top.nickname = "ST_t_top"
ST_t_top.plotname = "ST"
ST_t_top.weights = {"xsec":sampleDict[ST_t_top.nickname][1],"nevents":sampleDict[ST_t_top.nickname][3]}
ST_t_top.cuts = {"ST":"","STT":[["gen_match_4","==",5]],"STL":[["gen_match_4",">=",5]],"STJ":[["gen_match_4",">",5]]} 
ST_t_top.file = ST_t_top.nickname+"_2016.root"
ST_t_top.classification = sampleDict[ST_t_top.nickname][3]
ST_t_top.color = ["#CF8AC8"]    
HAA_processes[ST_t_top.nickname]=ST_t_top

EWKWMinus2Jets = Process()
EWKWMinus2Jets.nickname = "EWKWMinus2Jets"
EWKWMinus2Jets.plotname = "EWK"
EWKWMinus2Jets.weights = {"xsec":sampleDict[EWKWMinus2Jets.nickname][1],"nevents":sampleDict[EWKWMinus2Jets.nickname][3]}
EWKWMinus2Jets.cuts = {"EWK":""} 
EWKWMinus2Jets.file = EWKWMinus2Jets.nickname+"_2016.root"
EWKWMinus2Jets.classification = sampleDict[EWKWMinus2Jets.nickname][3]
EWKWMinus2Jets.color = ["#CF8AC8"]    
HAA_processes[EWKWMinus2Jets.nickname]=EWKWMinus2Jets

EWKWMinus2Jets_ext1 = Process()
EWKWMinus2Jets_ext1.nickname = "EWKWMinus2Jets_ext1"
EWKWMinus2Jets_ext1.plotname = "EWK"
EWKWMinus2Jets_ext1.weights = {"xsec":sampleDict[EWKWMinus2Jets_ext1.nickname][1],"nevents":sampleDict[EWKWMinus2Jets_ext1.nickname][3]}
EWKWMinus2Jets_ext1.cuts = {"EWK":""} 
EWKWMinus2Jets_ext1.file = EWKWMinus2Jets_ext1.nickname+"_2016.root"
EWKWMinus2Jets_ext1.classification = sampleDict[EWKWMinus2Jets_ext1.nickname][3]
EWKWMinus2Jets_ext1.color = ["#CF8AC8"]    
HAA_processes[EWKWMinus2Jets_ext1.nickname]=EWKWMinus2Jets_ext1

EWKWMinus2Jets_ext2 = Process()
EWKWMinus2Jets_ext2.nickname = "EWKWMinus2Jets_ext2"
EWKWMinus2Jets_ext2.plotname = "EWK"
EWKWMinus2Jets_ext2.weights = {"xsec":sampleDict[EWKWMinus2Jets_ext2.nickname][1],"nevents":sampleDict[EWKWMinus2Jets_ext2.nickname][3]}
EWKWMinus2Jets_ext2.cuts = {"EWK":""} 
EWKWMinus2Jets_ext2.file = EWKWMinus2Jets_ext2.nickname+"_2016.root"
EWKWMinus2Jets_ext2.classification = sampleDict[EWKWMinus2Jets_ext2.nickname][3]
EWKWMinus2Jets_ext2.color = ["#CF8AC8"]    
HAA_processes[EWKWMinus2Jets_ext2.nickname]=EWKWMinus2Jets_ext2

EWKWPlus2Jets = Process()
EWKWPlus2Jets.nickname = "EWKWPlus2Jets"
EWKWPlus2Jets.plotname = "EWK"
EWKWPlus2Jets.weights = {"xsec":sampleDict[EWKWPlus2Jets.nickname][1],"nevents":sampleDict[EWKWPlus2Jets.nickname][3]}
EWKWPlus2Jets.cuts = {"EWK":""} 
EWKWPlus2Jets.file = EWKWPlus2Jets.nickname+"_2016.root"
EWKWPlus2Jets.classification = sampleDict[EWKWPlus2Jets.nickname][3]
EWKWPlus2Jets.color = ["#CF8AC8"]    
HAA_processes[EWKWPlus2Jets.nickname]=EWKWPlus2Jets

EWKWPlus2Jets_ext1 = Process()
EWKWPlus2Jets_ext1.nickname = "EWKWPlus2Jets_ext1"
EWKWPlus2Jets_ext1.plotname = "EWK"
EWKWPlus2Jets_ext1.weights = {"xsec":sampleDict[EWKWPlus2Jets_ext1.nickname][1],"nevents":sampleDict[EWKWPlus2Jets_ext1.nickname][3]}
EWKWPlus2Jets_ext1.cuts = {"EWK":""} 
EWKWPlus2Jets_ext1.file = EWKWPlus2Jets_ext1.nickname+"_2016.root"
EWKWPlus2Jets_ext1.classification = sampleDict[EWKWPlus2Jets_ext1.nickname][3]
EWKWPlus2Jets_ext1.color = ["#CF8AC8"]    
HAA_processes[EWKWPlus2Jets_ext1.nickname]=EWKWPlus2Jets_ext1

EWKWPlus2Jets_ext2 = Process()
EWKWPlus2Jets_ext2.nickname = "EWKWPlus2Jets_ext2"
EWKWPlus2Jets_ext2.plotname = "EWK"
EWKWPlus2Jets_ext2.weights = {"xsec":sampleDict[EWKWPlus2Jets_ext2.nickname][1],"nevents":sampleDict[EWKWPlus2Jets_ext2.nickname][3]}
EWKWPlus2Jets_ext2.cuts = {"EWK":""} 
EWKWPlus2Jets_ext2.file = EWKWPlus2Jets_ext2.nickname+"_2016.root"
EWKWPlus2Jets_ext2.classification = sampleDict[EWKWPlus2Jets_ext2.nickname][3]
EWKWPlus2Jets_ext2.color = ["#CF8AC8"]    
HAA_processes[EWKWPlus2Jets_ext2.nickname]=EWKWPlus2Jets_ext2

EWKZ2JetsLL = Process()
EWKZ2JetsLL.nickname = "EWKZ2JetsLL"
EWKZ2JetsLL.plotname = "EWK"
EWKZ2JetsLL.weights = {"xsec":sampleDict[EWKZ2JetsLL.nickname][1],"nevents":sampleDict[EWKZ2JetsLL.nickname][3]}
EWKZ2JetsLL.cuts = {"EWK":""} 
EWKZ2JetsLL.file = EWKZ2JetsLL.nickname+"_2016.root"
EWKZ2JetsLL.classification = sampleDict[EWKZ2JetsLL.nickname][3]
EWKZ2JetsLL.color = ["#CF8AC8"]    
HAA_processes[EWKZ2JetsLL.nickname]=EWKZ2JetsLL

EWKZ2JetsLL_ext1 = Process()
EWKZ2JetsLL_ext1.nickname = "EWKZ2JetsLL_ext1"
EWKZ2JetsLL_ext1.plotname = "EWK"
EWKZ2JetsLL_ext1.weights = {"xsec":sampleDict[EWKZ2JetsLL_ext1.nickname][1],"nevents":sampleDict[EWKZ2JetsLL_ext1.nickname][3]}
EWKZ2JetsLL_ext1.cuts = {"EWK":""} 
EWKZ2JetsLL_ext1.file = EWKZ2JetsLL_ext1.nickname+"_2016.root"
EWKZ2JetsLL_ext1.classification = sampleDict[EWKZ2JetsLL_ext1.nickname][3]
EWKZ2JetsLL_ext1.color = ["#CF8AC8"]    
HAA_processes[EWKZ2JetsLL_ext1.nickname]=EWKZ2JetsLL_ext1

EWKZ2JetsLL_ext2 = Process()
EWKZ2JetsLL_ext2.nickname = "EWKZ2JetsLL_ext2"
EWKZ2JetsLL_ext2.plotname = "EWK"
EWKZ2JetsLL_ext2.weights = {"xsec":sampleDict[EWKZ2JetsLL_ext2.nickname][1],"nevents":sampleDict[EWKZ2JetsLL_ext2.nickname][3]}
EWKZ2JetsLL_ext2.cuts = {"EWK":""} 
EWKZ2JetsLL_ext2.file = EWKZ2JetsLL_ext2.nickname+"_2016.root"
EWKZ2JetsLL_ext2.classification = sampleDict[EWKZ2JetsLL_ext2.nickname][3]
EWKZ2JetsLL_ext2.color = ["#CF8AC8"]    
HAA_processes[EWKZ2JetsLL_ext2.nickname]=EWKZ2JetsLL_ext2

EWKZ2JetsNuNu = Process()
EWKZ2JetsNuNu.nickname = "EWKZ2JetsNuNu"
EWKZ2JetsNuNu.plotname = "EWK"
EWKZ2JetsNuNu.weights = {"xsec":sampleDict[EWKZ2JetsNuNu.nickname][1],"nevents":sampleDict[EWKZ2JetsNuNu.nickname][3]}
EWKZ2JetsNuNu.cuts = {"EWK":""} 
EWKZ2JetsNuNu.file = EWKZ2JetsNuNu.nickname+"_2016.root"
EWKZ2JetsNuNu.classification = sampleDict[EWKZ2JetsNuNu.nickname][3]
EWKZ2JetsNuNu.color = ["#CF8AC8"]    
HAA_processes[EWKZ2JetsNuNu.nickname]=EWKZ2JetsNuNu

EWKZ2JetsNuNu_ext1 = Process()
EWKZ2JetsNuNu_ext1.nickname = "EWKZ2JetsNuNu_ext1"
EWKZ2JetsNuNu_ext1.plotname = "EWK"
EWKZ2JetsNuNu_ext1.weights = {"xsec":sampleDict[EWKZ2JetsNuNu_ext1.nickname][1],"nevents":sampleDict[EWKZ2JetsNuNu_ext1.nickname][3]}
EWKZ2JetsNuNu_ext1.cuts = {"EWK":""} 
EWKZ2JetsNuNu_ext1.file = EWKZ2JetsNuNu_ext1.nickname+"_2016.root"
EWKZ2JetsNuNu_ext1.classification = sampleDict[EWKZ2JetsNuNu_ext1.nickname][3]
EWKZ2JetsNuNu_ext1.color = ["#CF8AC8"]    
HAA_processes[EWKZ2JetsNuNu_ext1.nickname]=EWKZ2JetsNuNu_ext1

EWKZ2JetsNuNu_ext2 = Process()
EWKZ2JetsNuNu_ext2.nickname = "EWKZ2JetsNuNu_ext2"
EWKZ2JetsNuNu_ext2.plotname = "EWK"
EWKZ2JetsNuNu_ext2.weights = {"xsec":sampleDict[EWKZ2JetsNuNu_ext2.nickname][1],"nevents":sampleDict[EWKZ2JetsNuNu_ext2.nickname][3]}
EWKZ2JetsNuNu_ext2.cuts = {"EWK":""} 
EWKZ2JetsNuNu_ext2.file = EWKZ2JetsNuNu_ext2.nickname+"_2016.root"
EWKZ2JetsNuNu_ext2.classification = sampleDict[EWKZ2JetsNuNu_ext2.nickname][3]
EWKZ2JetsNuNu_ext2.color = ["#CF8AC8"]    
HAA_processes[EWKZ2JetsNuNu_ext2.nickname]=EWKZ2JetsNuNu_ext2

TTGamma_1LT = Process()
TTGamma_1LT.nickname = "TTGamma_1LT"
TTGamma_1LT.plotname = "TT"
TTGamma_1LT.weights = {"xsec":sampleDict[TTGamma_1LT.nickname][1],"nevents":sampleDict[TTGamma_1LT.nickname][3]}
TTGamma_1LT.cuts = {"rare":""} 
TTGamma_1LT.file = TTGamma_1LT.nickname+"_2016.root"
TTGamma_1LT.classification = sampleDict[TTGamma_1LT.nickname][3]
TTGamma_1LT.color = ["#CF8AC8"]    
HAA_processes[TTGamma_1LT.nickname]=TTGamma_1LT

TTGamma_Hadr = Process()
TTGamma_Hadr.nickname = "TTGamma_Hadr"
TTGamma_Hadr.plotname = "TT"
TTGamma_Hadr.weights = {"xsec":sampleDict[TTGamma_Hadr.nickname][1],"nevents":sampleDict[TTGamma_Hadr.nickname][3]}
TTGamma_Hadr.cuts = {"rare":""} 
TTGamma_Hadr.file = TTGamma_Hadr.nickname+"_2016.root"
TTGamma_Hadr.classification = sampleDict[TTGamma_Hadr.nickname][3]
TTGamma_Hadr.color = ["#CF8AC8"]    
HAA_processes[TTGamma_Hadr.nickname]=TTGamma_Hadr

TTGamma_1LTbar = Process()
TTGamma_1LTbar.nickname = "TTGamma_1LTbar"
TTGamma_1LTbar.plotname = "TT"
TTGamma_1LTbar.weights = {"xsec":sampleDict[TTGamma_1LTbar.nickname][1],"nevents":sampleDict[TTGamma_1LTbar.nickname][3]}
TTGamma_1LTbar.cuts = {"rare":""} 
TTGamma_1LTbar.file = TTGamma_1LTbar.nickname+"_2016.root"
TTGamma_1LTbar.classification = sampleDict[TTGamma_1LTbar.nickname][3]
TTGamma_1LTbar.color = ["#CF8AC8"]    
HAA_processes[TTGamma_1LTbar.nickname]=TTGamma_1LTbar

TTGamma_2L = Process()
TTGamma_2L.nickname = "TTGamma_2L"
TTGamma_2L.plotname = "TT"
TTGamma_2L.weights = {"xsec":sampleDict[TTGamma_2L.nickname][1],"nevents":sampleDict[TTGamma_2L.nickname][3]}
TTGamma_2L.cuts = {"rare":""} 
TTGamma_2L.file = TTGamma_2L.nickname+"_2016.root"
TTGamma_2L.classification = sampleDict[TTGamma_2L.nickname][3]
TTGamma_2L.color = ["#CF8AC8"]    
HAA_processes[TTGamma_2L.nickname]=TTGamma_2L

TTHH = Process()
TTHH.nickname = "TTHH_ext1"
TTHH.plotname = "TT"
TTHH.weights = {"xsec":sampleDict[TTHH.nickname][1],"nevents":sampleDict[TTHH.nickname][3]}
TTHH.cuts = {"rare":""} 
TTHH.file = TTHH.nickname+"_2016.root"
TTHH.classification = sampleDict[TTHH.nickname][3]
TTHH.color = ["#CF8AC8"]    
HAA_processes[TTHH.nickname]=TTHH

TTTT = Process()
TTTT.nickname = "TTTT"
TTTT.plotname = "TT"
TTTT.weights = {"xsec":sampleDict[TTTT.nickname][1],"nevents":sampleDict[TTTT.nickname][3]}
TTTT.cuts = {"rare":""} 
TTTT.file = TTTT.nickname+"_2016.root"
TTTT.classification = sampleDict[TTTT.nickname][3]
TTTT.color = ["#CF8AC8"]    
HAA_processes[TTTT.nickname]=TTTT

TTTJ_ext1 = Process()
TTTJ_ext1.nickname = "TTTJ_ext1"
TTTJ_ext1.plotname = "TT"
TTTJ_ext1.weights = {"xsec":sampleDict[TTTJ_ext1.nickname][1],"nevents":sampleDict[TTTJ_ext1.nickname][3]}
TTTJ_ext1.cuts = {"rare":""} 
TTTJ_ext1.file = TTTJ_ext1.nickname+"_2016.root"
TTTJ_ext1.classification = sampleDict[TTTJ_ext1.nickname][3]
TTTJ_ext1.color = ["#CF8AC8"]    
HAA_processes[TTTJ_ext1.nickname]=TTTJ_ext1

TTWH_ext1 = Process()
TTWH_ext1.nickname = "TTWH_ext1"
TTWH_ext1.plotname = "TT"
TTWH_ext1.weights = {"xsec":sampleDict[TTWH_ext1.nickname][1],"nevents":sampleDict[TTWH_ext1.nickname][3]}
TTWH_ext1.cuts = {"rare":""} 
TTWH_ext1.file = TTWH_ext1.nickname+"_2016.root"
TTWH_ext1.classification = sampleDict[TTWH_ext1.nickname][3]
TTWH_ext1.color = ["#CF8AC8"]    
HAA_processes[TTWH_ext1.nickname]=TTWH_ext1

TTTW = Process()
TTTW.nickname = "TTTW_ext1"
TTTW.plotname = "TT"
TTTW.weights = {"xsec":sampleDict[TTTW.nickname][1],"nevents":sampleDict[TTTW.nickname][3]}
TTTW.cuts = {"rare":""} 
TTTW.file = TTTW.nickname+"_2016.root"
TTTW.classification = sampleDict[TTTW.nickname][3]
TTTW.color = ["#CF8AC8"]    
HAA_processes[TTTW.nickname]=TTTW

TTW = Process()
TTW.nickname = "TTWJetsToLNu_ext1"
TTW.plotname = "TT"
TTW.weights = {"xsec":sampleDict[TTW.nickname][1],"nevents":sampleDict[TTW.nickname][3]}
TTW.cuts = {"rare":""} 
TTW.file = TTW.nickname+"_2016.root"
TTW.classification = sampleDict[TTW.nickname][3]
TTW.color = ["#CF8AC8"]    
HAA_processes[TTW.nickname]=TTW

TTWext2 = Process()
TTWext2.nickname = "TTWJetsToLNu_ext2"
TTWext2.plotname = "TT"
TTWext2.weights = {"xsec":sampleDict[TTWext2.nickname][1],"nevents":sampleDict[TTWext2.nickname][3]}
TTWext2.cuts = {"rare":""} 
TTWext2.file = TTWext2.nickname+"_2016.root"
TTWext2.classification = sampleDict[TTWext2.nickname][3]
TTWext2.color = ["#CF8AC8"]    
HAA_processes[TTWext2.nickname]=TTWext2

TTWJetsToQQ = Process()
TTWJetsToQQ.nickname = "TTWJetsToQQ"
TTWJetsToQQ.plotname = "TT"
TTWJetsToQQ.weights = {"xsec":sampleDict[TTWJetsToQQ.nickname][1],"nevents":sampleDict[TTWJetsToQQ.nickname][3]}
TTWJetsToQQ.cuts = {"rare":""} 
TTWJetsToQQ.file = TTWJetsToQQ.nickname+"_2016.root"
TTWJetsToQQ.classification = sampleDict[TTWJetsToQQ.nickname][3]
TTWJetsToQQ.color = ["#CF8AC8"]    
HAA_processes[TTWJetsToQQ.nickname]=TTWJetsToQQ

TTWW = Process()
TTWW.nickname = "TTWW_ext1"
TTWW.plotname = "TT"
TTWW.weights = {"xsec":sampleDict[TTWW.nickname][1],"nevents":sampleDict[TTWW.nickname][3]}
TTWW.cuts = {"rare":""} 
TTWW.file = TTWW.nickname+"_2016.root"
TTWW.classification = sampleDict[TTWW.nickname][3]
TTWW.color = ["#CF8AC8"]    
HAA_processes[TTWW.nickname]=TTWW

TTWZ = Process()
TTWZ.nickname = "TTWZ_ext1"
TTWZ.plotname = "TT"
TTWZ.weights = {"xsec":sampleDict[TTWZ.nickname][1],"nevents":sampleDict[TTWZ.nickname][3]}
TTWZ.cuts = {"rare":""} 
TTWZ.file = TTWZ.nickname+"_2016.root"
TTWZ.classification = sampleDict[TTWZ.nickname][3]
TTWZ.color = ["#CF8AC8"]    
HAA_processes[TTWZ.nickname]=TTWZ

TTZH = Process()
TTZH.nickname = "TTZH_ext1"
TTZH.plotname = "TT"
TTZH.weights = {"xsec":sampleDict[TTZH.nickname][1],"nevents":sampleDict[TTZH.nickname][3]}
TTZH.cuts = {"rare":""} 
TTZH.file = TTZH.nickname+"_2016.root"
TTZH.classification = sampleDict[TTZH.nickname][3]
TTZH.color = ["#CF8AC8"]    
HAA_processes[TTZH.nickname]=TTZH

TTZToQQ = Process()
TTZToQQ.nickname = "TTZToQQ"
TTZToQQ.plotname = "TT"
TTZToQQ.weights = {"xsec":sampleDict[TTZToQQ.nickname][1],"nevents":sampleDict[TTZToQQ.nickname][3]}
TTZToQQ.cuts = {"rare":""} 
TTZToQQ.file = TTZToQQ.nickname+"_2016.root"
TTZToQQ.classification = sampleDict[TTZToQQ.nickname][3]
TTZToQQ.color = ["#CF8AC8"]    
HAA_processes[TTZToQQ.nickname]=TTZToQQ

TTZToLLNuNu_ext1 = Process()
TTZToLLNuNu_ext1.nickname = "TTZToLLNuNu_ext1"
TTZToLLNuNu_ext1.plotname = "TT"
TTZToLLNuNu_ext1.weights = {"xsec":sampleDict[TTZToLLNuNu_ext1.nickname][1],"nevents":sampleDict[TTZToLLNuNu_ext1.nickname][3]}
TTZToLLNuNu_ext1.cuts = {"rare":""} 
TTZToLLNuNu_ext1.file = TTZToLLNuNu_ext1.nickname+"_2016.root"
TTZToLLNuNu_ext1.classification = sampleDict[TTZToLLNuNu_ext1.nickname][3]
TTZToLLNuNu_ext1.color = ["#CF8AC8"]    
HAA_processes[TTZToLLNuNu_ext1.nickname]=TTZToLLNuNu_ext1

TTZToLLNuNu_ext2 = Process()
TTZToLLNuNu_ext2.nickname = "TTZToLLNuNu_ext2"
TTZToLLNuNu_ext2.plotname = "TT"
TTZToLLNuNu_ext2.weights = {"xsec":sampleDict[TTZToLLNuNu_ext2.nickname][1],"nevents":sampleDict[TTZToLLNuNu_ext2.nickname][3]}
TTZToLLNuNu_ext2.cuts = {"rare":""} 
TTZToLLNuNu_ext2.file = TTZToLLNuNu_ext2.nickname+"_2016.root"
TTZToLLNuNu_ext2.classification = sampleDict[TTZToLLNuNu_ext2.nickname][3]
TTZToLLNuNu_ext2.color = ["#CF8AC8"]    
HAA_processes[TTZToLLNuNu_ext2.nickname]=TTZToLLNuNu_ext2

TTZToLLNuNu_ext3 = Process()
TTZToLLNuNu_ext3.nickname = "TTZToLLNuNu_ext3"
TTZToLLNuNu_ext3.plotname = "TT"
TTZToLLNuNu_ext3.weights = {"xsec":sampleDict[TTZToLLNuNu_ext3.nickname][1],"nevents":sampleDict[TTZToLLNuNu_ext3.nickname][3]}
TTZToLLNuNu_ext3.cuts = {"rare":""} 
TTZToLLNuNu_ext3.file = TTZToLLNuNu_ext3.nickname+"_2016.root"
TTZToLLNuNu_ext3.classification = sampleDict[TTZToLLNuNu_ext3.nickname][3]
TTZToLLNuNu_ext3.color = ["#CF8AC8"]    
HAA_processes[TTZToLLNuNu_ext3.nickname]=TTZToLLNuNu_ext3

TTZZ = Process()
TTZZ.nickname = "TTZZ_ext1"
TTZZ.plotname = "TT"
TTZZ.weights = {"xsec":sampleDict[TTZZ.nickname][1],"nevents":sampleDict[TTZZ.nickname][3]}
TTZZ.cuts = {"rare":""} 
TTZZ.file = TTZZ.nickname+"_2016.root"
TTZZ.classification = sampleDict[TTZZ.nickname][3]
TTZZ.color = ["#CF8AC8"]    
HAA_processes[TTZZ.nickname]=TTZZ

WWTo1L1Nu2Q = Process()
WWTo1L1Nu2Q.nickname = "WWTo1L1Nu2Q"
WWTo1L1Nu2Q.plotname = "W"
WWTo1L1Nu2Q.weights = {"xsec":sampleDict[WWTo1L1Nu2Q.nickname][1],"nevents":sampleDict[WWTo1L1Nu2Q.nickname][3]}
WWTo1L1Nu2Q.cuts = {"rare":""} 
WWTo1L1Nu2Q.file = WWTo1L1Nu2Q.nickname+"_2016.root"
WWTo1L1Nu2Q.classification = sampleDict[WWTo1L1Nu2Q.nickname][3]
WWTo1L1Nu2Q.color = ["#CF8AC8"]    
HAA_processes[WWTo1L1Nu2Q.nickname]=WWTo1L1Nu2Q

WZTo1L1Nu2Q = Process()
WZTo1L1Nu2Q.nickname = "WZTo1L1Nu2Q"
WZTo1L1Nu2Q.plotname = "W"
WZTo1L1Nu2Q.weights = {"xsec":sampleDict[WZTo1L1Nu2Q.nickname][1],"nevents":sampleDict[WZTo1L1Nu2Q.nickname][3]}
WZTo1L1Nu2Q.cuts = {"rare":""} 
WZTo1L1Nu2Q.file = WZTo1L1Nu2Q.nickname+"_2016.root"
WZTo1L1Nu2Q.classification = sampleDict[WZTo1L1Nu2Q.nickname][3]
WZTo1L1Nu2Q.color = ["#CF8AC8"]    
HAA_processes[WZTo1L1Nu2Q.nickname]=WZTo1L1Nu2Q

WZG = Process()
WZG.nickname = "WZG"
WZG.plotname = "W"
WZG.weights = {"xsec":sampleDict[WZG.nickname][1],"nevents":sampleDict[WZG.nickname][3]}
WZG.cuts = {"W":""} 
WZG.file = WZG.nickname+"_2016.root"
WZG.classification = sampleDict[WZG.nickname][3]
WZG.color = ["#CF8AC8"]    
HAA_processes[WZG.nickname]=WZG

WZTo1L3Nu = Process()
WZTo1L3Nu.nickname = "WZTo1L3Nu"
WZTo1L3Nu.plotname = "W"
WZTo1L3Nu.weights = {"xsec":sampleDict[WZTo1L3Nu.nickname][1],"nevents":sampleDict[WZTo1L3Nu.nickname][3]}
WZTo1L3Nu.cuts = {"rare":""} 
WZTo1L3Nu.file = WZTo1L3Nu.nickname+"_2016.root"
WZTo1L3Nu.classification = sampleDict[WZTo1L3Nu.nickname][3]
WZTo1L3Nu.color = ["#CF8AC8"]    
HAA_processes[WZTo1L3Nu.nickname]=WZTo1L3Nu

WZTo2L2Q = Process()
WZTo2L2Q.nickname = "WZTo2L2Q"
WZTo2L2Q.plotname = "W"
WZTo2L2Q.weights = {"xsec":sampleDict[WZTo2L2Q.nickname][1],"nevents":sampleDict[WZTo2L2Q.nickname][3]}
WZTo2L2Q.cuts = {"rare":""} 
WZTo2L2Q.file = WZTo2L2Q.nickname+"_2016.root"
WZTo2L2Q.classification = sampleDict[WZTo2L2Q.nickname][3]
WZTo2L2Q.color = ["#CF8AC8"]    
HAA_processes[WZTo2L2Q.nickname]=WZTo2L2Q

WZTo2Q2Nu = Process()
WZTo2Q2Nu.nickname = "WZTo2Q2Nu"
WZTo2Q2Nu.plotname = "W"
WZTo2Q2Nu.weights = {"xsec":sampleDict[WZTo2Q2Nu.nickname][1],"nevents":sampleDict[WZTo2Q2Nu.nickname][3]}
WZTo2Q2Nu.cuts = {"rare":""} 
WZTo2Q2Nu.file = WZTo2Q2Nu.nickname+"_2016.root"
WZTo2Q2Nu.classification = sampleDict[WZTo2Q2Nu.nickname][3]
WZTo2Q2Nu.color = ["#CF8AC8"]    
HAA_processes[WZTo2Q2Nu.nickname]=WZTo2Q2Nu


WZTo3LNumll_pow = Process()
WZTo3LNumll_pow.nickname = "WZTo3LNumll_pow"
WZTo3LNumll_pow.plotname = "W"
WZTo3LNumll_pow.weights = {"xsec":sampleDict[WZTo3LNumll_pow.nickname][1],"nevents":sampleDict[WZTo3LNumll_pow.nickname][3]}
WZTo3LNumll_pow.cuts = {"rare":""} 
WZTo3LNumll_pow.file = WZTo3LNumll_pow.nickname+"_2016.root"
WZTo3LNumll_pow.classification = sampleDict[WZTo3LNumll_pow.nickname][3]
WZTo3LNumll_pow.color = ["#CF8AC8"]    
HAA_processes[WZTo3LNumll_pow.nickname]=WZTo3LNumll_pow

WZTo3LNumll_pow_ext1 = Process()
WZTo3LNumll_pow_ext1.nickname = "WZTo3LNumll_pow_ext1"
WZTo3LNumll_pow_ext1.plotname = "W"
WZTo3LNumll_pow_ext1.weights = {"xsec":sampleDict[WZTo3LNumll_pow_ext1.nickname][1],"nevents":sampleDict[WZTo3LNumll_pow_ext1.nickname][3]}
WZTo3LNumll_pow_ext1.cuts = {"rare":""} 
WZTo3LNumll_pow_ext1.file = WZTo3LNumll_pow_ext1.nickname+"_2016.root"
WZTo3LNumll_pow_ext1.classification = sampleDict[WZTo3LNumll_pow_ext1.nickname][3]
WZTo3LNumll_pow_ext1.color = ["#CF8AC8"]    
HAA_processes[WZTo3LNumll_pow_ext1.nickname]=WZTo3LNumll_pow_ext1


WZTo1L1Nu2Q = Process()
WZTo1L1Nu2Q.nickname = "WZTo1L1Nu2Q"
WZTo1L1Nu2Q.plotname = "W"
WZTo1L1Nu2Q.weights = {"xsec":sampleDict[WZTo1L1Nu2Q.nickname][1],"nevents":sampleDict[WZTo1L1Nu2Q.nickname][3]}
WZTo1L1Nu2Q.cuts = {"rare":""} 
WZTo1L1Nu2Q.file = WZTo1L1Nu2Q.nickname+"_2016.root"
WZTo1L1Nu2Q.classification = sampleDict[WZTo1L1Nu2Q.nickname][3]
WZTo1L1Nu2Q.color = ["#CF8AC8"]    
HAA_processes[WZTo1L1Nu2Q.nickname]=WZTo1L1Nu2Q

WWTo2L2Nu = Process()
WWTo2L2Nu.nickname = "WWTo2L2Nu"
WWTo2L2Nu.plotname = "W"
WWTo2L2Nu.weights = {"xsec":sampleDict[WWTo2L2Nu.nickname][1],"nevents":sampleDict[WWTo2L2Nu.nickname][3]}
WWTo2L2Nu.cuts = {"rare":""} 
WWTo2L2Nu.file = WWTo2L2Nu.nickname+"_2016.root"
WWTo2L2Nu.classification = sampleDict[WWTo2L2Nu.nickname][3]
WWTo2L2Nu.color = ["#CF8AC8"]    
HAA_processes[WWTo2L2Nu.nickname]=WWTo2L2Nu

WWTo4Q = Process()
WWTo4Q.nickname = "WWTo4Q"
WWTo4Q.plotname = "W"
WWTo4Q.weights = {"xsec":sampleDict[WWTo4Q.nickname][1],"nevents":sampleDict[WWTo4Q.nickname][3]}
WWTo4Q.cuts = {"rare":""} 
WWTo4Q.file = WWTo4Q.nickname+"_2016.root"
WWTo4Q.classification = sampleDict[WWTo4Q.nickname][3]
WWTo4Q.color = ["#CF8AC8"]    
HAA_processes[WWTo4Q.nickname]=WWTo4Q

WWToLNuQQ = Process()
WWToLNuQQ.nickname = "WWToLNuQQ"
WWToLNuQQ.plotname = "W"
WWToLNuQQ.weights = {"xsec":sampleDict[WWToLNuQQ.nickname][1],"nevents":sampleDict[WWToLNuQQ.nickname][3]}
WWToLNuQQ.cuts = {"rare":""} 
WWToLNuQQ.file = WWToLNuQQ.nickname+"_2016.root"
WWToLNuQQ.classification = sampleDict[WWToLNuQQ.nickname][3]
WWToLNuQQ.color = ["#CF8AC8"]    
HAA_processes[WWToLNuQQ.nickname]=WWToLNuQQ

WWToLNuQQ_ext1 = Process()
WWToLNuQQ_ext1.nickname = "WWToLNuQQ_ext1"
WWToLNuQQ_ext1.plotname = "W"
WWToLNuQQ_ext1.weights = {"xsec":sampleDict[WWToLNuQQ_ext1.nickname][1],"nevents":sampleDict[WWToLNuQQ_ext1.nickname][3]}
WWToLNuQQ_ext1.cuts = {"rare":""} 
WWToLNuQQ_ext1.file = WWToLNuQQ_ext1.nickname+"_2016.root"
WWToLNuQQ_ext1.classification = sampleDict[WWToLNuQQ_ext1.nickname][3]
WWToLNuQQ_ext1.color = ["#CF8AC8"]    
HAA_processes[WWToLNuQQ_ext1.nickname]=WWToLNuQQ_ext1

WZToLNu2Q_pow = Process()
WZToLNu2Q_pow.nickname = "WZToLNu2Q_pow"
WZToLNu2Q_pow.plotname = "W"
WZToLNu2Q_pow.weights = {"xsec":sampleDict[WZToLNu2Q_pow.nickname][1],"nevents":sampleDict[WZToLNu2Q_pow.nickname][3]}
WZToLNu2Q_pow.cuts = {"rare":""} 
WZToLNu2Q_pow.file = WZToLNu2Q_pow.nickname+"_2016.root"
WZToLNu2Q_pow.classification = sampleDict[WZToLNu2Q_pow.nickname][3]
WZToLNu2Q_pow.color = ["#CF8AC8"]    
HAA_processes[WZToLNu2Q_pow.nickname]=WZToLNu2Q_pow

WZTo3LNu = Process()
WZTo3LNu.nickname = "WZTo3LNu"
WZTo3LNu.plotname = "W"
WZTo3LNu.weights = {"xsec":sampleDict[WZTo3LNu.nickname][1],"nevents":sampleDict[WZTo3LNu.nickname][3]}
WZTo3LNu.cuts = {"rare":""} 
WZTo3LNu.file = WZTo3LNu.nickname+"_2016.root"
WZTo3LNu.classification = sampleDict[WZTo3LNu.nickname][3]
WZTo3LNu.color = ["#CF8AC8"]    
HAA_processes[WZTo3LNu.nickname]=WZTo3LNu

WZext = Process()
WZext.nickname = "WZ_ext1"
WZext.plotname = "TT"
WZext.weights = {"xsec":sampleDict[WZext.nickname][1],"nevents":sampleDict[WZext.nickname][3]}
WZext.cuts = {"rare":""} 
WZext.file = WZext.nickname+"_2016.root"
WZext.classification = sampleDict[WZext.nickname][3]
WZext.color = ["#CF8AC8"]    
HAA_processes[WZext.nickname]=WZext

WZ = Process()
WZ.nickname = "WZ"
WZ.plotname = "W"
WZ.weights = {"xsec":sampleDict[WZ.nickname][1],"nevents":sampleDict[WZ.nickname][3]}
WZ.cuts = {"rare":""} 
WZ.file = WZ.nickname+"_2016.root"
WZ.classification = sampleDict[WZ.nickname][3]
WZ.color = ["#CF8AC8"]    
HAA_processes[WZ.nickname]=WZ

WW = Process()
WW.nickname = "WW"
WW.plotname = "W"
WW.weights = {"xsec":sampleDict[WW.nickname][1],"nevents":sampleDict[WW.nickname][3]}
WW.cuts = {"rare":""} 
WW.file = WW.nickname+"_2016.root"
WW.classification = sampleDict[WW.nickname][3]
WW.color = ["#CF8AC8"]    
HAA_processes[WW.nickname]=WW

WW_ext1 = Process()
WW_ext1.nickname = "WW_ext1"
WW_ext1.plotname = "W"
WW_ext1.weights = {"xsec":sampleDict[WW_ext1.nickname][1],"nevents":sampleDict[WW_ext1.nickname][3]}
WW_ext1.cuts = {"rare":""} 
WW_ext1.file = WW_ext1.nickname+"_2016.root"
WW_ext1.classification = sampleDict[WW_ext1.nickname][3]
WW_ext1.color = ["#CF8AC8"]    
HAA_processes[WW_ext1.nickname]=WW_ext1

#WWext = Process()
#WWext.nickname = "WWext"
#WWext.plotname = "W"
#WWext.weights = {"xsec":sampleDict[WWext.nickname][1],"nevents":sampleDict[WWext.nickname][3]}
#WWext.cuts = {"rare":""} 
#WWext.file = WWext.nickname+"_2016.root"
#WWext.classification = sampleDict[WWext.nickname][3]
#WWext.color = ["#CF8AC8"]    
#HAA_processes[WWext.nickname]=WWext

#WWext1 = Process()
#WWext1.nickname = "WWext1"
#WWext1.plotname = "W"
#WWext1.weights = {"xsec":sampleDict[WWext1.nickname][1],"nevents":sampleDict[WWext1.nickname][3]}
#WWext1.cuts = {"rare":""} 
#WWext1.file = WWext1.nickname+"_2016.root"
#WWext1.classification = sampleDict[WWext1.nickname][3]
#WWext1.color = ["#CF8AC8"]    
#HAA_processes[WWext1.nickname]=WWext1

VBFHToTauTau = Process()
VBFHToTauTau.nickname = "VBFHToTauTau"
VBFHToTauTau.plotname = "rare"
VBFHToTauTau.weights = {"xsec":sampleDict[VBFHToTauTau.nickname][1],"nevents":sampleDict[VBFHToTauTau.nickname][3]}
VBFHToTauTau.cuts = {"rare":""} 
VBFHToTauTau.file = VBFHToTauTau.nickname+"_2016.root"
VBFHToTauTau.classification = sampleDict[VBFHToTauTau.nickname][3]
VBFHToTauTau.color = ["#CF8AC8"]    
HAA_processes[VBFHToTauTau.nickname]=VBFHToTauTau

WplusHToTauTau = Process()
WplusHToTauTau.nickname = "WplusHToTauTau"
WplusHToTauTau.plotname = "rare"
WplusHToTauTau.weights = {"xsec":sampleDict[WplusHToTauTau.nickname][1],"nevents":sampleDict[WplusHToTauTau.nickname][3]}
WplusHToTauTau.cuts = {"rare":""} 
WplusHToTauTau.file = WplusHToTauTau.nickname+"_2016.root"
WplusHToTauTau.classification = sampleDict[WplusHToTauTau.nickname][3]
WplusHToTauTau.color = ["#CF8AC8"]    
HAA_processes[WplusHToTauTau.nickname]=WplusHToTauTau

WminusHToTauTau = Process()
WminusHToTauTau.nickname = "WminusHToTauTau"
WminusHToTauTau.plotname = "rare"
WminusHToTauTau.weights = {"xsec":sampleDict[WminusHToTauTau.nickname][1],"nevents":sampleDict[WminusHToTauTau.nickname][3]}
WminusHToTauTau.cuts = {"rare":""} 
WminusHToTauTau.file = WminusHToTauTau.nickname+"_2016.root"
WminusHToTauTau.classification = sampleDict[WminusHToTauTau.nickname][3]
WminusHToTauTau.color = ["#CF8AC8"]    
HAA_processes[WminusHToTauTau.nickname]=WminusHToTauTau

#WWext = Process()
#WWext.nickname = "WWext1"
#WWext.plotname = "rare"
#WWext.weights = {"xsec":sampleDict[WWext.nickname][1],"nevents":sampleDict[WWext.nickname][3]}
#WWext.cuts = {"rare":""} 
#WWext.file = WWext.nickname+"_2016.root"
#WWext.classification = sampleDict[WWext.nickname][3]
#WWext.color = ["#CF8AC8"]    
#HAA_processes[WWext.nickname]=WWext

ggZH_HToNuNu_ZToLL = Process()
ggZH_HToNuNu_ZToLL.nickname = "ggZH_HToNuNu_ZToLL"
ggZH_HToNuNu_ZToLL.plotname = "rare"
ggZH_HToNuNu_ZToLL.weights = {"xsec":sampleDict[ggZH_HToNuNu_ZToLL.nickname][1],"nevents":sampleDict[ggZH_HToNuNu_ZToLL.nickname][3]}
ggZH_HToNuNu_ZToLL.cuts = {"rare":""} 
ggZH_HToNuNu_ZToLL.file = ggZH_HToNuNu_ZToLL.nickname+"_2016.root"
ggZH_HToNuNu_ZToLL.classification = sampleDict[ggZH_HToNuNu_ZToLL.nickname][3]
ggZH_HToNuNu_ZToLL.color = ["#CF8AC8"]    
HAA_processes[ggZH_HToNuNu_ZToLL.nickname]=ggZH_HToNuNu_ZToLL

ggZH_HToNuNu_ZToNuNu = Process()
ggZH_HToNuNu_ZToNuNu.nickname = "ggZH_HToNuNu_ZToNuNu"
ggZH_HToNuNu_ZToNuNu.plotname = "rare"
ggZH_HToNuNu_ZToNuNu.weights = {"xsec":sampleDict[ggZH_HToNuNu_ZToNuNu.nickname][1],"nevents":sampleDict[ggZH_HToNuNu_ZToNuNu.nickname][3]}
ggZH_HToNuNu_ZToNuNu.cuts = {"rare":""} 
ggZH_HToNuNu_ZToNuNu.file = ggZH_HToNuNu_ZToNuNu.nickname+"_2016.root"
ggZH_HToNuNu_ZToNuNu.classification = sampleDict[ggZH_HToNuNu_ZToNuNu.nickname][3]
ggZH_HToNuNu_ZToNuNu.color = ["#CF8AC8"]    
HAA_processes[ggZH_HToNuNu_ZToNuNu.nickname]=ggZH_HToNuNu_ZToNuNu

ggZH_HToTauTau_ZToQQ = Process()
ggZH_HToTauTau_ZToQQ.nickname = "ggZH_HToTauTau_ZToQQ"
ggZH_HToTauTau_ZToQQ.plotname = "rare"
ggZH_HToTauTau_ZToQQ.weights = {"xsec":sampleDict[ggZH_HToTauTau_ZToQQ.nickname][1],"nevents":sampleDict[ggZH_HToTauTau_ZToQQ.nickname][3]}
ggZH_HToTauTau_ZToQQ.cuts = {"rare":""} 
ggZH_HToTauTau_ZToQQ.file = ggZH_HToTauTau_ZToQQ.nickname+"_2016.root"
ggZH_HToTauTau_ZToQQ.classification = sampleDict[ggZH_HToTauTau_ZToQQ.nickname][3]
ggZH_HToTauTau_ZToQQ.color = ["#CF8AC8"]    
HAA_processes[ggZH_HToTauTau_ZToQQ.nickname]=ggZH_HToTauTau_ZToQQ


#adding extra weights to all files... 
for objkey in HAA_processes.keys():
    if not objkey=="data":
        HAA_processes[objkey].eventWeights = EventWeights
#        HAA_processes[objkey].cuts = {sampleDict[HAA_processes[objkey].nickname][0]:""}
        #HAA_processes[objkey].cuts["prompt"] = [["gen_match_3","!=",0],["gen_match_4","!=",0]] # doing this in other script 
        print HAA_processes[objkey].nickname+"_2016"
print "added SF weights"



del sf_MuonId
