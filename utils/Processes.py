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
from Parametrization import Process

#This is a list of Process objects
HAA_processes={}
HAA_processes_test={}

#Gathering extra information

#Gather the analysis datasets and info 
sampleDict = {}

with open("/afs/cern.ch/work/s/shigginb/cmssw/HAA/ZH_Run2_10_2_9/src/AnalysisVisualization/MCsamples_2016.csv")  as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        #[nickname]        = [category,xsec,numberOfEvents,finishedEvents,idk?,DASDataset]
        sampleDict[row[0]] = [row[1],row[2],row[3],row[4],row[5],row[6]]

#HAA_processes_auto = {}
#
#for nickname in sampleDict.keys():
#    HAA_processes_auto[nickname] = Process()
#    
#    HAA_processes_auto[nickname].nickname=nickname
#    HAA_processes_auto[nickname].plotname=nickname

#AToZhToLLTauTau_M-220_2016.root  HToAAToMuMuTauTau_M50_2016.root  TTW_2016.root
#AToZhToLLTauTau_M-240_2016.root  HToAAToMuMuTauTau_M55_2016.root  TTZH_2016.root
#AToZhToLLTauTau_M-280_2016.root  HToAAToMuMuTauTau_M60_2016.root  TTZZ_2016.root
#AToZhToLLTauTau_M-340_2016.root  ST_s_2016.root                   W1JetsToLNu_2016.root
#AToZhToLLTauTau_M-350_2016.root  ST_tW_antitop_2016.root          W2JetsToLNu_2016.root
#DY1JetsToLL_2016.root            ST_tW_top_2016.root              W2JetsToLNuext_2016.root
#DY2JetsToLL_2016.root            ST_t_antitop_2016.root           W3JetsToLNu_2016.root
#DY3JetsToLL_2016.root            ST_t_top_2016.root               W3JetsToLNuext_2016.root
#DY4JetsToLL_2016.root            TTGamma_1LT_2016.root            W4JetsToLNu_2016.root
#DYJetsToLLM10to50_2016.root      TTGamma_1LTbar_2016.root         W4JetsToLNuext1_2016.root
#DYJetsToLL_2016.root             TTGamma_2L_2016.root             W4JetsToLNuext2_2016.root
#DYJetsToLLext_2016.root          TTGamma_Hadr_2016.root           WJetsToLNu_2016.root
#EWKWMinus2Jets_2016.root         TTHH_2016.root                   WJetsToLNuext_2016.root
#EWKWPlus2Jets_2016.root          TTJets_1LT_2016.root             WW_2016.root
#EWKZ2JetsLL_2016.root            TTJets_1LTbar_2016.root          WWext_2016.root
#EWKZ2JetsNuNu_2016.root          TTJets_DiLept_2016.root          WZext_2016.root
#HToAAToMuMuTauTau_M15_2016.root  TTTJ_2016.root                   WminusHToTauTau_2016.root
#HToAAToMuMuTauTau_M20_2016.root  TTTT_2016.root                   WplusHToTauTau_2016.root
#HToAAToMuMuTauTau_M25_2016.root  TTTW_2016.root                   ZHToTauTau_2016.root
#HToAAToMuMuTauTau_M30_2016.root  TTWH_2016.root                   ZZ_2016.root
#HToAAToMuMuTauTau_M35_2016.root  TTWJetsToQQ_2016.root            ZZext_2016.root
#HToAAToMuMuTauTau_M40_2016.root  TTWW_2016.root
#HToAAToMuMuTauTau_M45_2016.root  TTWZ_2016.root


a40test = Process()
a40test.nickname = "HToAAToMuMuTauTau_M40"
a40test.plotname = "a40"
a40test.weights = {"xsec":1,"nevents":250000,"theoryXsec":(31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a40test.cuts = {"a40":""}
a40test.file = "HToAAToMuMuTauTau_M40_2016.root"
a40test.classification = sampleDict[a40test.nickname][3]
a40test.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes_test[a40test.nickname]=a40test
    
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
from Categories import ffbasecuts
FF = Process()
FF.nickname = "FF"
FF.plotname = "FF"
FF.weights = {"xsec":1.0} 
FF.cuts = {"FF":"","FF_1":[ffbasecuts[0],["mediumId_3","<",1],["isGlobal_3","<",1],["isTracker_3","<",1],["idDeepTau2017v2p1VSjet_4",">=",2.]],"FF_2":[ffbasecuts[0],["mediumId_3",">=",1],[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]],["idDeepTau2017v2p1VSjet_4","<",2.]],"FF_12":[ffbasecuts[0]]}
FF.file = "data.root"
FF.classification = "FF"
#FF.color = [ROOT.kBlack]    #do i need a list here??? to embed it
HAA_processes[FF.nickname]=FF

a15 = Process()
a15.nickname = "HToAAToMuMuTauTau_M15"
a15.plotname = "a15"
a15.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a15.cuts = {"a15":""}
a15.file = "HToAAToMuMuTauTau_M15_2016.root"
a15.classification = sampleDict[a15.nickname][3]
a15.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a15.nickname]=a15

a20 = Process()
a20.nickname = "HToAAToMuMuTauTau_M20"
a20.plotname = "a20"
a20.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a20.cuts = {"a20":""}
a20.file = "HToAAToMuMuTauTau_M20_2016.root"
a20.classification = sampleDict[a20.nickname][3]
a20.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a20.nickname]=a20

a25 = Process()
a25.nickname = "HToAAToMuMuTauTau_M25"
a25.plotname = "a25"
a25.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a25.cuts = {"a25":""}
a25.file = "HToAAToMuMuTauTau_M25_2016.root"
a25.classification = sampleDict[a25.nickname][3]
a25.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a25.nickname]=a25

a30 = Process()
a30.nickname = "HToAAToMuMuTauTau_M30"
a30.plotname = "a30"
a30.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a30.cuts = {"a30":""}
a30.file = "HToAAToMuMuTauTau_M30_2016.root"
a30.classification = sampleDict[a30.nickname][3]
a30.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a30.nickname]=a30

a35 = Process()
a35.nickname = "HToAAToMuMuTauTau_M35"
a35.plotname = "a35"
a35.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a35.cuts = {"a35":""}
a35.file = "HToAAToMuMuTauTau_M35_2016.root"
a35.classification = sampleDict[a35.nickname][3]
a35.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a35.nickname]=a35

a40 = Process()
a40.nickname = "HToAAToMuMuTauTau_M40"
a40.plotname = "a40"
a40.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a40.cuts = {"a40":""}
a40.file = "HToAAToMuMuTauTau_M40_2016.root"
a40.classification = sampleDict[a40.nickname][3]
a40.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a40.nickname]=a40

a45 = Process()
a45.nickname = "HToAAToMuMuTauTau_M45"
a45.plotname = "a45"
a45.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a45.cuts = {"a45":""}
a45.file = "HToAAToMuMuTauTau_M45_2016.root"
a45.classification = sampleDict[a45.nickname][3]
a45.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a45.nickname]=a45

a50 = Process()
a50.nickname = "HToAAToMuMuTauTau_M50"
a50.plotname = "a50"
a50.cuts = {"a50":""}
a50.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a50.file = "HToAAToMuMuTauTau_M50_2016.root"
a50.classification = sampleDict[a50.nickname][3]
a50.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a50.nickname]=a50

a55 = Process()
a55.nickname = "HToAAToMuMuTauTau_M55"
a55.plotname = "a55"
a55.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a55.cuts = {"a55":""}
a55.file = "HToAAToMuMuTauTau_M55_2016.root"
a55.classification = sampleDict[a55.nickname][3]
a55.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a55.nickname]=a55

a60 = Process()
a60.nickname = "HToAAToMuMuTauTau_M60"
a60.plotname = "a60"
a60.weights = {"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a60.cuts = {"a60":""}
a60.file = "HToAAToMuMuTauTau_M60_2016.root"
a60.classification = sampleDict[a60.nickname][0]
a60.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a60.nickname]=a60

GluGluToContinToZZTo2mu2tau = Process()
GluGluToContinToZZTo2mu2tau.nickname = "GluGluToContinToZZTo2mu2tau"
GluGluToContinToZZTo2mu2tau.plotname = "GluGluToContinToZZTo2mu2tau"
GluGluToContinToZZTo2mu2tau.weights = {"xsec":sampleDict[GluGluToContinToZZTo2mu2tau.nickname][1],"nevents":sampleDict[GluGluToContinToZZTo2mu2tau.nickname][3]} 
GluGluToContinToZZTo2mu2tau.cuts = {"GluGluToContinToZZTo2mu2tau":""}
GluGluToContinToZZTo2mu2tau.file = GluGluToContinToZZTo2mu2tau.nickname+"_2016.root"
GluGluToContinToZZTo2mu2tau.classification = sampleDict[GluGluToContinToZZTo2mu2tau.nickname][3]
GluGluToContinToZZTo2mu2tau.color = ["#65E114"]    #do i need a list here??? to embed it
HAA_processes[GluGluToContinToZZTo2mu2tau.nickname]=GluGluToContinToZZTo2mu2tau

ttZ = Process()
ttZ.nickname = "ttZ"
ttZ.plotname = "ttZ"
ttZ.weights = {"xsec":sampleDict[ttZ.nickname][1],"nevents":sampleDict[ttZ.nickname][3]} 
ttZ.cuts = {"ttZ":""}
ttZ.file = ttZ.nickname+"_2016.root"
ttZ.classification = sampleDict[ttZ.nickname][3]
ttZ.color = ["#65E114"]    #do i need a list here??? to embed it
HAA_processes[ttZ.nickname]=ttZ

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

ZZext = Process()
ZZext.nickname = "ZZext"
ZZext.plotname = "ZZ"
ZZext.weights = {"xsec":sampleDict[ZZext.nickname][1],"nevents":sampleDict[ZZext.nickname][3]} 
ZZext.cuts = {"ZZ":""}
ZZext.file = ZZext.nickname+"_2016.root"
ZZext.classification = sampleDict[ZZext.nickname][3]
ZZext.color = ["#13E2FE"]    #do i need a list here??? to embed it
HAA_processes[ZZext.nickname]=ZZext

DY = Process()
DY.nickname = "DYJetsToLL"
DY.plotname = "ZTT"
DY.weights = {"xsec":sampleDict[DY.nickname][1],"nevents":sampleDict[DY.nickname][3]} #higgs xsec [pb] * 2hdm type Branching ratio
DY.cuts = {"Z":"","ZTT":[["gen_match_4","==",5]],"ZL":[["gen_match_4",">=",5]],"ZJ":[["gen_match_4",">",5]]} #higgs xsec [pb] * 2hdm type Branching ratio
DY.file = DY.nickname+"_2016.root"
DY.classification = sampleDict[DY.nickname][3]
DY.color = ["#CF8AC8"]    #do i need a list here??? to embed it
HAA_processes[DY.nickname]=DY

DYext = Process()
DYext.nickname = "DYJetsToLLext"
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

WJets = Process()
WJets.nickname = "WJetsToLNu"
WJets.plotname = "W"
WJets.weights = {"xsec":sampleDict[WJets.nickname][1],"nevents":sampleDict[WJets.nickname][3]}
WJets.cuts = {"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]} 
WJets.file = WJets.nickname+"_2016.root"
WJets.classification = sampleDict[WJets.nickname][3]
WJets.color = ["#CF8AC8"]    
HAA_processes[WJets.nickname]=WJets

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
W2Jetsext.nickname = "W2JetsToLNuext"
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
W3Jetsext.nickname = "W3JetsToLNuext"
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
W4Jetsext1.nickname = "W4JetsToLNuext1"
W4Jetsext1.plotname = "W"
W4Jetsext1.weights = {"xsec":sampleDict[W4Jetsext1.nickname][1],"nevents":sampleDict[W4Jetsext1.nickname][3]}
W4Jetsext1.cuts = {"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]} 
W4Jetsext1.file = W4Jetsext1.nickname+"_2016.root"
W4Jetsext1.classification = sampleDict[W4Jetsext1.nickname][3]
W4Jetsext1.color = ["#CF8AC8"]    
HAA_processes[W4Jetsext1.nickname]=W4Jetsext1

W4Jetsext2 = Process()
W4Jetsext2.nickname = "W4JetsToLNuext2"
W4Jetsext2.plotname = "W"
W4Jetsext2.weights = {"xsec":sampleDict[W4Jetsext2.nickname][1],"nevents":sampleDict[W4Jetsext2.nickname][3]}
W4Jetsext2.cuts = {"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]} 
W4Jetsext2.file = W4Jetsext2.nickname+"_2016.root"
W4Jetsext2.classification = sampleDict[W4Jetsext2.nickname][3]
W4Jetsext2.color = ["#CF8AC8"]    
HAA_processes[W4Jetsext2.nickname]=W4Jetsext2

TTJets_1LT = Process()
TTJets_1LT.nickname = "TTJets_1LT"
TTJets_1LT.plotname = "TT"
TTJets_1LT.weights = {"xsec":sampleDict[TTJets_1LT.nickname][1],"nevents":sampleDict[TTJets_1LT.nickname][3]}
TTJets_1LT.cuts = {"TT":"","TTT":[["gen_match_4","==",5]],"TTL":[["gen_match_4",">=",5]],"TTJ":[["gen_match_4",">",5]]} 
TTJets_1LT.file = TTJets_1LT.nickname+"_2016.root"
TTJets_1LT.classification = sampleDict[TTJets_1LT.nickname][3]
TTJets_1LT.color = ["#CF8AC8"]    
HAA_processes[TTJets_1LT.nickname]=TTJets_1LT

TTJets_1LTbar = Process()
TTJets_1LTbar.nickname = "TTJets_1LTbar"
TTJets_1LTbar.plotname = "T"
TTJets_1LTbar.weights = {"xsec":sampleDict[TTJets_1LTbar.nickname][1],"nevents":sampleDict[TTJets_1LTbar.nickname][3]}
TTJets_1LTbar.cuts = {"TT":"","TTT":[["gen_match_4","==",5]],"TTL":[["gen_match_4",">=",5]],"TTJ":[["gen_match_4",">",5]]} 
TTJets_1LTbar.file = TTJets_1LTbar.nickname+"_2016.root"
TTJets_1LTbar.classification = sampleDict[TTJets_1LTbar.nickname][3]
TTJets_1LTbar.color = ["#CF8AC8"]    
HAA_processes[TTJets_1LTbar.nickname]=TTJets_1LTbar

TTJets_DiLept = Process()
TTJets_DiLept.nickname = "TTJets_DiLept"
TTJets_DiLept.plotname = "T"
TTJets_DiLept.weights = {"xsec":sampleDict[TTJets_DiLept.nickname][1],"nevents":sampleDict[TTJets_DiLept.nickname][3]}
TTJets_DiLept.cuts = {"TT":"","TTT":[["gen_match_4","==",5]],"TTL":[["gen_match_4",">=",5]],"TTJ":[["gen_match_4",">",5]]} 
TTJets_DiLept.file = TTJets_DiLept.nickname+"_2016.root"
TTJets_DiLept.classification = sampleDict[TTJets_DiLept.nickname][3]
TTJets_DiLept.color = ["#CF8AC8"]    
HAA_processes[TTJets_DiLept.nickname]=TTJets_DiLept

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
ST_tW_antitop.nickname = "ST_tW_antitop"
ST_tW_antitop.plotname = "ST"
ST_tW_antitop.weights = {"xsec":sampleDict[ST_tW_antitop.nickname][1],"nevents":sampleDict[ST_tW_antitop.nickname][3]}
ST_tW_antitop.cuts = {"ST":"","STT":[["gen_match_4","==",5]],"STL":[["gen_match_4",">=",5]],"STJ":[["gen_match_4",">",5]]} 
ST_tW_antitop.file = ST_tW_antitop.nickname+"_2016.root"
ST_tW_antitop.classification = sampleDict[ST_tW_antitop.nickname][3]
ST_tW_antitop.color = ["#CF8AC8"]    
HAA_processes[ST_tW_antitop.nickname]=ST_tW_antitop

ST_tW_top = Process()
ST_tW_top.nickname = "ST_tW_top"
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

EWKWPlus2Jets = Process()
EWKWPlus2Jets.nickname = "EWKWPlus2Jets"
EWKWPlus2Jets.plotname = "EWK"
EWKWPlus2Jets.weights = {"xsec":sampleDict[EWKWPlus2Jets.nickname][1],"nevents":sampleDict[EWKWPlus2Jets.nickname][3]}
EWKWPlus2Jets.cuts = {"EWK":""} 
EWKWPlus2Jets.file = EWKWPlus2Jets.nickname+"_2016.root"
EWKWPlus2Jets.classification = sampleDict[EWKWPlus2Jets.nickname][3]
EWKWPlus2Jets.color = ["#CF8AC8"]    
HAA_processes[EWKWPlus2Jets.nickname]=EWKWPlus2Jets


EWKZ2JetsLL = Process()
EWKZ2JetsLL.nickname = "EWKZ2JetsLL"
EWKZ2JetsLL.plotname = "EWK"
EWKZ2JetsLL.weights = {"xsec":sampleDict[EWKZ2JetsLL.nickname][1],"nevents":sampleDict[EWKZ2JetsLL.nickname][3]}
EWKZ2JetsLL.cuts = {"EWK":""} 
EWKZ2JetsLL.file = EWKZ2JetsLL.nickname+"_2016.root"
EWKZ2JetsLL.classification = sampleDict[EWKZ2JetsLL.nickname][3]
EWKZ2JetsLL.color = ["#CF8AC8"]    
HAA_processes[EWKZ2JetsLL.nickname]=EWKZ2JetsLL

EWKZ2JetsNuNu = Process()
EWKZ2JetsNuNu.nickname = "EWKZ2JetsNuNu"
EWKZ2JetsNuNu.plotname = "EWK"
EWKZ2JetsNuNu.weights = {"xsec":sampleDict[EWKZ2JetsNuNu.nickname][1],"nevents":sampleDict[EWKZ2JetsNuNu.nickname][3]}
EWKZ2JetsNuNu.cuts = {"EWK":""} 
EWKZ2JetsNuNu.file = EWKZ2JetsNuNu.nickname+"_2016.root"
EWKZ2JetsNuNu.classification = sampleDict[EWKZ2JetsNuNu.nickname][3]
EWKZ2JetsNuNu.color = ["#CF8AC8"]    
HAA_processes[EWKZ2JetsNuNu.nickname]=EWKZ2JetsNuNu

#Rare ??? man there are many... 
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
TTHH.nickname = "TTHH"
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

TTTJ = Process()
TTTJ.nickname = "TTTJ"
TTTJ.plotname = "TT"
TTTJ.weights = {"xsec":sampleDict[TTTJ.nickname][1],"nevents":sampleDict[TTTJ.nickname][3]}
TTTJ.cuts = {"rare":""} 
TTTJ.file = TTTJ.nickname+"_2016.root"
TTTJ.classification = sampleDict[TTTJ.nickname][3]
TTTJ.color = ["#CF8AC8"]    
HAA_processes[TTTJ.nickname]=TTTJ

TTTW = Process()
TTTW.nickname = "TTTW"
TTTW.plotname = "TT"
TTTW.weights = {"xsec":sampleDict[TTTW.nickname][1],"nevents":sampleDict[TTTW.nickname][3]}
TTTW.cuts = {"rare":""} 
TTTW.file = TTTW.nickname+"_2016.root"
TTTW.classification = sampleDict[TTTW.nickname][3]
TTTW.color = ["#CF8AC8"]    
HAA_processes[TTTW.nickname]=TTTW

TTW = Process()
TTW.nickname = "TTW"
TTW.plotname = "TT"
TTW.weights = {"xsec":sampleDict[TTW.nickname][1],"nevents":sampleDict[TTW.nickname][3]}
TTW.cuts = {"rare":""} 
TTW.file = TTW.nickname+"_2016.root"
TTW.classification = sampleDict[TTW.nickname][3]
TTW.color = ["#CF8AC8"]    
HAA_processes[TTW.nickname]=TTW

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
TTWW.nickname = "TTWW"
TTWW.plotname = "TT"
TTWW.weights = {"xsec":sampleDict[TTWW.nickname][1],"nevents":sampleDict[TTWW.nickname][3]}
TTWW.cuts = {"rare":""} 
TTWW.file = TTWW.nickname+"_2016.root"
TTWW.classification = sampleDict[TTWW.nickname][3]
TTWW.color = ["#CF8AC8"]    
HAA_processes[TTWW.nickname]=TTWW

TTWZ = Process()
TTWZ.nickname = "TTWZ"
TTWZ.plotname = "TT"
TTWZ.weights = {"xsec":sampleDict[TTWZ.nickname][1],"nevents":sampleDict[TTWZ.nickname][3]}
TTWZ.cuts = {"rare":""} 
TTWZ.file = TTWZ.nickname+"_2016.root"
TTWZ.classification = sampleDict[TTWZ.nickname][3]
TTWZ.color = ["#CF8AC8"]    
HAA_processes[TTWZ.nickname]=TTWZ

TTZH = Process()
TTZH.nickname = "TTZH"
TTZH.plotname = "TT"
TTZH.weights = {"xsec":sampleDict[TTZH.nickname][1],"nevents":sampleDict[TTZH.nickname][3]}
TTZH.cuts = {"rare":""} 
TTZH.file = TTZH.nickname+"_2016.root"
TTZH.classification = sampleDict[TTZH.nickname][3]
TTZH.color = ["#CF8AC8"]    
HAA_processes[TTZH.nickname]=TTZH

TTZZ = Process()
TTZZ.nickname = "TTZZ"
TTZZ.plotname = "TT"
TTZZ.weights = {"xsec":sampleDict[TTZZ.nickname][1],"nevents":sampleDict[TTZZ.nickname][3]}
TTZZ.cuts = {"rare":""} 
TTZZ.file = TTZZ.nickname+"_2016.root"
TTZZ.classification = sampleDict[TTZZ.nickname][3]
TTZZ.color = ["#CF8AC8"]    
HAA_processes[TTZZ.nickname]=TTZZ

WW = Process()
WW.nickname = "WW"
WW.plotname = "TT"
WW.weights = {"xsec":sampleDict[WW.nickname][1],"nevents":sampleDict[WW.nickname][3]}
WW.cuts = {"rare":""} 
WW.file = WW.nickname+"_2016.root"
WW.classification = sampleDict[WW.nickname][3]
WW.color = ["#CF8AC8"]    
HAA_processes[WW.nickname]=WW

WWext = Process()
WWext.nickname = "WWext"
WWext.plotname = "TT"
WWext.weights = {"xsec":sampleDict[WWext.nickname][1],"nevents":sampleDict[WWext.nickname][3]}
WWext.cuts = {"rare":""} 
WWext.file = WWext.nickname+"_2016.root"
WWext.classification = sampleDict[WWext.nickname][3]
WWext.color = ["#CF8AC8"]    
HAA_processes[WWext.nickname]=WWext

WZext = Process()
WZext.nickname = "WZext"
WZext.plotname = "TT"
WZext.weights = {"xsec":sampleDict[WZext.nickname][1],"nevents":sampleDict[WZext.nickname][3]}
WZext.cuts = {"rare":""} 
WZext.file = WZext.nickname+"_2016.root"
WZext.classification = sampleDict[WZext.nickname][3]
WZext.color = ["#CF8AC8"]    
HAA_processes[WZext.nickname]=WZext
#runII legacy SMHTT
SMHTT = {}
#data 
SMHTT["muDATA.root"] = ["data_obs"]


#signals
SMHTT["ggH125.root"]=["ggH"]
SMHTT["vbfH125.root"]=["qqH"]
#SMHTT["ZH"]="ggH125.root"
#SMHTT["WH"]="WH125.root"

#backgrounds 
SMHTT["ZJETS.root"]=["ZTT","ZL"]
SMHTT["FF.root"]=["jetFakes"]
SMHTT["DiBoson.root"]=["VVL","VVT"]
SMHTT["TT.root"]=["TTL","TTT"]
#SMHTT["W"]="WJETS.root"
