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


    

a15 = Process()
a15.nickname = "HToAAToMuMuTauTau_M15"
a15.plotname = "a15"
a15.weights = {"xsec":1,"nevents":250000,"theoryXsec":(16*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a15.cuts = {"a15":""}
a15.file = "HToAAToMuMuTauTau_M15_2016.root"
a15.classification = sampleDict[a15.nickname][3]
a15.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a15.nickname]=a15

a20 = Process()
a20.nickname = "HToAAToMuMuTauTau_M20"
a20.plotname = "a20"
a20.weights = {"xsec":1,"nevents":250000,"theoryXsec":(16*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a20.cuts = {"a20":""}
a20.file = "HToAAToMuMuTauTau_M20_2016.root"
a20.classification = sampleDict[a20.nickname][3]
a20.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a20.nickname]=a20

a25 = Process()
a25.nickname = "HToAAToMuMuTauTau_M25"
a25.plotname = "a25"
a25.weights = {"xsec":1,"nevents":250000,"theoryXsec":(16*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a25.cuts = {"a25":""}
a25.file = "HToAAToMuMuTauTau_M25_2016.root"
a25.classification = sampleDict[a25.nickname][3]
a25.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a25.nickname]=a25

a30 = Process()
a30.nickname = "HToAAToMuMuTauTau_M30"
a30.plotname = "a30"
a30.weights = {"xsec":1,"nevents":250000,"theoryXsec":(16*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a30.cuts = {"a30":""}
a30.file = "HToAAToMuMuTauTau_M30_2016.root"
a30.classification = sampleDict[a30.nickname][3]
a30.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a30.nickname]=a30

a35 = Process()
a35.nickname = "HToAAToMuMuTauTau_M35"
a35.plotname = "a35"
a35.weights = {"xsec":1,"nevents":250000,"theoryXsec":(16*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a35.cuts = {"a35":""}
a35.file = "HToAAToMuMuTauTau_M35_2016.root"
a35.classification = sampleDict[a35.nickname][3]
a35.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a35.nickname]=a35

a40 = Process()
a40.nickname = "HToAAToMuMuTauTau_M40"
a40.plotname = "a40"
a40.weights = {"xsec":1,"nevents":250000,"theoryXsec":(16*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a40.cuts = {"a40":""}
a40.file = "HToAAToMuMuTauTau_M40_2016.root"
a40.classification = sampleDict[a40.nickname][3]
a40.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a40.nickname]=a40

a45 = Process()
a45.nickname = "HToAAToMuMuTauTau_M45"
a45.plotname = "a45"
a45.weights = {"xsec":1,"nevents":250000,"theoryXsec":(16*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a45.cuts = {"a45":""}
a45.file = "HToAAToMuMuTauTau_M45_2016.root"
a45.classification = sampleDict[a45.nickname][3]
a45.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a45.nickname]=a45

a50 = Process()
a50.nickname = "HToAAToMuMuTauTau_M50"
a50.plotname = "a50"
a50.cuts = {"a50":""}
a50.weights = {"xsec":1,"nevents":250000,"theoryXsec":(16*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a50.file = "HToAAToMuMuTauTau_M50_2016.root"
a50.classification = sampleDict[a50.nickname][3]
a50.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a50.nickname]=a50

a55 = Process()
a55.nickname = "HToAAToMuMuTauTau_M55"
a55.plotname = "a55"
a55.weights = {"xsec":1,"nevents":250000,"theoryXsec":(16*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a55.cuts = {"a55":""}
a55.file = "HToAAToMuMuTauTau_M55_2016.root"
a55.classification = sampleDict[a55.nickname][3]
a55.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a55.nickname]=a55

a60 = Process()
a60.nickname = "HToAAToMuMuTauTau_M60"
a60.plotname = "a60"
a60.weights = {"xsec":1,"nevents":250000,"theoryXsec":(16*0.00005)} #higgs xsec [pb] * 2hdm type Branching ratio
a60.cuts = {"a60":""}
a60.file = "HToAAToMuMuTauTau_M60_2016.root"
a60.classification = sampleDict[a60.nickname][0]
a60.color = [ROOT.kRed]    #do i need a list here??? to embed it
HAA_processes[a60.nickname]=a60

ZZ = Process()
ZZ.nickname = "ZZ"
ZZ.plotname = "ZZ"
ZZ.weights = {"xsec":sampleDict[ZZ.nickname][1],"nevents":sampleDict[ZZ.nickname][3]} 
ZZ.cuts = {"ZZ":""}
ZZ.file = ZZ.nickname+"_2016.root"
ZZ.classification = sampleDict[ZZ.nickname][0]
ZZ.color = ["#13E2FE"]    #do i need a list here??? to embed it
HAA_processes[ZZ.nickname]=ZZ

ZZext = Process()
ZZext.nickname = "ZZext"
ZZext.plotname = "ZZ"
ZZext.weights = {"xsec":sampleDict[ZZext.nickname][1],"nevents":sampleDict[ZZext.nickname][3]} 
ZZext.cuts = {"ZZ":""}
ZZext.file = ZZext.nickname+"_2016.root"
ZZext.classification = sampleDict[ZZext.nickname][0]
ZZext.color = ["#13E2FE"]    #do i need a list here??? to embed it
HAA_processes[ZZext.nickname]=ZZext

DY = Process()
DY.nickname = "DYJetsToLL"
DY.plotname = "ZTT"
DY.weights = {"xsec":sampleDict[DY.nickname][1],"nevents":sampleDict[DY.nickname][3]} #higgs xsec [pb] * 2hdm type Branching ratio
DY.cuts = {"Z":"","ZTT":"gen_match_4==5","ZL":"gen_match_4>=5","ZJ":"gen_match_4>5"} #higgs xsec [pb] * 2hdm type Branching ratio
DY.file = DY.nickname+"_2016.root"
DY.classification = sampleDict[DY.nickname][3]
DY.color = ["#CF8AC8"]    #do i need a list here??? to embed it
HAA_processes[DY.nickname]=DY

DYext = Process()
DYext.nickname = "DYJetsToLLext"
DYext.plotname = "ZTT"
DYext.weights = {"xsec":sampleDict[DYext.nickname][1],"nevents":sampleDict[DYext.nickname][3]} #higgs xsec [pb] * 2hdm type Branching ratio
DYext.cuts = {"Z":"","ZTT":"gen_match_4==5","ZL":"gen_match_4>=5","ZJ":"gen_match_4>5"} #higgs xsec [pb] * 2hdm type Branching ratio
DYext.file = DYext.nickname+"_2016.root"
DYext.classification = sampleDict[DYext.nickname][3]
DYext.color = ["#CF8AC8"]    #do i need a list here??? to embed it
HAA_processes[DYext.nickname]=DYext

DY1 = Process()
DY1.nickname = "DY1JetsToLL"
DY1.plotname = "ZTT"
DY1.weights = {"xsec":sampleDict[DY1.nickname][1],"nevents":sampleDict[DY1.nickname][3]} #higgs xsec [pb] * 2hdm type Branching ratio
DY1.cuts = {"Z":"","ZTT":"gen_match_4==5","ZL":"gen_match_4>=5","ZJ":"gen_match_4>5"} #higgs xsec [pb] * 2hdm type Branching ratio
DY1.file = DY1.nickname+"_2016.root"
DY1.classification = sampleDict[DY1.nickname][3]
DY1.color = ["#CF8AC8"]    #do i need a list here??? to embed it
HAA_processes[DY1.nickname]=DY1

DY2 = Process()
DY2.nickname = "DY2JetsToLL"
DY2.plotname = "ZTT"
DY2.weights = {"xsec":sampleDict[DY2.nickname][1],"nevents":sampleDict[DY2.nickname][3]} #higgs xsec [pb] * 2hdm type Branching ratio
DY2.cuts = {"Z":"","ZTT":"gen_match_4==5","ZL":"gen_match_4>=5","ZJ":"gen_match_4>5"} #higgs xsec [pb] * 2hdm type Branching ratio
DY2.file = DY2.nickname+"_2016.root"
DY2.classification = sampleDict[DY2.nickname][3]
DY2.color = ["#CF8AC8"]    #do i need a list here??? to embed it
HAA_processes[DY2.nickname]=DY2

DY3 = Process()
DY3.nickname = "DY3JetsToLL"
DY3.plotname = "ZTT"
DY3.weights = {"xsec":sampleDict[DY3.nickname][1],"nevents":sampleDict[DY3.nickname][3]} #higgs xsec [pb] * 2hdm type Branching ratio
DY3.cuts = {"Z":"","ZTT":"gen_match_4==5","ZL":"gen_match_4>=5","ZJ":"gen_match_4>5"} #higgs xsec [pb] * 2hdm type Branching ratio
DY3.file = DY3.nickname+"_2016.root"
DY3.classification = sampleDict[DY3.nickname][3]
DY3.color = ["#CF8AC8"]    #do i need a list here??? to embed it
HAA_processes[DY3.nickname]=DY3

DY4 = Process()
DY4.nickname = "DY4JetsToLL"
DY4.plotname = "ZTT"
DY4.weights = {"xsec":sampleDict[DY4.nickname][1],"nevents":sampleDict[DY4.nickname][3]} #higgs xsec [pb] * 2hdm type Branching ratio
DY4.cuts = {"Z":"","ZTT":"gen_match_4==5","ZL":"gen_match_4>=5","ZJ":"gen_match_4>5"} #higgs xsec [pb] * 2hdm type Branching ratio
DY4.file = DY4.nickname+"_2016.root"
DY4.classification = sampleDict[DY4.nickname][3]
DY4.color = ["#CF8AC8"]    #do i need a list here??? to embed it
HAA_processes[DY4.nickname]=DY4

DYM = Process()
DYM.nickname = "DYJetsToLLM10to50"
DYM.plotname = "ZTT"
DYM.weights = {"xsec":sampleDict[DYM.nickname][1],"nevents":sampleDict[DYM.nickname][3]} #higgs xsec [pb] * 2hdm type Branching ratio
DYM.cuts = {"Z":"","ZTT":"gen_match_4==5","ZL":"gen_match_4>=5","ZJ":"gen_match_4>5"} #higgs xsec [pb] * 2hdm type Branching ratio
DYM.file = DYM.nickname+"_2016.root"
DYM.classification = sampleDict[DYM.nickname][3]
DYM.color = ["#CF8AC8"]    #do i need a list here??? to embed it
HAA_processes[DYM.nickname]=DYM




#old method
HAA = {}
#signals
HAA["HToAAToMuMuTauTau_M15_2016.root"]=["a15"]
HAA["HToAAToMuMuTauTau_M20_2016.root"]=["a20"]
HAA["HToAAToMuMuTauTau_M25_2016.root"]=["a25"]
HAA["HToAAToMuMuTauTau_M30_2016.root"]=["a30"]
HAA["HToAAToMuMuTauTau_M35_2016.root"]=["a35"]
HAA["HToAAToMuMuTauTau_M40_2016.root"]=["a40"]
HAA["HToAAToMuMuTauTau_M45_2016.root"]=["a45"]
HAA["HToAAToMuMuTauTau_M50_2016.root"]=["a50"]
HAA["HToAAToMuMuTauTau_M55_2016.root"]=["a55"]
HAA["HToAAToMuMuTauTau_M60_2016.root"]=["a60"]

#backgrounds 
HAA["ZJETS.root"]=["Z","ZJ","ZL","ZTT"]
HAA["WJETS.root"]=["W","WJ","WL"]
HAA["WZ.root"]=["WZ","WZJ","WZL"]
HAA["ZZ.root"]=["ZZ"]
#HAA["EWK"]=["EWK"]
HAA["TT.root"]=["TTL","TTT","TTJ"]







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
