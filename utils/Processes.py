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
from Parametrization import Process

HAA = {}
#data 
#HAA["muDATA.root"] = ["data_obs"]
#This is a list of Process objects
#HAA={}
#a20 = Process()
#a20.name = "a20"
#a20.file = "all_HToAAToMuMuTauTau_M20_001.root"
#a20.classification = "signal"
#HAA[a20.file] = a20

#signals
HAA["all_HToAAToMuMuTauTau_M20_001.root"]=["a20"]
#HAA["ggH125.root"]=["ggH"]
#HAA["vbfH125.root"]=["qqH"]
#HAA["ZH"]="ggH125.root"
#HAA["WH"]="WH125.root"

#backgrounds 
#HAA["ZJETS.root"]=["ZTT","ZL"]
#HAA["FF.root"]=["jetFakes"]
#HAA["DiBoson.root"]=["VVL","VVT"]
#HAA["TT.root"]=["TTL","TTT"]
#HAA["W"]="WJETS.root"







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
