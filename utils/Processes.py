#########################
#Author: Sam Higginbotham
'''

* File Name : Processes.py

* Purpose : Stores dictionaries for analysis signal and background files 

* Creation Date : 04-02-2020

* Last Modified :

'''
#########################

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
