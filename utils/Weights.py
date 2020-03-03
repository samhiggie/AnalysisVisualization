#########################
#Author: Sam Higginbotham
'''

* File Name :Weights.py

* Purpose :Contains the mapping between weights and processes

* Creation Date : 05-02-2020

* Last Modified :

'''
#########################
CommonWeights={}
CommonWeights["lumi"]=["59740"]
CommonWeights["mcweight"]=["__WEIGHT__*GENWEIGHT*puweight"]


SMHTTWeights={}
SMHTTWeights["data_obs"]=[""]
SMHTTWeights["TTT"]=CommonWeights["mcweight"]
SMHTTWeights["TTL"]=CommonWeights["mcweight"]
SMHTTWeights["ZTT"]=CommonWeights["mcweight"]
SMHTTWeights["ZTL"]=CommonWeights["mcweight"]
SMHTTWeights["ZL"]=CommonWeights["mcweight"]

