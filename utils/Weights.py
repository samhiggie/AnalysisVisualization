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
#CommonWeights["mcweight"]=["__WEIGHT__*GENWEIGHT*puweight"]
CommonWeights["mcweight"]=["LHEweight"]


SMHTTWeights={}
SMHTTWeights["data_obs"]=[""]
SMHTTWeights["TTT"]=CommonWeights["mcweight"]
SMHTTWeights["TTL"]=CommonWeights["mcweight"]
SMHTTWeights["ZTT"]=CommonWeights["mcweight"]
SMHTTWeights["ZTL"]=CommonWeights["mcweight"]
SMHTTWeights["ZL"]=CommonWeights["mcweight"]

#Where do I place energy scale shifts? 
HAAWeights={}
HAAWeights["data_obs"]=[""]
HAAWeights["TTT"]=CommonWeights["mcweight"]
HAAWeights["TTL"]=CommonWeights["mcweight"]
HAAWeights["ZTT"]=CommonWeights["mcweight"]
HAAWeights["ZTL"]=CommonWeights["mcweight"]
HAAWeights["ZL"]=CommonWeights["mcweight"]
HAAWeights["TES"]=["(((pt_4*1.011)>40&&decayMode_4==0)||((pt_4*0.995)>40&&decayMode_4==1)||((pt_4*0.988)>40&&decayMode_4==10))&&(((pt_3*1.011)>50&&decayMode_3==0)||((pt_3*0.995)>50&&decayMode_3==1)||((pt_3*0.988)>50&&decayMode_3==10))"]
