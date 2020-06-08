#########################
#Author: Sam Higginbotham
'''

* File Name :functions.py 

* Purpose : The Categories in the analysis... typically a set of extra cuts that define a phase space region of interest and stores the variables of interest in the fit and visualization

* Creation Date : 05-02-2020

* Last Modified :

'''
#########################
import ROOT
#load these functions into a dictionary so we can call on them in the event loop

functs={}
minus = lambda evt,x: getattr(evt,x[0],None) - getattr(evt,x[1],None)
functs["minus"]=minus 
plus = lambda evt,x: getattr(evt,x[0],None) + getattr(evt,x[1],None)
functs["plus"]=plus
def multi(evt,x): 
    product = 1.0
    for i,j in enumerate(x):
        product = product * getattr(evt,x[i],None)
    #if product>0:
    #    print evt.evt,x[0],x[1],"  multi   ",product
    return product 
functs["multi"] = multi

def TLorentz(evt,x):
    vector = ROOT.TLorentzVector()
    if type(x[3])==float:
	    vector.SetPtEtaPhiM(getattr(evt,x[0],None),getattr(evt,x[1],None),getattr(evt,x[2],None),x[3])  #set the mass of the TLorentz Vector directly
    else: 
	    vector.SetPtEtaPhiM(getattr(evt,x[0],None),getattr(evt,x[1],None),getattr(evt,x[2],None),getattr(evt,x[3],None))
    vectorList.append(vector[0])
    vectorList.append(vector[1])
    vectorList.append(vector[2])
    vectorList.append(vector[3])
    return vectorList
functs["TLorentz"] = TLorentz

def TLorentzMT(vector):
    return vector.Mt()
functs["TLorentzMT"] = TLorentzMT

def addTLorentz(allvectors):
    sumvector = ROOT.TLorentzVector()
    for vector in allvectors:
        sumvector = sumvector + vector
    return sumvector
functs["addTLorentz"] = addTLorentz

