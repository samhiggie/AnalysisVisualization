# !/usr/bin/env python
#########################
#Author: Sam Higginbotham
'''

* File Name : MakeDataCards_HAA.py

* Purpose : For Datacard creation. Root file containing histograms that can be used with Combine. This will skim the events 

* Creation Date : 04-02-2020

* Last Modified :

'''
#########################
import sys
import ROOT
import numpy as np
#import argsparse
from array import array
import itertools
import operator
import csv


#GLOBAL VARIABLES AND FUNCTIONS!


#Setting up operators for cut string iterator ... very useful
ops = { "==": operator.eq, "!=": operator.eq, ">": operator.gt, "<": operator.lt, ">=": operator.ge, "<=": operator.le, "band": operator.and_,"bor":operator.or_}

def subArrayLogic(evt,array):
    boo=ops[array[1]](getattr(evt,array[0],None),array[2])
    return boo

def cutStringBool(evt,cuts):
    survive=False
    #print "cuts",cuts
    for cut in cuts:

        #recursive case ... keep er going...  
        if type(cut[0][0])==type([0]):
            sruvive=cutStringBool(evt,cuts)

        #print "single cut",cut,"  len ",len(cut)
        #Here there are cut types ... everything will be product ands... but these are the subtypes 
        # "OR" so we can AND(OR)AND(OR) 
        if cut[0][0]=="OR":
            tempbool=False
            for oe in range(1,len(cut)):
                #print cut[oe]
                if ops[cut[oe][1]](getattr(evt,cut[oe][0],None),cut[oe][2]):
                    tempbool=True
                    break
                    #print "passed trigger "
                
            survive=tempbool    
            if not survive:
                break
            else: continue
        # "EQT" does a particular set of variables when applied in a function pass a cut?
        if cut[0][0]=="EQT":
            #print cut[0]
            for i,var in enumerate(cut[1]):
                if cut[2]=="mult":
                    if i==0:
                        tempvar = 1.0
                    tempvar = tempvar * getattr(evt,var,None)
                if cut[2]=="div":
                    if(i==0 and getattr(evt,var,None)!=0):
                        tempvar = 1.0
                    tempvar = tempvar / getattr(evt,var,None)
                if cut[2]=="add":
                    if(i==0):
                        tempvar = 0.0
                    tempvar = tempvar+getattr(evt,var,None)
                if cut[2]=="sub":
                    if(i==0):
                        tempvar = 0.0
                    tempvar = tempvar - getattr(evt,var,None)

            if ops[cut[3]](tempvar,cut[4]):
                survive=True
            if not survive:
                break
            else: continue 
        # "IF" statments in order to make cuts after variables in an event pass a selection criteria ... "THEN" included only for asethetics 
        if cut[0][0]=="IF":
            tempbool=True
            for i,state in enumerate(cut[1]):
                #if not ops[state[1]](getattr(evt,state[0],None),state[2]):
                #print state,ops[state[1]](getattr(evt,state[0],None),state[2])
                if ops[state[1]](getattr(evt,state[0],None),state[2]):
                    tempbool=True
                else:
                    tempbool=False
                    break
            if tempbool==True:
                for i,state in enumerate(cut[3]):
                    #print state,ops[state[1]](getattr(evt,state[0],None),state[2]),getattr(evt,state[0],None)
                    if ops[state[1]](getattr(evt,state[0],None),state[2]):
                        survive=True 
                    else:
                        survive=False 
                        break
            else:
                #print "no if statement soo no cuts"
                continue 

            if not survive:
                break
            else: continue 

        # Default CASE! does it pass the cut?
        # should  provide this as separate function??? because the absg does exist above... 
        else:
            statement = cut[1]
            #print statement
            if statement=="absg":
                if ops["<"](getattr(evt,cut[0],None),(-1 * cut[2])): 
                    survive=True
                if ops[">"](getattr(evt,cut[0],None),cut[2]): 
                    survive=True 
                else: 
                    survive=False
                    break
                    #print "failed ",cut
            if statement=="absl":
                if (ops[">"](getattr(evt,cut[0],None),(-1 * cut[2])) and ops["<"](getattr(evt,cut[0],None),cut[2])): 
                    survive=True
                else: 
                    survive=False
                    break
                    
                    #print "failed ",cut
            else:  
                if ops[statement](getattr(evt,cut[0],None),cut[2]):
                    survive=True
                    #print "passed cut",cut
                else:
                    survive=False
                    #return survive
                    break
                    #print "failed ",cut
        
    return survive 

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

def runSVFit(entry, channel) :
		  
    measuredMETx = entry.met*cos(entry.metphi)
    measuredMETy = entry.met*sin(entry.metphi)

    #define MET covariance
    covMET = TMatrixD(2,2)
    covMET[0][0] = entry.metcov00
    covMET[1][0] = entry.metcov10
    covMET[0][1] = entry.metcov01
    covMET[1][1] = entry.metcov11


    #self.kUndefinedDecayType, self.kTauToHadDecay,  self.kTauToElecDecay, self.kTauToMuDecay = 0, 1, 2, 3
    if channel == 'et' :
	measTau1 = ROOT.MeasuredTauLepton(kTauToElecDecay, entry.pt_3, entry.eta_3, entry.phi_3, 0.000511) 
    elif channel == 'mt' :
	measTau1 = ROOT.MeasuredTauLepton(kTauToMuDecay, entry.pt_3, entry.eta_3, entry.phi_3, 0.106) 
    elif channel == 'tt' :
	measTau1 = ROOT.MeasuredTauLepton(kTauToHadDecay, entry.pt_3, entry.eta_3, entry.phi_3, entry.m_3)
		    
    if channel != 'em' :
	measTau2 = ROOT.MeasuredTauLepton(kTauToHadDecay, entry.pt_4, entry.eta_4, entry.phi_4, entry.m_4)

    if channel == 'em' :
	measTau1 = ROOT.MeasuredTauLepton(kTauToElecDecay,  entry.pt_3, entry.eta_3, entry.phi_3, 0.000511)
	measTau2 = ROOT.MeasuredTauLepton(kTauToMuDecay, entry.pt_4, entry.eta_4, entry.phi_4, 0.106)

    VectorOfTaus = ROOT.std.vector('MeasuredTauLepton')
    instance = VectorOfTaus()
    instance.push_back(measTau1)
    instance.push_back(measTau2)

    FMTT = ROOT.FastMTT()
    FMTT.run(instance, measuredMETx, measuredMETy, covMET)
    ttP4 = FMTT.getBestP4()
    return ttP4.M(), ttP4.Mt() 

def createUnrolledHistogram(cat,numvar,filelist,process,variable,weight,filedict):

    filedict[variable[0]].cd()
    filedict[variable[0]].mkdir(cat.name)
    filedict[variable[0]].cd(cat.name)

    #creating primary bins
    primaryBins = np.asarray(cat.binning[numvar][0]) #always the first binning set 
    #creating secondary bin list
    secondaryBins = []

    for bins in range(1,len(cat.binning[numvar])):
        secondaryBins.append(cat.binning[numvar][bins])

    binlist = returnBinList(secondaryBins)
    combos = returnPerms(binlist)
    #print "Primary Bins  ",primaryBins
    #print "Combinations for secondary bins  ",combos
    procut = HAA_ProCuts.get(process)
    histList = []
    urnum = 0 
    for bincuts in combos:
        htmp=ROOT.TH1D("htemp"+str(process)+str(urnum),"htemp"+str(process)+str(urnum),len(primaryBins)-1,primaryBins)
        print "working on variable   ",str(variable[0]),"for htemp"+str(process)+str(urnum)
        unrollCut=""

        for extraVarnum in range(1,len(variable)):
            unrollCut += "&&"+str(variable[extraVarnum])+">"+str(bincuts[extraVarnum-1][0])+"&&"+str(variable[extraVarnum])+"<"+str(bincuts[extraVarnum-1][1])
        if procut:
            allcuts = "("+str(cat.cuts)+"&&"+str(cat.cuts["preselection"])+"&&"+str(cat.cuts["trigger"])+"&&"+str(HAA_ProCuts[process])+unrollCut+")"
            #print "entries after cut  ",tree.GetEntries(allcuts)
            tree.Draw(str(variable[0])+">>"+"htemp"+str(process)+str(urnum),allcuts+weight)
        else:
            allcuts = "("+str(cat.cuts)+"&&"+str(cat.cuts["preselection"])+"&&"+str(cat.cuts["trigger"])+unrollCut+")"
            #print "entries after cut  ",tree.GetEntries(allcuts)
            tree.Draw(str(variable[0])+">>"+"htemp"+str(process)+str(urnum),allcuts+weight)
         
        histList.append(htmp)
        htmp.Write(htmp.GetName(),ROOT.TObject.kOverwrite)
        urnum+=1
        #print "htmp entries  ",htmp.GetEntries()
        #print allcuts
        htmp=0
    #next create a list of histograms and try to chain them together bin by bin. 
    masterbins = cat.binning[numvar][0]*len(combos)
    print masterbins
    masterUnroll = ROOT.TH1D(str(process)+"_unrolled",str(process)+"_unrolled",len(masterbins),0,len(masterbins))
    masterBinNum=0 
    for hist in histList:
        for bin in range(0,hist.GetNbinsX()+1):
            masterUnroll.SetBinContent(masterBinNum,hist.GetBinContent(bin))
            #print "bin    ",masterBinNum,"    bin content ",hist.GetBinContent(bin)
            masterUnroll.SetBinError(masterBinNum,hist.GetBinError(bin))
            masterBinNum+=1
           
    
    
    masterUnroll.Write(masterUnroll.GetName(),ROOT.TObject.kOverwrite)
    return


def calculateHistos(functs,tree,HAA_Inc_mmmt,allcats,processObj,nickname,histodict,weightstring,commonweight):
            

    #do I need to use branches?
    #need a temp variable per new variable to fill each event
    newVarVals={}
    for var in HAA_Inc_mmmt.newvariables.keys():
        newVarVals[var]=0.0

    
    for ievt, evt in enumerate(tree):
        #copy the standard weights... will add scale factors to this later 

        #do stuffs with  tree  
        if ievt % 1000 == 0:
            print "processing event ",ievt
        
        #this won't work... how do I update variables per event?? just make a new tree and shit?? but that is basically nanoAOD!
        #evt.pt_1 = evt.pt_1 +30.0
        #new variable mll - mtt
        #mllmtt = evt.mll - evt.mtt 
        #mllmtt = getattr(evt,"mll",None) - getattr(evt,"m_vis",None)
        #newVars["mll-mtt"]= getattr(evt,"mll",None) - getattr(evt,"m_vis",None)
        #newVars["mll-mtt"]= fun["-"][getattr(evt,"mll",None),getattr(evt,"m_vis",None)]

        weightfinal = commonweight
        #find weight from string
        wtstring = 1.0
        for wts in weightstring:
            wtstring = wtstring * abs(getattr(evt,wts[0],None))
        #weightfinal = weightfinal * wtstring
        #fill histograms 

        for cat in allcats:

            for process in processObj.cuts.keys():
                #Cuts
                procut = processObj.cuts[process]

                #Weights addition to common
                weightDict = processObj.weights
                weightfinal = commonweight
                for scalefactor in weightDict.keys():
                    if scalefactor == "nevents":
                        #weightfinal =  weightfinal+"*"+"("+"1/"+str(weightDict[scalefactor])+")"
                        weightfinal =  weightfinal * (1 / float(weightDict[scalefactor]))
                    else:
                        #weightfinal =  weightfinal+"*"+str(weightDict[scalefactor])
                        weightfinal =  weightfinal * float(weightDict[scalefactor])
                    
                #combining all the cuts 
                cuts = []
                for cuttype in cat.cuts.keys():
                    for cut in cat.cuts[cuttype]:
                        cuts.append(cut)
                #cuts.append(procut)
                #print "allcuts   ",cuts

                #Operator expansion cut string
                survive=cutStringBool(evt,cuts)

                if process=="data_obs":
                    weightfinal = 1.0   #don't weight the data!!
                
                if survive==True:

                    #fill the new variables 
                    for var in newVarVals.keys():
                        #print cat.newvariables[var][1]
                        #print var+":"+cat.name+":"+process
                        newhistodict[var+":"+cat.name+":"+process].Fill(functs[cat.newvariables[var][0]](evt,cat.newvariables[var][1]),float(weightfinal))

                    #fill the current variables
                    #for variable in cat.variables:
                    for variableHandle in cat.vars.keys():
                        #move this outside variables
      
                        #obtaining the right variable... WHAT about changed variables in the event!!?
                        val = getattr(evt,cat.vars[variableHandle][0],None)
                        filedict[variableHandle].cd()
                        filedict[variableHandle].cd(cat.name)
                        #Fix the weight problem!!!
                        if not val==None:
                            histodict[variableHandle+":"+cat.name+":"+process].Fill(float(val),float(weightfinal))
                        #if process=="data_obs":
                        #    histodict[variable[0]+":"+cat.name+":"+process].Fill(float(val))
                        #else:
                        #    histodict[variable[0]+":"+cat.name+":"+process].Fill(float(val),float(weightfinal))
                        #histodict[variable[0]+":"+cat.name+":"+process].Fill(float(val))
                


    return

if __name__ == "__main__":
    


    #Structure for plotting variables 
    from utils.Parametrization import Category
    from utils.Parametrization import Process


    #Structure for mapping between root files and processes  
    from utils.Processes import HAA_processes
    #for testing only... one cat a15
    #from utils.Processes import HAA_processes_test
    from utils.ProcessCuts import HAA_ProCuts

    #Structure for mapping between processes and weights
    from utils.Weights import CommonWeights
    from utils.Weights import HAAWeights

    #Structure for mapping processes and extra distributions 
    #from utils.Systematics import systematics

    #parser = argsparse.ArgumentParser(description="make datacards")
    #parser.add_arguement('--CategoryFiles',nargs="+",help="Select the files containing the categories for the datacards")
    #args = parser.parse_args()
    #loader = 
    categories = []

    #gather all the analysis categories
    from utils.Categories import allcats
    from utils.Categories import HAA_Inc_mmmt
    
    #gather functions for computing variables in the event loop
    from utils.functions import functs


    #This is where the plotting takes place!
    #categories.append(HAA_Inc_mmmt) 
    #categories = allcats

    #Gather the Analysis Files
    #dir = "/eos/user/s/shigginb/HAA_ntuples/March_2020/"
    dir = "/eos/home-s/shigginb/HAA_ntuples/March_2020/"
    #filelist = HAA
    filelist = {}

    #for testing only... one cat a15
    #HAA_processes=HAA_processes_test

    for proObj in HAA_processes.keys():
       filelist[proObj]=HAA_processes[proObj].file
    print "The files ",filelist
    treename = "Events"

    #Gather the analysis datasets and info 
    sampleDict = {}
 
    with open("MCsamples_2016.csv")  as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            sampleDict[row[0]] = [row[1],row[2],row[3],row[4],row[5],row[6]]

    print "initializing histograms"
    

    filedict = {}
    histodict = {}
    newhistodict = {}
    #cat=HAA_Inc_mmmt

    numvar=0
    for variableHandle in HAA_Inc_mmmt.vars.keys():
    
        #variable=fullvariable.split(":") # for the unrolled cases

        #filedict[variable[0]]=ROOT.TFile.Open(cat.name+"_"+str(variable[0])+".root","RECREATE")
        filedict[variableHandle]=ROOT.TFile.Open("out/"+str(variableHandle)+".root","RECREATE")
        #print "on variableHandle ",variableHandle
        filedict[variableHandle].cd()
        for cat in allcats:
            #print "Working on cat ",cat.name
            filedict[variableHandle].mkdir(cat.name)
            filedict[variableHandle].cd(cat.name)

            for nickname in filelist.keys():

                processObj = HAA_processes[nickname]

                for process in processObj.cuts.keys():


                    if filedict[variableHandle].Get(cat.name).GetListOfKeys().Contains(str(process)):
                        continue
                    else:
                        #bins = cat.binning[numvar]
                        bins = cat.vars[variableHandle][1]
                        if type(bins[0])==list:
                            histodict[variableHandle+":"+cat.name+":"+process] = ROOT.TH1D(str(process),str(process),bins[0][0],bins[0][1],bins[0][2])
                            histodict[variableHandle+":"+cat.name+":"+process].Write(str(process),ROOT.TObject.kOverwrite)
                        else:
                            tmpbin = np.asarray(bins)
                            histodict[variableHandle+":"+cat.name+":"+process] = ROOT.TH1D(str(process),str(process),len(tmpbin)-1,tmpbin)
                            histodict[variableHandle+":"+cat.name+":"+process].Write(str(process),ROOT.TObject.kOverwrite)
        numvar=numvar+1

    #gathering new varibles
    numvar=0
    for variable in HAA_Inc_mmmt.newvariables.keys():

        #filedict[variable]=ROOT.TFile.Open(cat.name+"_"+str(variable)+".root","RECREATE")
        filedict[variable]=ROOT.TFile.Open("out/"+str(variable)+".root","RECREATE")
        #print "on variable ",variable
        filedict[variable].cd()
        for cat in allcats:
            filedict[variable].mkdir(cat.name)
            filedict[variable].cd(cat.name)

            for nickname in filelist.keys():

                processObj = HAA_processes[nickname]

                for process in processObj.cuts.keys():

                    if filedict[variable].Get(cat.name).GetListOfKeys().Contains(str(process)):
                        continue
                    else:
                        tmpbin = np.asarray(cat.newvariablesbins[numvar])
                        #print variable+":"+cat.name+":"+process
                        newhistodict[variable+":"+cat.name+":"+process] = ROOT.TH1D(str(process),str(process),len(tmpbin)-1,tmpbin)
                        newhistodict[variable+":"+cat.name+":"+process].Write(str(process),ROOT.TObject.kOverwrite)
        numvar=numvar+1

    #gather extra global variables or weights
    
    weight = CommonWeights["lumi"][0]
    weightstring = CommonWeights["string"]
   
    #weight = weight+"*"+CommonWeights["mcweight"][0]

    #Calculate the scale factors and fill the histograms 
    #print sampleDict.keys()
    for nickname in filelist.keys():

        fin = ROOT.TFile.Open(dir+filelist[nickname],"read")
        print "working on file ",fin.GetName()
        tree = fin.Get(treename)
        #newtree=tree.CloneTree()
        print "tree entries ",tree.GetEntries()

        #processes will be a dictionary of key: process name value: cuts
        #processes = HAA_processes[nickname].cuts
        processObj = HAA_processes[nickname]

        calculateHistos(functs,tree,HAA_Inc_mmmt,allcats,processObj,nickname,histodict,weightstring,weight)
        
     
                
    print "writing the histograms to the file"
    for varcatPro in histodict.keys():
        #print "file ",filedict[varcatPro.split(":")[0]].GetName()," contents ",filedict[varcatPro.split(":")[0]].ls()
        print "file ",varcatPro.split(":")[0]," writing ",varcatPro,"  final entries ",histodict[varcatPro].GetEntries()
        filedict[varcatPro.split(":")[0]].cd()
        filedict[varcatPro.split(":")[0]].cd(varcatPro.split(":")[1])
        histodict[varcatPro].Write(histodict[varcatPro].GetName(),ROOT.TObject.kOverwrite)
    print "writing the histograms to the file"
    for varcatPro in newhistodict.keys():
        #print "file ",filedict[varcatPro.split(":")[0]].GetName()," contents ",filedict[varcatPro.split(":")[0]].ls()
        print "file ",varcatPro.split(":")[0]," writing ",varcatPro,"  final entries ",newhistodict[varcatPro].GetEntries()
        filedict[varcatPro.split(":")[0]].cd()
        filedict[varcatPro.split(":")[0]].cd(varcatPro.split(":")[1])
        newhistodict[varcatPro].Write(newhistodict[varcatPro].GetName(),ROOT.TObject.kOverwrite)
            

