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
import os
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

            survive = ops[cut[3]](tempvar,cut[4])
            #survive=True
            if not survive:
                #print cut[0],"     ",cut[1][0],"  ",getattr(evt,cut[1][0],None),"   ",cut[1][1],"   ",getattr(evt,cut[1][1],None),"  multicharge  ",tempvar,"   bool   ",ops[cut[3]](tempvar,cut[4])
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


def calculateHistos(functs,tree,HAA_Inc_mmmt,allcats,processObj,nickname,histodict,weightstring,commonweight,datadrivenPackage):
            

    newVarVals={}
    for var in HAA_Inc_mmmt.newvariables.keys():
        newVarVals[var]=0.0

    
    for ievt, evt in enumerate(tree):

        #do stuffs with  tree  
        if ievt % 1000 == 0:
            print "processing event ",ievt

        weightfinal = commonweight
        #find weight from string
        wtstring = 1.0
        for wts in weightstring:
            wtstring = wtstring * abs(getattr(evt,wts[0],None))

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
                        weightfinal =  weightfinal * (1 / float(weightDict[scalefactor]))
                    else:
                        weightfinal =  weightfinal * float(weightDict[scalefactor])
                    
                #combining all the cuts ... what about the new variables??
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


                if ((process=="FF_1" or process=="FF_2" or process=="FF_12") and (datadrivenPackage["bool"])):

                    ffLeg1Bool=cutStringBool(evt,processObj.cuts["FF_1"])
                    ffLeg2Bool=cutStringBool(evt,processObj.cuts["FF_2"])
                    ffLeg12Bool=cutStringBool(evt,processObj.cuts["FF_12"])

                    if ffLeg1Bool: 
                       #weightfinal = datadrivenPackage["fakerate1"].GetBinContent(datadrivenPackage["fakerate1"].FindBin(evt.pt_3))/(1.0000000001 - datadrivenPackage["fakerate1"].GetBinContent(datadrivenPackage["fakerate1"].FindBin(evt.pt_3)))
                        #weightfinal = datadrivenPackage["fitrate1"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate1"].GetParameter(0)) #const for now
                        #weightfinal = 0.025/(1-0.025) # from previous paper const
                        weightfinal = datadrivenPackage["pseudofit1"].Eval(evt.pt_3)/(1.0000001-datadrivenPackage["pseudofit1"].Eval(evt.pt_3)) #const for now
                    if ffLeg2Bool: 
                       #weightfinal = datadrivenPackage["fakerate2"].GetBinContent(datadrivenPackage["fakerate2"].FindBin(evt.pt_4))/(1.0000000001 - datadrivenPackage["fakerate2"].GetBinContent(datadrivenPackage["fakerate2"].FindBin(evt.pt_4)))
                        #weightfinal = datadrivenPackage["fitrate2"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate2"].GetParameter(0)) #const for now
                        #weightfinal = .17/(1-0.17) # from previous paper const
                        weightfinal = datadrivenPackage["pseudofit2"].Eval(evt.pt_4)/(1.0000001-datadrivenPackage["pseudofit2"].Eval(evt.pt_4)) #const for now
                    if ffLeg12Bool: 
                       #weightfinal = (datadrivenPackage["fakerate1"].GetBinContent(datadrivenPackage["fakerate1"].FindBin(evt.pt_3))/(1.0000000001 - datadrivenPackage["fakerate1"].GetBinContent(datadrivenPackage["fakerate1"].FindBin(evt.pt_3))))*(datadrivenPackage["fakerate2"].GetBinContent(datadrivenPackage["fakerate2"].FindBin(evt.pt_4))/(1.0000000001 - datadrivenPackage["fakerate2"].GetBinContent(datadrivenPackage["fakerate2"].FindBin(evt.pt_4))))
                        #weightfinal = (datadrivenPackage["fitrate1"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate1"].GetParameter(0)))*(datadrivenPackage["fitrate2"].GetParameter(0)/(1.0000001-datadrivenPackage["fitrate2"].GetParameter(0))) #const for now
                        #weightfinal = (0.025/(1-0.025))*(.17/(1-0.17))
                        weightfinal = (datadrivenPackage["pseudofit1"].Eval(evt.pt_3)/(1.0000001-datadrivenPackage["pseudofit1"].Eval(evt.pt_3)))*(datadrivenPackage["pseudofit2"].Eval(evt.pt_4)/(1.0000001-datadrivenPackage["pseudofit2"].Eval(evt.pt_4))) #const for now

                    if(ffLeg1Bool or ffLeg2Bool or ffLeg12Bool): 
                        for var in newVarVals.keys():
                            #print cat.newvariables[var][1]
                            #print var+":"+cat.name+":"+process
                            temp = functs[cat.newvariables[var][0]](evt,cat.newvariables[var][2])

                            newhistodict[var+":"+cat.name+":"+process].Fill(functs[cat.newvariables[var][0]](evt,cat.newvariables[var][2]),float(weightfinal))
                            if type(temp)==float:
                                newhistodict[var+":"+cat.name+":"+process].Fill(functs[cat.newvariables[var][0]](evt,cat.newvariables[var][2]),float(weightfinal))
                            if type(temp)==list:
                                print "for tlorentz objects"
                                for var,ivar in enumerate(temp):
                                    newhistodict[var+":"+cat.name+":"+process].Fill(functs[cat.newvariables[var][0]](evt,cat.newvariables[var][2]),float(weightfinal))

                        for variableHandle in cat.vars.keys():
                            val = getattr(evt,cat.vars[variableHandle][0],None)
                            filedict[variableHandle].cd()
                            filedict[variableHandle].cd(cat.name)
                            if not val==None:
                                histodict[variableHandle+":"+cat.name+":"+process].Fill(float(val),float(weightfinal))

                    continue  # for FF we don't want to compare to regular cuts


                if ((process!="FF" or process!="FF_1" or process!="FF_2" or process!="FF_12")and (survive==True)):

                    #fill the new variables 
                    for var in newVarVals.keys():
                        #print cat.newvariables[var][1]
                        #print var+":"+cat.name+":"+process
                        temp = functs[cat.newvariables[var][0]](evt,cat.newvariables[var][2])

                        newhistodict[var+":"+cat.name+":"+process].Fill(functs[cat.newvariables[var][0]](evt,cat.newvariables[var][2]),float(weightfinal))
                        if type(temp)==float:
                            newhistodict[var+":"+cat.name+":"+process].Fill(functs[cat.newvariables[var][0]](evt,cat.newvariables[var][2]),float(weightfinal))
                        if type(temp)==list:
                            print "for tlorentz objects"
                            for var,ivar in enumerate(temp):
                                newhistodict[var+":"+cat.name+":"+process].Fill(functs[cat.newvariables[var][0]](evt,cat.newvariables[var][2]),float(weightfinal))

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
    import argparse

    parser = argparse.ArgumentParser(description="This file generates root files containing Histograms ... files in utils contain selections and settings")
    #parser.add_arguement('--CategoryFiles',nargs="+",help="Select the files containing the categories for the datacards")
    parser.add_argument("-o",  "--outname", default="",  help="postfix string")
    parser.add_argument("-dd",  "--datadriven", default=False,action='store_true',  help="Use DataDriven Method")
    args = parser.parse_args()


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
    
    os.mkdir("out"+str(args.outname))
    filedict = {}
    histodict = {}
    newhistodict = {}
    #cat=HAA_Inc_mmmt

    numvar=0
    for variableHandle in HAA_Inc_mmmt.vars.keys():
    
        #variable=fullvariable.split(":") # for the unrolled cases

        #filedict[variable[0]]=ROOT.TFile.Open(cat.name+"_"+str(variable[0])+".root","RECREATE")
        filedict[variableHandle]=ROOT.TFile.Open("out"+str(args.outname)+"/"+str(variableHandle)+".root","RECREATE")
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
        filedict[variable]=ROOT.TFile.Open("out"+str(args.outname)+"/"+str(variable)+".root","RECREATE")
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
                        #tmpbin = np.asarray(cat.newvariablesbins[numvar])
                        bins = HAA_Inc_mmmt.newvariables[variable][1]
                        tmpbin = np.asarray(bins)
                        #print variable+":"+cat.name+":"+process
                        if  HAA_Inc_mmmt.newvariables[variable][-1]=="multiob":
                            for var,ivar in enumerate(temp):
                                newhistodict[variable+"_"+str(ivar)+":"+cat.name+":"+process] = ROOT.TH1D(str(process),str(process),len(tmpbin[ivar])-1,tmpbin[ivar])
                                newhistodict[variable+"_"+str(ivar)+":"+cat.name+":"+process].Write(str(process),ROOT.TObject.kOverwrite)
                        else:
                            newhistodict[variable+":"+cat.name+":"+process] = ROOT.TH1D(str(process),str(process),len(tmpbin)-1,tmpbin)
                            newhistodict[variable+":"+cat.name+":"+process].Write(str(process),ROOT.TObject.kOverwrite)
                        #newhistodict[variable+":"+cat.name+":"+process] = ROOT.TH1D(str(process),str(process),len(tmpbin)-1,tmpbin)
                        #newhistodict[variable+":"+cat.name+":"+process].Write(str(process),ROOT.TObject.kOverwrite)
        numvar=numvar+1

    #gather extra global variables or weights
    
    weight = CommonWeights["lumi"][0]
    weightstring = CommonWeights["string"]
   
    #weight = weight+"*"+CommonWeights["mcweight"][0]

    #Calculate the scale factors and fill the histograms 
    #print sampleDict.keys()

    datadriven = args.datadriven
    datadrivenPackage={}
    datadrivenPackage["bool"]=datadriven
    if datadriven:
        #load necessary histograms
            ff_file_3 = ROOT.TFile.Open("FFhistos/pt_3.root")
            ff_file_4 = ROOT.TFile.Open("FFhistos/pt_4.root")
            #sstight_3 = ff_file_3.Get("mmmt_inclusive_FF_tight_SS/data_obs")
            #sstight_4 = ff_file_4.Get("mmmt_inclusive_FF_tight_SS/data_obs")
            #ssloose_3 = ff_file_3.Get("mmmt_inclusive_FF_loose_SS/data_obs")
            #ssloose_4 = ff_file_4.Get("mmmt_inclusive_FF_loose_SS/data_obs")
            #osloose_3 = ff_file_3.Get("mmmt_inclusive_FF_loose_OS/data_obs")
            #osloose_4 = ff_file_4.Get("mmmt_inclusive_FF_loose_OS/data_obs")
            ss_1_tight = ff_file_3.Get("mmmt_inclusive_FF_SS_1_tight/data_obs")
            ss_1_loose = ff_file_3.Get("mmmt_inclusive_FF_SS_1_loose/data_obs")
            ss_2_tight = ff_file_4.Get("mmmt_inclusive_FF_SS_2_tight/data_obs")
            ss_2_loose = ff_file_4.Get("mmmt_inclusive_FF_SS_2_loose/data_obs")
            
            #f_1 = sstight_3.Clone()
            f_1= ss_1_tight.Clone()
            #f_1.Divide(ssloose_3)
            f_1.Divide(ss_1_loose)
            f_1.SetName("FakeRateMuLeg")
            f_2 = ss_2_tight.Clone()
            f_2.Divide(ss_2_loose)
            f_2.SetName("FakeRateTauLeg")

            #tf_1 = ROOT.TF1("tf_1","[0]*expo[1]",f_1.GetXaxis().GetXmin(),f_1.GetXaxis().GetXmax()) 
            #tf_2 = ROOT.TF1("tf_2","[0]*expo[1]",f_2.GetXaxis().GetXmin(),f_2.GetXaxis().GetXmax()) 
            tf_1 = ROOT.TF1("tf_1","[0]",f_1.GetXaxis().GetXmin(),f_1.GetXaxis().GetXmax()) 
            tf_2 = ROOT.TF1("tf_2","[0]",f_2.GetXaxis().GetXmin(),f_2.GetXaxis().GetXmax()) 

            fakemeasurefile = ROOT.TFile.Open("FFhistos/fakemeasure.root","RECREATE")
            fakemeasurefile.cd()
            f_1.Fit("tf_1")
            f_2.Fit("tf_2")
            f_1.Write(f_1.GetName(),ROOT.TObject.kOverwrite)
            f_2.Write(f_2.GetName(),ROOT.TObject.kOverwrite)
            tf_1.Write(tf_1.GetName(),ROOT.TObject.kOverwrite)
            tf_2.Write(tf_2.GetName(),ROOT.TObject.kOverwrite)
            datadrivenPackage["fakerate1"]=f_1
            datadrivenPackage["fakerate2"]=f_2
            datadrivenPackage["fitrate1"]=tf_1
            datadrivenPackage["fitrate2"]=tf_2
            applyf = ROOT.TFile.Open("FFhistos/ffmeas.root","READ")
            applyf.cd()
            pg1 = applyf.Get("mupoints")
            pg2 = applyf.Get("taupoints")
            pf1 = ROOT.TF1("pf1","[0]*exp([1]*x)+[2]",0,100) 
            pf2 = ROOT.TF1("pf2","[0]*exp([1]*x)+[2]",0,100) 
            pg1.Fit("pf1")
            pg2.Fit("pf2")
            datadrivenPackage["pseudofit1"]=pf1
            datadrivenPackage["pseudofit2"]=pf2

    for nickname in filelist.keys():

        fin = ROOT.TFile.Open(dir+filelist[nickname],"read")
        print "working on file ",fin.GetName()
        tree = fin.Get(treename)
        #newtree=tree.CloneTree()
        print "tree entries ",tree.GetEntries()

        #processes will be a dictionary of key: process name value: cuts
        #processes = HAA_processes[nickname].cuts
        processObj = HAA_processes[nickname]

        calculateHistos(functs,tree,HAA_Inc_mmmt,allcats,processObj,nickname,histodict,weightstring,weight,datadrivenPackage)
        
     
                
    print "writing the histograms to the file"
    for varcatPro in histodict.keys():
        #print "file ",filedict[varcatPro.split(":")[0]].GetName()," contents ",filedict[varcatPro.split(":")[0]].ls()
        #print "file ",varcatPro.split(":")[0]," writing ",varcatPro,"  final entries ",histodict[varcatPro].GetEntries()
        filedict[varcatPro.split(":")[0]].cd()
        filedict[varcatPro.split(":")[0]].cd(varcatPro.split(":")[1])
        histodict[varcatPro].Write(histodict[varcatPro].GetName(),ROOT.TObject.kOverwrite)
    print "writing the histograms to the file"
    for varcatPro in newhistodict.keys():
        #print "file ",filedict[varcatPro.split(":")[0]].GetName()," contents ",filedict[varcatPro.split(":")[0]].ls()
        #print "file ",varcatPro.split(":")[0]," writing ",varcatPro,"  final entries ",newhistodict[varcatPro].GetEntries()
        filedict[varcatPro.split(":")[0]].cd()
        filedict[varcatPro.split(":")[0]].cd(varcatPro.split(":")[1])
        newhistodict[varcatPro].Write(newhistodict[varcatPro].GetName(),ROOT.TObject.kOverwrite)
            

