# !/usr/bin/env python
#########################
#Author: Sam Higginbotham
'''

* File Name : MakeDataCards.py

* Purpose : For Datacard creation. Root file containing histograms that can be used with Combine.

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
ops = { "==": operator.eq, ">": operator.gt, "<": operator.lt, ">=": operator.ge, "<=": operator.le, "band": operator.and_,"bor":operator.or_}

def cutStringBool(evt,cuts):
    survive=False
    #print "cuts",cuts
    for cut in cuts:
        #print "single cut",cut,"  len ",len(cut)
        if cut[0][0]=="OR":
            tempbool=False
            for oe in range(1,len(cut)):
                #print cut[oe]
                if ops[cut[oe][1]](getattr(evt,cut[oe][0],None),cut[oe][2]):
                    tempbool=True
                    #print "passed trigger "
                
            survive=tempbool    
        else:
            statement = cut[1]
            if statement=="absg":
                if ops["<"](getattr(evt,cut[0],None),(-1 * cut[2])): 
                    survive=True
                if ops[">"](getattr(evt,cut[0],None),cut[2]): 
                    survive=True 
                else: 
                    survive=False
                    #print "failed ",cut
            if statement=="absl":
                if (ops[">"](getattr(evt,cut[0],None),(-1 * cut[2])) and ops["<"](getattr(evt,cut[0],None),cut[2])): 
                    survive=True
                else: 
                    survive=False
                    #print "failed ",cut
            else:  
                if ops[statement](getattr(evt,cut[0],None),cut[2]):
                    survive=True
                    #print "passed cut",cut
                    continue 
                else:
                    survive=False
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

def createHistogram(cat,numvar,filelist,processObj,process,variable,weight,filedict):
    filedict[variable[0]].cd()
    filedict[variable[0]].mkdir(cat.name)
    filedict[variable[0]].cd(cat.name)

    bins = cat.binning[numvar]
    tmpbin = np.asarray(cat.binning[numvar])
   
    if filedict[variable[0]].Get(cat.name).GetListOfKeys().Contains(str(process)):
        tmp = ROOT.TH1D("tmp","tmp",len(tmpbin)-1,tmpbin)
        real = filedict[variable[0]].Get(cat.name).Get(str(process))
    else:
        tmp = ROOT.TH1D(str(process),str(process),len(tmpbin)-1,tmpbin)

    #Cuts
    procut = processObj.cuts[process]

    if procut=="":
        allcuts = "("+str(cat.cuts["categoryCuts"])+"&&"+str(cat.cuts["preselection"])+"&&"+str(cat.cuts["trigger"])+")"
        #print "Process doesn't have extra cut?" 
    else:
        allcuts = "("+str(cat.cuts["categoryCuts"])+"&&"+str(cat.cuts["preselection"])+"&&"+str(cat.cuts["trigger"])+"&&"+str(procut)+")"

    #Weights addition to common
    weightfinal = weight
    weightDict = processObj.weights
    for scalefactor in weightDict.keys():
        if scalefactor == "xsec":
            weightfinal =  weightfinal+"*"+str(weightDict[scalefactor])
        elif scalefactor == "nevents":
            weightfinal =  weightfinal+"*"+"("+"1/"+str(weightDict[scalefactor])+")"
        else: 
            weightfinal = weightfinal+"*"+str(weightDict[scalefactor])
    
    #capping off the weight for saftey in draw method
    weightfinal = "("+weightfinal+")"
    
    print process,"  the cuts and weight ",allcuts+"*"+weightfinal

    #tell of histogram already exits in the root output file - if so then add... 
    if filedict[variable[0]].Get(cat.name).GetListOfKeys().Contains(str(process)):
        tmp = ROOT.TH1D("tmp","tmp",len(tmpbin)-1,tmpbin)
        tree.Draw(str(variable[0])+">>tmp",allcuts+"*"+weightfinal)
        print "adding temp histo to process ... will gain ",tmp.GetEntries(),"   entries"
        real.Add(tmp)
        real.Write(real.GetName(),ROOT.TObject.kOverwrite)
        
    else:
        tmp = ROOT.TH1D(str(process),str(process),len(tmpbin)-1,tmpbin)
        tree.Draw(str(variable[0])+">>"+str(process),allcuts+"*"+weightfinal)
        tmp.Write(tmp.GetName(),ROOT.TObject.kOverwrite)

    return


def calculateHistos(tree,allcats,processObj,nickname,histodict,commonweight):

    #do I need to use branches?

    
    for ievt, evt in enumerate(tree):
        #copy the standard weights... will add scale factors to this later 

        #do stuffs with  tree  
        if ievt % 1000 == 0:
            print "processing event ",ievt
        
        #fill histograms
        #FillHistograms(ievt,evt,allcats,filedict)
        
        
        for cat in allcats:
            for variable in cat.variables:
                for process in processObj.cuts.keys():
                    #Cuts
                    procut = processObj.cuts[process]

                    #Weights addition to common
                    weightfinal = commonweight
                    weightDict = processObj.weights
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

                    #Operator expansion cut string
                    survive=cutStringBool(evt,cuts)
                    
                  
      
                    if survive==True:
                        #obtaining the right variable...
                        val = getattr(evt,variable[0],None)
                        filedict[variable[0]].cd()
                        filedict[variable[0]].cd(cat.name)
                        histodict[variable[0]+":"+cat.name+":"+process].Fill(float(val),float(weightfinal))


    return

if __name__ == "__main__":
    


    #Structure for plotting variables 
    from utils.Parametrization import Category
    from utils.Parametrization import Process


    #Structure for mapping between root files and processes  
    from utils.Processes import HAA
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


    #This is where the plotting takes place!
    #categories.append(HAA_Inc_mmmt) 
    #categories = allcats

    #Gather the Analysis Files
    dir = "/eos/user/s/shigginb/HAA_ntuples/March_2020/"
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
    for cat in allcats:
        numvar=0
        for variable in cat.variables:

            #filedict[variable[0]]=ROOT.TFile.Open(cat.name+"_"+str(variable[0])+".root","RECREATE")
            filedict[variable[0]]=ROOT.TFile.Open(str(variable[0])+".root","RECREATE")
            #print "on variable ",variable
            filedict[variable[0]].cd()
            filedict[variable[0]].mkdir(cat.name)
            filedict[variable[0]].cd(cat.name)

            for nickname in filelist.keys():

                processObj = HAA_processes[nickname]

                for process in processObj.cuts.keys():


                    if filedict[variable[0]].Get(cat.name).GetListOfKeys().Contains(str(process)):
                        continue
                    else:
                        bins = cat.binning[numvar]
                        tmpbin = np.asarray(cat.binning[numvar])
                        histodict[variable[0]+":"+cat.name+":"+process] = ROOT.TH1D(str(process),str(process),len(tmpbin)-1,tmpbin)
                        histodict[variable[0]+":"+cat.name+":"+process].Write(str(process),ROOT.TObject.kOverwrite)
            numvar=numvar+1

    #gather extra global variables or weights
    
    weight = CommonWeights["lumi"][0]
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

        calculateHistos(tree,allcats,processObj,nickname,histodict,weight)
        
     
        #write new ttree??
        
    

    #for cat in categories:
    #    print "working on category",cat.name
    #     
    #    ##parms = Params(category.varlist,category.cutlist) 
    #    #need to make a chain fo TFile or something... 
    #    filedict = {}

    #    #Creating output files
    #    for variables in cat.variables:
    #        if len(cat.variables)>1:
    #            filedict[variables[0]]=ROOT.TFile.Open(cat.name+"_"+str(variables[0])+".root","RECREATE")
    #        else:
    #            filedict[variables]=ROOT.TFile.Open(cat.name+"_"+str(variables)+".root","RECREATE")
    #
    #    #for file in filelist.keys():
    #    for nickname in filelist.keys():

    #        fin = ROOT.TFile.Open(dir+filelist[nickname],"read")
    #        tree = fin.Get(treename)
    #        print "tree entries ",tree.GetEntries()

    #        #processes will be a dictionary of key: process name value: cuts
    #        #processes = HAA_processes[nickname].cuts
    #        processObj = HAA_processes[nickname]

    #        tree = calculateHistos(tree,cat,processObj,nickname)

    #        #for process in filelist[file]:

    #        for process in processObj.cuts.keys():
    #            
    #            numvar=0
    #            weight = CommonWeights["lumi"][0]
    #            weight = weight+"*"+CommonWeights["mcweight"][0]
    #            #weight = "*("+weight+")"
    #            #weight = "*"+weight+"*"+SMHTTWeights["TES"][0]
    #            if process=="data_obs":
    #                weight = ""
    #             
    #            for variable in cat.variables:
    #                print "on variable ",variable
    #                #if len(variable)>1:
    #                    #createUnrolledHistogram(cat,numvar,filelist,processObj,variable,weight,filedict)
    #                #else:
    #                createHistogram(cat,numvar,filelist,processObj,process,variable,weight,filedict)
    #                numvar=numvar+1

    #for varcatPro in histodict.keys():
    #    print varcatPro," final entries ",histodict[varcatPro].GetEntries()
                
    print "writing the histograms to the file"
    for varcatPro in histodict.keys():
        #print "file ",filedict[varcatPro.split(":")[0]].GetName()," contents ",filedict[varcatPro.split(":")[0]].ls()
        print "file ",varcatPro.split(":")[0]," writing ",varcatPro,"  final entries ",histodict[varcatPro].GetEntries()
        filedict[varcatPro.split(":")[0]].cd()
        filedict[varcatPro.split(":")[0]].cd(varcatPro.split(":")[1])
        histodict[varcatPro].Write(histodict[varcatPro].GetName(),ROOT.TObject.kOverwrite)
            

    #for cat in allcats:
    #    numvar=0
    #    for variable in cat.variables:
    #        filedict[variable[0]].cd(cat.name)
    #        print "Writing ",variable[0]+":"+process,"  final entries ",histodict[variable[0]+":"+process].GetEntries()
    #        histodict[variable[0]+":"+process].Write(str(process),ROOT.TObject.kOverwrite)
    #                    
    #        numvar=numvar+1

