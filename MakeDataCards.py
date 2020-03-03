#########################
#Author: Sam Higginbotham
'''

* File Name : MakeDataCards.py

* Purpose : For Datacard creation. Root file containing histograms that can be used with Combine.

* Creation Date : 04-02-2020

* Last Modified :

'''
#########################
import ROOT
import numpy as np
#import argsparse
from array import array
import itertools

#Structure for plotting variables 
from utils.Parametrization import Params


#Structure for mapping between root files and processes  
from utils.Processes import SMHTT
from utils.ProcessCuts import ProCuts

#Structure for mapping between processes and weights
from utils.Weights import CommonWeights
from utils.Weights import SMHTTWeights

#Structure for mapping processes and extra distributions 
#from utils.FunctionDefinitions.Systematics import systematics

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


if __name__ == "__main__":

    #parser = argsparse.ArgumentParser(description="make datacards")
    #parser.add_arguement('--CategoryFiles',nargs="+",help="Select the files containing the categories for the datacards")
    #args = parser.parse_args()
    #loader = 
    categories = []
    #gather all the analysis categories
    from utils.Categories import Inclusive
    from utils.Categories import vbftest

    #categories.append(Inclusive) 
    categories.append(vbftest) 

    #Gather the Analysis Files
    dir = "reduce/"
    filelist = SMHTT
    treename = "muTauEventTree/eventTree"

    for cat in categories:
        print "working on category",cat.name
         
        ##parms = Params(category.varlist,category.cutlist) 
        #need to make a chain fo TFile or something... 
        filedict = {}

        for variables in cat.variables:
            if len(cat.variables)>1:
                filedict[variables[0]]=ROOT.TFile.Open(cat.name+"_"+str(variables[0])+".root","RECREATE")
            else:
                filedict[variables]=ROOT.TFile.Open(cat.name+"_"+str(variables)+".root","RECREATE")
    
        for file in filelist.keys():

            fin = ROOT.TFile.Open(dir+file,"read")
            tree = fin.Get(treename)
            print "tree entries ",tree.GetEntries()


            for process in filelist[file]:
                print "working on ",process,"    in file   ",file
                numvar=0
                innervar=0
                weight = CommonWeights["lumi"][0]
                weight = "*"+weight+"*"+CommonWeights["weights"][0]
                if process=="data_obs":
                    weight = ""
                for variables in cat.variables:
                    if len(cat.variables)>1:
                        filedict[variables[0]].cd()
                        filedict[variables[0]].mkdir(cat.name)
                        filedict[variables[0]].cd(cat.name)
                    else:
                        filedict[variables].cd()
                        filedict[variables].mkdir(cat.name)
                        filedict[variables].cd(cat.name)
                    print "working on variable", variables," numvar ",numvar

                    if len(variables)>1:


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
                        procut = ProCuts.get(process)
                        histList = []
                        urnum = 0 
                        for bincuts in combos:
                            htmp=ROOT.TH1D("htemp"+str(process)+str(urnum),"htemp"+str(process)+str(urnum),len(primaryBins)-1,primaryBins)
                            print "working on variable   ",str(variables[0]),"for htemp"+str(process)+str(urnum)
                            unrollCut=""

                            for extraVarnum in range(1,len(variables)):
                                unrollCut += "&&"+str(variables[extraVarnum])+">"+str(bincuts[extraVarnum-1][0])+"&&"+str(variables[extraVarnum])+"<"+str(bincuts[extraVarnum-1][1])
                            if procut:
                                allcuts = str(cat.cuts)+"&&"+str(ProCuts["preselection"])+"&&"+str(ProCuts["trigger"])+"&&"+str(ProCuts[process])+unrollCut
                                #print "entries after cut  ",tree.GetEntries(allcuts)
                                tree.Draw(str(variables[0])+">>"+"htemp"+str(process)+str(urnum),allcuts+weight)
                            else:
                                allcuts = str(cat.cuts)+"&&"+str(ProCuts["preselection"])+"&&"+str(ProCuts["trigger"])+unrollCut
                                #print "entries after cut  ",tree.GetEntries(allcuts)
                                tree.Draw(str(variables[0])+">>"+"htemp"+str(process)+str(urnum),allcuts+weight)
                             
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

                    else:

                        bins = cat.binning[numvar]
                        print "working on single variable", variables 
                        #tmpbin = np.array('d',cat.binning)
                        tmpbin = np.asarray(cat.binning)
                        print bins
                        tmp = ROOT.TH1D(str(process),str(process),len(tmpbin)-1,tmpbin)
                        procut = ProCuts.get(process)
                        if procut:
                            allcuts = str(cat.cuts)+"&&"+str(ProCuts["preselection"])+"&&"+str(ProCuts["trigger"])+"&&"+str(ProCuts[process])
                        else:
                            allcuts = str(cat.cuts)+"&&"+str(ProCuts["preselection"])+"&&"+str(ProCuts["trigger"])
                            print "Process doesn't have extra cut?" 
                        print "The cuts ", allcuts
                        tree.Draw(str(variables[0])+">>"+str(process),allcuts+weight)
                        tmp.Write(tmp.GetName(),ROOT.TObject.kOverwrite)
                    numvar=numvar+1
