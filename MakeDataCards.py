#########################
#Author: Sam Higginbotham
'''

* File Name : MakeDataCards.py

* Purpose : For Datacard creation. Root file containing histograms that can be used with Combine. Motivation working on HTT with Andrew Loeliger

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
    treename = "mutauEventTree/eventTree"

    for cat in categories:
        print "working on category",cat.name

         
        ##parms = Params(category.varlist,category.cutlist) 
    
        for file in filelist.keys():

            fin = ROOT.TFile.Open(dir+file,"read")
            tree = fin.Get(treename)

            for process in filelist[file]:
                print "working on ",process,"    in file   ",file
                numvar=0
                innervar=0
                for variables in cat.variables:
                    print "working on variable", variables," numvar ",numvar
                    if len(variables)>1:
                        #creating primary bins
                        primaryBins = cat.binning[numvar][0] #always the first binning set 
                        #creating secondary bin list
                        secondaryBins = []

                        for bins in range(1,len(cat.binning[numvar])):
                            secondaryBins.append(cat.binning[numvar][bins])
                    
                        binlist = returnBinList(secondaryBins)
                        combos = returnPerms(binlist)
                        #print "Primary Bins  ",primaryBins
                        #print "Combinations for secondary bins  ",combos
                        procut = ProCuts.get(process)
                        
                        for bincuts in combos:
                            unrollCut=""
                            for extraVarnum in range(1,len(variables)):
                                unrollCut += "&&"+str(variables[extraVarnum])+">"+str(bincuts[extraVarnum-1][0])+"&&"+str(variables[extraVarnum])+"<"+str(bincuts[extraVarnum-1][1])
                            if procut:
                                allcuts = str(cat.cuts)+"&&"+str(ProCuts["preselection"])+"&&"+str(ProCuts["trigger"])+"&&"+str(ProCuts[process])+unrollCut
                            else:
                                allcuts = str(cat.cuts)+"&&"+str(ProCuts["preselection"])+"&&"+str(ProCuts["trigger"])+unrollCut
                            print allcuts
                        

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
                    numvar=numvar+1

                    #tree.Draw(str(var)+">>"+str(process),allcuts)
                    #fileout = ROOT.TFile.Open(cat.name+"_"+str(var)+".root")
                    #fileout.cd()
                    #fileout.mkdir(cat.name)
                    #fileout.cd(cat.name)
                    #fileout.Write(tmp,ROOT.TOverWrite())
