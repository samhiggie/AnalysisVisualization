#########################
#Author: Sam Higginbotham
'''

* File Name :

* Purpose :

* Creation Date : 05-02-2020

* Last Modified :

'''
#########################
import ROOT
from AnalysisVisualization.Parametrization import Params

vbftest = Params()
vbftest.name = "vbf"
vbftest.variables = ["m_sv","mjj"]
vbftest.cuts = ["mt>50.0","againstMuonLoose_1>0.5"]   #stuff like this
vbftest.binning = [[50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0],[0.0,300.0,500.0]]
