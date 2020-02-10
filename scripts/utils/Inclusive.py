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

Inclusive = Params()
Inclusive.name = "inclusive"
Inclusive.variables = ["m_sv"]
Inclusive.cuts = ["mt>50.0","againstMuonLoose_1>0.5"]   #stuff like this
Inclusive.binning = [50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0]




