#########################
#Author: Sam Higginbotham
'''

* File Name :Categories.py 

* Purpose : The Categories in the analysis... typically a set of extra cuts that define a phase space region of interest

* Creation Date : 05-02-2020

* Last Modified :

'''
#########################
import ROOT
from Parametrization import Params

Inclusive = Params()
Inclusive.name = "inclusive"
Inclusive.variables = [["m_sv"]]
Inclusive.cuts = "mt_1<50.0"   #stuff like this
Inclusive.binning = [50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0]


vbftest = Params()
vbftest.name = "vbf"
vbftest.variables = [["m_sv","mjj"],["njets"]]
vbftest.cuts = "mt_1>50.0&&njets>=2"   #stuff like this
vbftest.binning = [[[50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0],[0.0,50.0,100.0,200.0,300.0,500.0]],[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0]]


vbftestmela = Params()
vbftestmela.name = "vbfmela"
vbftestmela.variables = [["m_sv","mjj","mela"],["njets"]]
vbftestmela.cuts = "mt_1>50.0&&njets>=2"   #stuff like this
vbftestmela.binning = [[[50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0],[0.0,50.0,100.0,200.0,300.0,500.0],[0.0,4.0,5.0]],[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0]]

