#########################
#Author: Sam Higginbotham
'''

* File Name : Parametrization.py

* Purpose : Loads the variables that will be plotted against one another for further analysis  

* Creation Date : 04-02-2020

* Last Modified : Development ongoing - plot the same variable in multiple ways ... "extravariabels"

'''
#########################

class Params():
    def __init__(self):
        self.name = []
        self.variables = []
        self.cuts = {}
        self.binning = []
        self.extraplots = {}
        
class Process():
    def __init__(self):
        self.name = ""
        self.file = ""
        self.classification = ""
