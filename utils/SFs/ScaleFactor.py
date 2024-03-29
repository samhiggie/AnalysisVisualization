       
import ROOT as root
import sys,math
from array import array

       
global etaLabel
       
class SFs():  
    def __init__(self):
        # initialize your global vars instead as
        global eff_dataH #= root.std.map("string", root.TGraphAsymmErrors)()
        global eff_mcH #= root.std.map("string", root.TGraphAsymmErrors)()
        global inputRootFile

    def ScaleFactor(self,inputFile) :
	self.inputRootFile = str(inputFile)
        self.eff_dataH = root.std.map("string", root.TGraphAsymmErrors)()   
        self.eff_mcH = root.std.map("string", root.TGraphAsymmErrors)()
        EtaBins=[]
        if "Mu" in str(inputFile) :   EtaBins=["Lt0p9", "0p9to1p2","1p2to2p1","Gt2p1"]
        if "El" in str(inputFile) :   EtaBins=["Lt1p0","1p0to1p48","1p48to1p65","1p65to2p1", "Gt2p1"]
       

        print 'Opening ScaleFactors file', inputFile
        self.fileIn = root.TFile(self.inputRootFile,"read") 
	print '====================================================='
        self.fileIn.ls()  
        HistoBaseName = "ZMassEta"
        self.etaBinsH = self.fileIn.Get("etaBinsH")
        nEtaBins = int(self.etaBinsH.GetNbinsX())
        #print "EtaBins...........",nEtaBins, len(EtaBins)
        for iBin in range (0, nEtaBins) : 
            etaLabel = EtaBins[iBin]
            GraphName = HistoBaseName+etaLabel+"_Data"
       	    #print GraphName,etaLabel
       
       	    self.eff_dataH[etaLabel]=self.fileIn.Get(str(GraphName))
            self.SetAxisBins(self.eff_dataH[etaLabel])
       	 
       	    GraphName = HistoBaseName+etaLabel+"_MC"
       	    self.eff_mcH[etaLabel]=self.fileIn.Get(str(GraphName))
            self.SetAxisBins(self.eff_mcH[etaLabel])
       	 
            #print "some checks.................",self.eff_mcH[etaLabel].GetXaxis().GetNbins(),self.eff_dataH[etaLabel].GetXaxis().GetNbins(),"etaLabel",etaLabel,self.eff_dataH[etaLabel].GetN(),"eff_mcH.GetN()",self.eff_mcH[etaLabel].GetN()
            #print "just get some value",self.eff_mcH[etaLabel].GetX()[5],"for etaLabel",etaLabel

        

    def SetAxisBins(self,graph):
        NPOINTS = graph.GetN() 
        AXISBINSS = array('d')
        for i in range(0, NPOINTS):  
            AXISBINSS.append (graph.GetX()[i] - graph.GetErrorXlow(i))

        AXISBINSS.append( graph.GetX()[NPOINTS-1] + graph.GetErrorXhigh(NPOINTS-1))

        graph.GetXaxis().Set(int(NPOINTS), AXISBINSS)

    def get_ScaleFactor(self,pt, eta) :
        efficiency_data = self.get_EfficiencyData(pt, eta)
        efficiency_mc = self.get_EfficiencyMC(pt, eta)
        if  efficiency_mc != 0. : 
            SF = float(efficiency_data)/float(efficiency_mc)
        else  : 	
            SF=1.
    
        #print "ScaleFactor::get_ScaleFactor(double pt, double eta) Scale Factor set to",SF,self.efficiency_data,self.efficiency_mc
        return SF	
    
    
    def get_EfficiencyMC(self,pt, eta) :
    
        label = str(self.FindEtaLabel(eta))
        binNumber = self.etaBinsH.GetXaxis().FindFixBin(eta)
        label = self.etaBinsH.GetXaxis().GetBinLabel(binNumber)
        ptbin = self.FindPtBin("mc", label, pt)
        Eta = math.fabs(eta)
        label=label.replace("Eta","")
        #print "inside eff_mc pt",pt,"eta",eta,"label",label,"N",self.eff_mcH[label].GetN(),"ptbin",ptbin
     
  
        if ptbin == -99 : eff =1
        else  : eff= self.eff_mcH[label].GetY()[ptbin-1]
    
        if eff > 1.  : eff = -1 
        if eff < 0 : eff = 0. 
        return eff
    
    def get_EfficiencyData(self,pt, eta) :
    
        label = self.FindEtaLabel(eta)
        binNumber = self.etaBinsH.GetXaxis().FindFixBin(eta)
        label = self.etaBinsH.GetXaxis().GetBinLabel(binNumber)
        #print self.eff_dataH
        ptbin = self.FindPtBin("data", label, pt)
        Eta = math.fabs(eta)
        label=label.replace("Eta","")
        #print "inside eff_data pt",pt,"eta",eta,"label",label,"N",self.eff_dataH[label].GetN(),"ptbin",ptbin
    
        if ptbin == -99 : eff =1
        else  : eff= self.eff_dataH[label].GetY()[ptbin-1]
        #print "inside eff_data",eff
    
        if eff > 1.  : eff = -1 
        if eff < 0 : eff = 0. 
    
        return eff

    def FindPtBin(self, whatisit,EtaLabel, Pt) :
        
        eff_map=self.eff_mcH
        EtaLabel=EtaLabel.replace("Eta","")
        if whatisit == "mc" : eff_map=self.eff_mcH #EtaLabel = "ZMass"+EtaLabel+"_MC"
        if whatisit == "data" : eff_map=self.eff_dataH #EtaLabel = "ZMass"+EtaLabel+"_MC"
        #if whatisit == "data" : EtaLabel = "ZMass"+EtaLabel+"_Data"
        #print "and inside FintPtBin",EtaLabel,self.eff_mcH[EtaLabel].GetN(),EtaLabel,eff_map[EtaLabel].GetN()

        Npoints = eff_map[EtaLabel].GetN()
        ptMAX, ptMIN = 0,0
        try : ptMAX = (eff_map[EtaLabel].GetX()[Npoints-1])+(eff_map[EtaLabel].GetErrorXhigh(Npoints-1))
        except IndexError : return -99

        try : ptMIN = (eff_map[EtaLabel].GetX()[0])-(eff_map[EtaLabel].GetErrorXlow(0))
        except IndexError : return -99
        #print Npoints, "Npoints for  ===============>",eff_map[EtaLabel].GetN(),EtaLabel,Pt,eff_map[EtaLabel].GetN(),whatisit,ptMAX,ptMIN,eff_map[EtaLabel].GetXaxis().FindBin(Pt)
        if Pt >= ptMAX : return Npoints
        elif Pt < ptMIN :
            return -99
        else : return eff_map[EtaLabel].GetXaxis().FindFixBin(Pt)
    
    
    def FindEtaLabel(self,Eta) : 
         
        Eta = math.fabs(Eta)
        binNumber = self.etaBinsH.GetXaxis().FindFixBin(Eta)
        EtaLabel = self.etaBinsH.GetXaxis().GetBinLabel(binNumber)
    	
        #print "inside FindEtaLabel",EtaLabel
    	
        return EtaLabel

#sf = SFs()
#sf.ScaleFactor("Muon_IsoMu27.root")
