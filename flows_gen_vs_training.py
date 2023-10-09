import ROOT
from scipy import stats
import numpy
import sys
from math import * 

# define the output file
out_file = ROOT.TFile('flows_histograms.root', 'recreate')

#add data from input ROOT file to chain
chain = ROOT.TChain('tree')
chain.Add(sys.argv[1]) 

#def dicts associating label indices to data
ksvals={}
pmasses={}
omasses={}

#find max label index
nindex=int(chain.GetMaximum('labelIndex'))


for n in range(0,nindex+1): 

 genphi=[]
 genomega=[]
 trainphi=[]
 trainomega=[]
 labels=[]
  
 for event in chain:
  if event.labelIndex[0]==n:
   labels.append(list(event.label))
   if event.isGen[0]==0:
    trainphi.append(event.data[0])
    trainomega.append(event.data[1])
   elif event.isGen[0]==1:
    genphi.append(event.data[0])
    genomega.append(event.data[1]) 
    
 phiall=genphi+trainphi
 omegaall=genomega+trainomega  
 phiupper= round(1.1 * max(phiall) ,2)
 nbinphi=round(sqrt(len(phiall)),0)   
 omegaupper= round(1.1 * max(omegaall),3)
 nbinomega=round(sqrt(len(omegaall)),0)
 
 genphiarr=numpy.array(genphi)
 trainphiarr=numpy.array(trainphi) 
 genomegaarr=numpy.array(genomega) 
 trainomegaarr=numpy.array(trainomega)
  
 phikstest = stats.kstest(genphiarr,trainphiarr)
 omegakstest = stats.kstest(genomegaarr, trainomegaarr)

 ksvals.update({n:[phikstest.pvalue,omegakstest.pvalue,phikstest.statistic,omegakstest.statistic]}) 
  
 histogenphi=ROOT.TH1D('gen_phi_histo_' + str(n),'Generated Phi mass, ' + str(round(labels[0][0],2)) + ' GeV',int(nbinphi),0,phiupper)
 histotrainphi=ROOT.TH1D('train_phi_histo_' + str(n),'Training Phi mass, ' +str(round(labels[0][0],2))+' GeV',int(nbinphi),0,phiupper)

 histogenomega=ROOT.TH1D('gen_omega_histo_' + str(n),'Generated omega mass, ' + str(round(labels[0][1],3)) + ' GeV',int(nbinomega),0,omegaupper)
 histotrainomega=ROOT.TH1D('train_omega_histo_' + str(n),'Training omega mass, '+ str(round(labels[0][1],3))+' GeV',int(nbinomega),0,omegaupper) 

 pmasses.update({n:str(round(labels[0][0],2))})
 omasses.update({n:str(round(labels[0][1],3))})
 
 for a in genphi:
  histogenphi.Fill(a) 
 
 for b in genomega:
  histogenomega.Fill(b)
  
 for c in trainphi:
  histotrainphi.Fill(c)

 for d in trainomega:
  histotrainomega.Fill(d)
  
 out_file.Write()
 
#define modified chi-square test 
def chisquare(gen,train):
 chilist=[]
 num=gen.GetNbinsX()
 gen.Sumw2(ROOT.kFALSE)
 gen.SetBinErrorOption(ROOT.TH1.kPoisson)
 train.Sumw2(ROOT.kFALSE)
 train.SetBinErrorOption(ROOT.TH1.kPoisson)  
  
 for m in range(0,num):
  gerrlow=gen.GetBinErrorLow(m)
  gerrup=gen.GetBinErrorUp(m)
  terrlow=train.GetBinErrorLow(m)
  terrup=train.GetBinErrorUp(m)  
  gen.SetMarkerStyle(m)
  train.SetMarkerStyle(m)
  
  if train.GetBinContent(m)>3:
   if gen.GetBinContent(m) > train.GetBinContent(m):
    chilist.append(pow(gen.GetBinContent(m)-train.GetBinContent(m),2)/(pow(terrup,2) +pow(gerrlow,2)))
   elif train.GetBinContent(m) > gen.GetBinContent(m):
    chilist.append(pow(gen.GetBinContent(m)-train.GetBinContent(m),2)/(pow(terrlow,2) +pow(gerrup,2)))
   elif gen.GetBinContent(m) == train.GetBinContent(m):
    terrs=[]
    gerrs=[]
    terrs.append(terrlow)
    terrs.append(terrup)
    gerrs.append(gerrlow)
    gerrs.append(gerrup)
    tmax=max(terrs)
    gmax=max(gerrs)
    chilist.append(pow(gen.GetBinContent(m)-train.GetBinContent(m),2)/(pow(tmax,2) +pow(gmax,2)))
  else:
   pass
   
 chisquared=sum(chilist)
 return chisquared
 
 
for o in range(0,nindex+1):
 phigen=out_file.Get('gen_phi_histo_' +str(o))
 phitrain=out_file.Get('train_phi_histo_' +str(o))
 omegagen=out_file.Get('gen_omega_histo_' +str(o))
 omegatrain=out_file.Get('train_omega_histo_' +str(o))
 

 chi2phi=chisquare(phigen,phitrain)
 canvas1=ROOT.TCanvas('phi_plots_'+str(o),'phi_plots_'+str(o),700,500)
 xmax = 0.9
 ymax = 0.9
 xmin = 0.7
 ymin = 0.7
 leg = ROOT.TLegend(xmax, ymax, xmin, ymin)
 leg.AddEntry(phigen,str(o), 'L')
 leg.AddEntry('Chi Squared','Chi-squared test: '+str(round(chi2phi,2)),'L')
 leg.AddEntry('KS test','KS test p-value: ' + str(round(ksvals[o][0],7)),'L')
 leg.AddEntry('KS test','KS test statistic: ' + str(round(ksvals[o][2],7)),'L')
 leg.Draw()
 pgmax=phigen.GetBinContent(phigen.GetMaximumBin())
 ptmax=phitrain.GetBinContent(phitrain.GetMaximumBin())
 histmaxphi=max([pgmax,ptmax])
 phigen.SetMaximum(histmaxphi*1.2)
 phitrain.SetMaximum(histmaxphi*1.2)
 
 phigen.SetTitle('Phi mass, ' + pmasses[o] + ' GeV')
 phigen.Draw('') # 'E'
 phitrain.SetLineColor(ROOT.kRed)
 phitrain.Draw('same')
 leg.Draw('same')
 canvas1.Write('phi_plots_'+str(o),0,0)
 
 
 chi2omega=chisquare(omegagen,omegatrain)
 leg2 = ROOT.TLegend(xmax, ymax, xmin, ymin)
 leg2.AddEntry(omegagen,str(o), 'L')
 leg2.AddEntry('Chi Squared','Chi-squared test: '+str(round(chi2omega,2)),'L')
 leg2.AddEntry('KS test','KS test p-value: ' + str(round(ksvals[o][1],7)),'L')
 leg2.AddEntry('KS test','KS test statistic: ' + str(round(ksvals[o][3],7)),'L')
 ogmax=omegagen.GetBinContent(omegagen.GetMaximumBin())
 otmax=omegatrain.GetBinContent(omegatrain.GetMaximumBin())
 histmaxomega=max([ogmax,otmax])
 omegagen.SetMaximum(histmaxomega*1.2)
 omegatrain.SetMaximum(histmaxomega*1.2)
 
 canvas2=ROOT.TCanvas('omega_plots_'+str(o),'omega_plots_'+str(o),700,500)
 omegagen.SetTitle('Omega mass, ' + omasses[o] + ' GeV')
 omegagen.Draw('')
 omegatrain.SetLineColor(ROOT.kRed)
 omegatrain.Draw('same')
 leg2.Draw('same')
 canvas2.Write('omega_plots_'+str(o),0,0)
 
out_file.cd()
out_file.Close()
