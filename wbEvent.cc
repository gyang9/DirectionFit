#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>

#include "wbEvent.hh"
#include "TMath.h"
#include "RooArgList.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "TString.h"
#include <TVectorD.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFrame.h>
#include <TFile.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1D.h>
#include <TVector3.h>

using namespace std;

wbEvent::wbEvent (const char* name)
{
  
}

std::vector<wbEvent::wbHit> wbEvent::sample_photons( std::vector<wbHit> list, bool smearVtx , bool smearTime , bool smearQE , bool addNoise  ){ 
  // Ben Land's code at the end of the script
  cout<<"currently there is no smearing applied .. "<<endl;
  for (wbHit i: list){
    if(smearVtx){}
    if(smearTime){}
    if(smearQE){}
    if(addNoise){}
  } 
  return list;
}

wbPDF* wbEvent::create2DPDFs(std::vector<wbHit> list, double prompt_cut, double wavelength_cut, bool doCharge, bool doCos, int nbins){
  
  double c = 300/1.333; //mm/ns
  
  int count = 0;
  wbPDF* wbpdf = new wbPDF("_wbpdf");
  wbpdf->SetTimePDFBinning(30,-10,20);
  if (!doCos)
    wbpdf->SetThetaPDFBinning(nbins,-3.14,3.14);
  else 
    wbpdf->SetThetaPDFBinning(nbins,-1,1);

  cout<<"creating 2D pdf.. total wbhit number "<<list.size()<<endl;
  for (wbHit i: list){
    //cout<<"prompt cut "<<prompt_cut<<endl;
    if (i.t > prompt_cut) continue;
    double x     = i.trX;
    double y     = i.trY;
    double z     = i.trZ;
    double t     = i.trTime;
    double theta = i.trTheta;
    double phi   = i.trPhi;
    double charge= i.charge;

    TVector3 P (i.x - x, i.y - y, i.z - z);
    double D = TMath::Sqrt(TMath::Power(i.x - x, 2) + TMath::Power(i.y - y, 2) + TMath::Power(i.z - z, 2));
    double T = i.t - t;
    double tresid = T;
    //double tresid = T - D/c;
    TVector3 dvec (TMath::Cos(phi)*TMath::Sin(theta), TMath::Sin(phi)*TMath::Sin(theta), TMath::Cos(theta));
    double costheta;
    if (doCos)
      costheta = TMath::Cos(dvec.Angle(P));
    else
      costheta = dvec.Angle(P);

    double no_wavelength_cut = wavelength_cut;
    if (doCharge){
      wbpdf->GetTimePDF()->TH1F::Fill(tresid, i.charge);
      wbpdf->GetThetaPDF()->TH1F::Fill(costheta, i.charge);
      wbpdf->GetTimeThetaPDF()->TH2F::Fill(tresid, costheta, i.charge);
    }
    else {
      wbpdf->GetTimePDF()->TH1F::Fill(tresid);
      wbpdf->GetThetaPDF()->TH1F::Fill(costheta);
      wbpdf->GetTimeThetaPDF()->TH2F::Fill(tresid, costheta);
    }
    count ++;
  }

  return wbpdf;
}

wbPDF* wbEvent::createPDFs(std::vector<wbHit> list, double prompt_cut, double wavelength_cut, bool doCharge, bool doCos, int nbins){
  // Ben Land's code at the end of the script
  double c = 300/1.333; //mm/ns

  int count = 0;
  wbPDF* wbpdf = new wbPDF("_wbpdf");
  wbpdf->SetTimePDFBinning(30,-10,20);
  if (!doCos)
    wbpdf->SetThetaPDFBinning(nbins,-3.14,3.14);
  else 
    wbpdf->SetThetaPDFBinning(nbins,-1,1);
  cout<<"creating pdf.. total wbhit number "<<list.size()<<endl;
  for (wbHit i: list){
    //cout<<"prompt cut "<<prompt_cut<<endl;
    if (i.t > prompt_cut) continue;
/*	  
    double x     = tPara[count].at(0);
    double y     = tPara[count].at(1);
    double z     = tPara[count].at(2);
    double t     = tPara[count].at(3);
    double theta = tPara[count].at(4);
    double phi   = tPara[count].at(5);
    count++;
*/
    double x     = i.trX;
    double y     = i.trY;
    double z     = i.trZ;
    double t     = i.trTime;
    double theta = i.trTheta;
    double phi   = i.trPhi;     

    //cout<<"element "<<count<<"   x y z t theta phi "<<x<<" "<<y<<" "<<z<<" "<<theta<<" "<<phi<<endl;
    //cout<<"i.x  i.y  i.z  i.t "<<i.x<<" "<<i.y<<" "<<i.z<<" "<<i.t<<endl;

    TVector3 P (i.x - x, i.y - y, i.z - z);
    double D = TMath::Sqrt(TMath::Power(i.x - x, 2) + TMath::Power(i.y - y, 2) + TMath::Power(i.z - z, 2));
    double T = i.t - t;
    //double tresid = T - D/c;
    double tresid = T;
    TVector3 dvec (TMath::Cos(phi)*TMath::Sin(theta), TMath::Sin(phi)*TMath::Sin(theta), TMath::Cos(theta));
    
    double costheta;
    if (doCos)
      costheta= TMath::Cos(dvec.Angle(P));
    else
      costheta = dvec.Angle(P);

    double no_wavelength_cut = wavelength_cut;

    //wbpdf->SetTimePDFBinning(20,-10,10);
    //wbpdf->SetThetaPDFBinning(48,-6,6);
    if (doCharge){
      wbpdf->GetTimePDF()->TH1F::Fill(tresid, i.charge);
      wbpdf->GetThetaPDF()->TH1F::Fill(costheta, i.charge);    
    }
    else{
      wbpdf->GetTimePDF()->TH1F::Fill(tresid);
      wbpdf->GetThetaPDF()->TH1F::Fill(costheta);
    }
    count ++;
/*	    
    if coordinate is not None: #coordinate variable should be the truth information
        x,y,z,t,theta,phi = coordinate
        pos = np.asarray([x,y,z])
        P = xyz - pos
        D = np.sqrt(np.sum(np.square(P),axis=1))
        T = times - t
        tresid = T - D/c
        dvec = np.asarray([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])
        costheta = np.sum(P*dvec,axis=1)/D
        if prompt_cut is None:
            return tresid,costheta
        else:
            mask = tresid<prompt_cut
            return tresid,costheta[mask]	  
*/  
  
  }

  return wbpdf;
}

std::vector<wbPDF*> wbEvent::createPMTPDFs(std::vector<wbHit> list, double prompt_cut, double wavelength_cut, bool doCharge, bool doCos, int nbins){
  // Ben Land's code at the end of the script
  double c = 300/1.333; //mm/ns

  int count = 0;

  std::vector<wbPDF*> wbpdf(500); 
  for (int ipmt = 0; ipmt< 500;ipmt++){
    wbpdf[ipmt] = new wbPDF("_wbpdf");
    wbpdf[ipmt]->SetPMTPDFBinning(nbins,0,3.14,nbins,0,6.28);
  }
  cout<<"creating pdf.. total wbhit number "<<list.size()<<endl;
  for (wbHit i: list){
    //cout<<"prompt cut "<<prompt_cut<<endl;
    if (i.t > prompt_cut) continue;
    double x     = i.trX;
    double y     = i.trY;
    double z     = i.trZ;
    int pmtid = i.pmtid;
    double t     = i.trTime;
    double theta = i.trTheta;
    double phi   = i.trPhi;

    //cout<<"element "<<count<<" out of "<<list.size()<<" pmt id "<<pmtid<<" theta and phi "<<theta<<" "<<phi<<endl;
    //cout<<"element "<<count<<"   x y z t theta phi "<<x<<" "<<y<<" "<<z<<" "<<theta<<" "<<phi<<endl;
    //cout<<"i.x  i.y  i.z  i.t "<<i.x<<" "<<i.y<<" "<<i.z<<" "<<i.t<<endl;

    TVector3 P (i.x - x, i.y - y, i.z - z);
    double D = TMath::Sqrt(TMath::Power(i.x - x, 2) + TMath::Power(i.y - y, 2) + TMath::Power(i.z - z, 2));
    double T = i.t - t;
    //double tresid = T - D/c;
    double tresid = T;
    TVector3 dvec (TMath::Cos(phi)*TMath::Sin(theta), TMath::Sin(phi)*TMath::Sin(theta), TMath::Cos(theta));

    if (doCos)
      int nothing = 1;
    else
      int nothing = 1;

    double no_wavelength_cut = wavelength_cut;

    if (doCharge){
      wbpdf[pmtid]->GetPMTPDF()->TH2F::Fill(theta, phi, i.charge);
    }
    else{
      double geoWei = 2*TMath::Pi()* (1-TMath::Cos(abs(theta)) ) - 2*TMath::Pi()* (1-TMath::Cos(abs(theta)-0.01) ); 
      wbpdf[pmtid]->GetPMTPDF()->TH2F::Fill(theta, phi, 1./geoWei);
    }
    count ++;
  }
  return wbpdf;
}

void wbEvent::SetHitList(std::vector<std::vector<double> > num){
  cout<<"total hit size "<<num.size()<<endl;
  for (int i=0;i<num.size();i++){
    //for (int j=0;j<10;j++){
    //  cout<<"element "<<j<<"  value  "<<num[i].at(j)<<endl;
    //}
    wbEvent::wbHit hit;
    hit.x = num[i].at(0);
    hit.y = num[i].at(1);
    hit.z = num[i].at(2);
    hit.t = num[i].at(3);
    hit.trX = num[i].at(4);
    hit.trY = num[i].at(5);
    hit.trZ = num[i].at(6);
    hit.trTime = num[i].at(7);
    hit.trTheta = num[i].at(8);
    hit.trPhi = num[i].at(9);
    hit.charge = num[i].at(10);
    _hitlist.push_back(hit);
  }
}

void wbEvent::SetHitListPMT(std::vector<std::vector<double> > num){
  cout<<"total hit size "<<num.size()<<endl;
  for (int i=0;i<num.size();i++){
    //for (int j=0;j<10;j++){
    //  cout<<"element "<<j<<"  value  "<<num[i].at(j)<<endl;
    //}
    wbEvent::wbHit hit;
    hit.x = num[i].at(0);
    hit.y = num[i].at(1);
    hit.pmtid = num[i].at(2);
    hit.t = num[i].at(3);
    hit.trX = num[i].at(4);
    hit.trY = num[i].at(5);
    hit.trZ = num[i].at(6);
    hit.trTime = num[i].at(7);
    hit.trTheta = num[i].at(8);
    hit.trPhi = num[i].at(9);
    hit.charge = num[i].at(10);
    _hitlist.push_back(hit);
  }
}

wbEvent ::~wbEvent ()
{;}

  
