#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <math.h>
#include <TRandom2.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TROOT.h>
#include <bits/stdc++.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <math.h>
#include <TRandom2.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TLegend.h>

float SD(std::vector <double> vec){

  double sum=0.;
  double standev=0.;
  double mean=0.;

  for (size_t ii=0;ii<vec.size();ii++){
    sum+=vec[ii];
  }

  mean=sum/vec.size();

  for (size_t jj=0;jj<vec.size();jj++){

    standev+=pow(vec[jj]-mean,2);
    
  }

  return sqrt(standev/vec.size());
  
}

std::pair< std::pair<double,double> , std::pair<double,double> > ave_points(std::vector<double> volt,std::vector<double> time){//first pair voltage data and error, second pair is time data and error

  double ttot;//total time
  double tave;//time average
  double terr;//time error
  double vave;//storing average for voltage
  double verr;//non-time error
  

   for (size_t ii=0;ii<volt.size();ii++){

    ttot+=time[ii];//adding up time
    vave=vave+volt[ii];

   }
  
  
  tave = ttot/time.size();//time average
  terr = (time[time.size()-1] - time[0])/2.0;
  vave = vave/25.0; //std::accumulate(volt.begin(),volt.end(),0);
  verr = SD(volt);

  //vave = hist->GetMean();//average of non time data
  //verr = hist->GetRMS();

  std::pair< std::pair<double,double> , std::pair<double,double> > ret;

  ret.first.first = vave;
  ret.first.second = verr;
  ret.second.first = tave;
  ret.second.second = terr;

  ttot=0;

  return ret;

}

int main(int argc,char* argv[]){

  TGraphErrors* graph = new TGraphErrors();

  TGraph* raw_graph = new TGraph(argv[1],"%lg,%lg");

  int counter = 0;//initialise counter for averaging data

  std::pair< std::pair<double,double> , std::pair<double,double> > aves;

  std::vector<double> tvec;//initialise vector for storing average time
  std::vector<double> tvec_err;//initialise vector for storing time error
  std::vector<double> vvec;//initialise vector for storing average voltage
  std::vector<double> vvec_err;//initialise vector for storing voltage error

  std::vector<double> v_calc;//initialise vector for calculating voltage average and error
  std::vector<double> t_calc;//initialise vector for calculating time average and error

  for (int ii=0;ii<raw_graph->GetN();ii++){

    v_calc.push_back(raw_graph->GetY()[ii]);
    t_calc.push_back(raw_graph->GetX()[ii]);

    counter=counter+1;//increment counter

    if (counter == 25){

                     aves = ave_points(v_calc,t_calc);
		     
		     vvec.push_back(aves.first.first);//add average and errors to end of vectors
		     vvec_err.push_back(aves.first.second);
		     tvec.push_back(aves.second.first);
		     tvec_err.push_back(aves.second.second);

		     counter=0;//reset counter
		     
		     v_calc.clear();
		     t_calc.clear();
		     
      }
    
  }

  for (size_t kk=0;kk<vvec.size();kk++){

    graph->SetPoint(kk,tvec[kk],vvec[kk]);
    graph->SetPointError(kk,tvec_err[kk],vvec_err[kk]);

  }

  ///////EVAL BOARD 1 DATA/////////

  auto c1 = new TCanvas("c1","test",1000,1000);
  //TF1 *f1 = new TF1("f1","[0]*x+[1]",-5,5);
  //graph->Fit(f1);
  //double ndf=f1->GetNDF();
  //double chi2=f1->GetChisquare();
  //double r_chi2 = chi2/ndf;
  graph->SetMarkerStyle(21);
  graph->SetMarkerColor(1);
  graph->SetTitle("Oscilloscope Readout;Time (s);Voltage (V)");
  graph->Draw("AP");
  c1->SaveAs("test.pdf");
  //std::cout<<"Reduced Chi2 is: ";
  //std::cout<<r_chi2;
  //std::cout<<"\n";
  
  
  
    

  return 0;
}

  
  