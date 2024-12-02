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

std::vector<double> sq_fit(TGraphErrors* graph_to_fit, int charge){

  TF1 *f1 = new TF1("f1","[0]",-5E-6,-1E-6);//fit for the square wave
  TF1 *f2 = new TF1("f2","[0]",1E-6,5E-6);//fit for the square wave
  //fit to graph
  graph_to_fit->Fit(f1,"R");
  graph_to_fit->Fit(f2,"R");
  //Chi2 of both fits
  double ndf_1=f1->GetNDF();
  double chi2_1=f1->GetChisquare();
  double r_chi2_1 = chi2_1/ndf_1;
  
  double ndf_2=f2->GetNDF();
  double chi2_2=f2->GetChisquare();
  double r_chi2_2 = chi2_2/ndf_2;

  std::vector<double> out;
  
  //print out chi2 and values of fit
  std::cout<<"The value of f1 is "<<f1->GetParameter(0)<<" +/- "<<f1->GetParError(0)<<" and the chi2 is "<<r_chi2_1<<"\n";
  std::cout<<"The value of f2 is "<<f2->GetParameter(0)<<" +/-" <<f2->GetParError(0)<<" and the chi2 is "<<r_chi2_2<<"\n";

  if (charge==0){

  //calculate signal voltage
  double vlt = abs((f1->GetParameter(0)-f2->GetParameter(0)));
  double vlt_err = sqrt(pow(f1->GetParError(0),2)+pow(f2->GetParError(0),2));

  out.push_back(vlt);
  out.push_back(vlt_err);

  }


  if (charge==1){

    double vlt = abs((f1->GetParameter(0)-f2->GetParameter(0)));
    double vlt_err = sqrt(pow(f1->GetParError(0),2)+pow(f2->GetParError(0),2));

    double crg = (3.3)*vlt;
    double crg_err = crg*sqrt(pow(vlt_err/vlt,2)+pow(0.5/3.3,2));

    out.push_back(crg);
    out.push_back(crg_err);

  }

  return out;
  
}


TGraphErrors* process(TGraph* process_graph){

  TGraphErrors* graph = new TGraphErrors();

  int counter = 0;//initialise counter for averaging data

  std::pair< std::pair<double,double> , std::pair<double,double> > aves;

  std::vector<double> tvec;//initialise vector for storing average time
  std::vector<double> tvec_err;//initialise vector for storing time error
  std::vector<double> vvec;//initialise vector for storing average voltage
  std::vector<double> vvec_err;//initialise vector for storing voltage error

  std::vector<double> v_calc;//initialise vector for calculating voltage average and error
  std::vector<double> t_calc;//initialise vector for calculating time average and error
 
  for (int ii=0;ii<process_graph->GetN();ii++){

    v_calc.push_back(process_graph->GetY()[ii]);
    t_calc.push_back(process_graph->GetX()[ii]);

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

  return graph;
  
}

int main(int argc,char* argv[]){

  std::vector<double> ret;

  float res0;
  float res_err0;
  
  float res1;
  float res_err1;

  TGraphErrors* graph_0 = new TGraphErrors();//graph for the preamp peak
  TGraphErrors* graph_1 = new TGraphErrors();//graph for the square wave

  TGraph* raw_graph_0 = new TGraph(argv[1],"%lg,%lg");//raw graph for the preamp peak
  TGraph* raw_graph_1 = new TGraph(argv[2],"%lg,%lg");//raw graph for the square wave

  graph_0 = process(raw_graph_0);
  graph_1 = process(raw_graph_1);

  ret = sq_fit(graph_0,0);

  res0=ret[0];
  res_err0=ret[1];

  ret = sq_fit(graph_1,1);

  res1=ret[0];
  res_err1=ret[1];

  std::ofstream myfile;
  myfile.open ("test.csv",std::ios::app);
  myfile << res1<<","<<res0<<","<<res_err1<<","<<res_err0<<"\n";
  myfile.close();

  //auto c1 = new TCanvas("c1","test",1000,1000);
  
  //graph->SetMarkerStyle(21);
  //graph->SetMarkerColor(1);
  //graph->SetTitle("Oscilloscope Readout;Time (s);Voltage (V)");
  //graph->Draw("AP");
  //c1->SaveAs("test.pdf");
  
  return 0;
}

  
  
