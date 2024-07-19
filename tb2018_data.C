/*
chi2_max[MCP1/2] does not effect the resolution

 */


#include <iostream>
#include <fstream>
#include <iomanip>
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TObject.h"
#include "TPostScript.h"
#include "TRandom.h"

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "TMinuit.h"
#include "TMath.h"

using namespace std;


void tb2018_data() {
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(.08);
  gStyle->SetPadGridX(0); //bool input yes/no, must give an input
  gStyle->SetPadGridY(0);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.11); 
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.02);
  
  gStyle->SetOptFit(101);
  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineStyle(11);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.95);
  gStyle->SetStatFontSize(.040);
  gStyle->SetStatW(0.40);
  gStyle->SetStatH(0.20);
  gStyle->SetOptLogy(0);
  gStyle->SetOptTitle(1);
  gStyle->SetLabelSize(0.070,"XYZ");
  gStyle->SetTitleSize(0.070,"XYZ");
  gStyle->SetTitleOffset(0.70,"XYZ"); 
  gStyle->SetNdivisions(406, "XYZ");
  
  char rootfiles[100];
  char datafile[200];
  char name[100];
  char namex[100];

  const double lowtmcp1=23.0, lowtmcp2=67.0, lowtecal=190.0; //160.0; //190.0;
  const double hghtmcp1=34.0, hghtmcp2=78.0, hghtecal=260.0; //290.0; //260.0;
  const double totalen=4800.0; //Maximum enerry range, this might be given from inu for various energy
  const double pedwidth=2.3; //Average pedestal values of all crystals
  const int nsigma=10; //Look for signals about nsigma
  const int noption=5; //Energy measurement of only seed, 2, 4, 9 and 25 crystals
  
  const int nhitsmx=20;
  int x1n_clusters, x2n_clusters, y1n_clusters, y2n_clusters, n_tracks;
  float x1X[nhitsmx],  x2X[nhitsmx],  y1Y[nhitsmx], y2Y[nhitsmx];
  vector<float> *trkx =0, *trky = 0;   
  
  int              nhits;
  UInt_t           trg;
  Int_t           ped;
  Int_t           phys;
  Int_t           laser;

  Int_t           A1;
  Int_t           A2;
  Int_t           A3;
  Int_t           A4;
  Int_t           A5;
  Int_t           B1;
  Int_t           B2;
  Int_t           B3;
  Int_t           B4;
  Int_t           B5;
  Int_t           C1;
  Int_t           C2;
  Int_t           C3;
  Int_t           C4;
  Int_t           C5;
  Int_t           D1;
  Int_t           D2;
  Int_t           D3;
  Int_t           D4;
  Int_t           D5;
  Int_t           E1;
  Int_t           E2;
  Int_t           E3;
  Int_t           E4;
  Int_t           E5;
  Int_t           MCP1;
  Int_t           MCP2;
  Int_t           CLK;
  Int_t           CFD;
  Int_t           LED;
  ULong64_t       index;
  UInt_t          n_channels;
  UInt_t          n_timetypes;
  Float_t         gain[28];   //[n_channels]
  Float_t         pedestal[28];   //[n_channels]
  Float_t         b_charge[28];   //[n_channels]
  Float_t         b_slope[28];   //[n_channels]
  Float_t         b_rms[28];   //[n_channels]
  Float_t         time[56];   //[n_timetypes]
  Float_t         time_chi2[56];   //[n_timetypes]
  Float_t         time_error[56];   //[n_timetypes]
  Float_t         time_slope[56];   //[n_timetypes]
  Float_t         period[28];   //[n_channels]
  Float_t         maximum[28];   //[n_channels]
  Float_t         time_maximum[28];   //[n_channels]
  Float_t         amp_max[28];   //[n_channels]
  Float_t         amp_max_crystal[25];   //[n_channels]
  Float_t         posX_crystal[25];   //[n_channels]
  Float_t         posY_crystal[25];   //[n_channels]
  Float_t         time_max[28];   //[n_channels]
  Float_t         chi2_max[28];   //[n_channels]
  Float_t         charge_tot[28];   //[n_channels]
  Float_t         charge_sig[28];   //[n_channels]
  Float_t         fit_ampl[28];   //[n_channels]
  Float_t         fit_time[28];   //[n_channels]
  Float_t         fit_terr[28];   //[n_channels]
  Float_t         fit_chi2[28];   //[n_channels]
  Float_t         fit_period[28];   //[n_channels]
  Float_t         fit_ampl_scint[28];   //[n_channels]
  Float_t         fit_time_scint[28];   //[n_channels]
  Float_t         fit_ampl_spike[28];   //[n_channels]
  Float_t         fit_time_spike[28];   //[n_channels]
  Float_t         fit_chi2_scint_plus_spike[28];   //[n_channels]
  Bool_t          fit_converged_scint_plus_spike[28];   //[n_channels]
  Float_t         ampl_calib[28];   //[n_channels]
  Float_t         time_calib[28];   //[n_channels]

  Float_t        E1x1;
  Float_t        amp_max_C2;
  Float_t        E3x3;
  Float_t        E5x5;
  Float_t        R13;
  Float_t        R15;
  Float_t        R35;
  Float_t        dt_C2_MCP1;
  Float_t        dt_C2_MCP2;
  Float_t        dt_MCP1_MCP2; 
  Float_t        dt_C2_C1;
  Float_t        dt_C2_C3;
  Float_t        dt_C2_B2;
  Float_t        dt_C2_D2;
  Float_t        dt_C2_B3;
  Float_t        dt_C2_D3;
  Float_t        dt_C2_B1;
  Float_t        dt_C2_D1;
  Float_t        t_C2_C1;
  Float_t        t_C2_C3;
  Float_t        gain_C2;
  Float_t        gain_MCP1;
  Float_t        gain_MCP2;
  Float_t        dt_C2_C1_phase_correct;
  Float_t        dt_C2_C3_phase_correct;
  Float_t        dt_C2_B1_phase_correct;
  Float_t        dt_C2_B2_phase_correct;
  Float_t        dt_C2_B3_phase_correct;
  Float_t        dt_C2_D1_phase_correct;
  Float_t        dt_C2_D2_phase_correct;
  Float_t        dt_C2_D3_phase_correct;

  Float_t        dt_C3_MCP1;
  Float_t        dt_C3_MCP2;
  Float_t        dt_C3_C4_phase_correct;
  Float_t        dt_C3_C2_phase_correct;
  Float_t        dt_C3_B4_phase_correct;
  Float_t        dt_C3_B2_phase_correct;
  Float_t        dt_C3_B3_phase_correct;
  Float_t        dt_C3_D4_phase_correct;
  Float_t        dt_C3_D2_phase_correct;
  Float_t        dt_C3_D3_phase_correct;

  Float_t        h1x;
  Float_t        h1y;
  Float_t        h2x;
  Float_t        h2y;

  UInt_t          run;
  UInt_t          spill;
  UInt_t          event;
   
  Float_t trackX[30],  trackY[30];

  // int const nfiles = 4;
  // string enrg[nfiles] = {"100","150","200","250"};
  // string target_crystal = "C32x2";

  int const nfiles = 1;
  string enrg[nfiles] = {"pedestal"};
  string target_crystal = "pedestal";

  //string intag = "VFEfix_syncfix_v3_files";
  //string outtag = "VFEfixe_syncfix_v2";
  
  string intag = "digisync_VFE_fix_v2_files";                                                                       
  string outtag = "digisync_VFE_fix_test_files";

  // string intag = "digisync_VFE_fix_v2_files";                                                                       
  // string outtag = "dig_VFE_files";

  for (int ifl=0; ifl<nfiles; ifl++) {

    sprintf(datafile,"./merge_rootfiles/%s_%s_%s.root", enrg[ifl].c_str(), target_crystal.c_str(),intag.c_str());
    //    sprintf(datafile, "/data/gobinda/anal/ecal/tb2018/%i/%i.root", irun, ifl+1);

    TFile* fileIn = new TFile(datafile, "read");
    TTree *T1X = (TTree*)fileIn->Get("h1X");
    TTree *T1Y = (TTree*)fileIn->Get("h1Y");
    TTree *T2X = (TTree*)fileIn->Get("h2X");
    TTree *T2Y = (TTree*)fileIn->Get("h2Y");
    TTree *Ttrk = (TTree*)fileIn->Get("track_tree");

    TTree *Thtrg = (TTree*)fileIn->Get("trg");
    Thtrg->SetBranchAddress("PHYS", &phys);
    Thtrg->SetBranchAddress("PED", &ped);
    Thtrg->SetBranchAddress("LASER", &laser);
    Thtrg->SetBranchAddress("trg", &trg);


    TTree *Th4 = (TTree*)fileIn->Get("h4");
    Th4->SetBranchAddress("run", &run);
    Th4->SetBranchAddress("spill", &spill);
    Th4->SetBranchAddress("event", &event);
 
    T1X->SetBranchAddress("n_clusters", &x1n_clusters);   
    T1X->SetBranchAddress("clusters.X_", x1X);   

    T2X->SetBranchAddress("n_clusters", &x2n_clusters);   
    T2X->SetBranchAddress("clusters.X_", x2X);   

    T1Y->SetBranchAddress("n_clusters", &y1n_clusters);   
    T1Y->SetBranchAddress("clusters.Y_", y1Y);   

    T2Y->SetBranchAddress("n_clusters", &y2n_clusters);   
    T2Y->SetBranchAddress("clusters.Y_", y2Y);   

    Ttrk->SetBranchAddress("n_tracks", &n_tracks);   
    Ttrk->SetBranchAddress("X", &trkx);   
    Ttrk->SetBranchAddress("Y", &trky);   
    
    TTree *Tin = (TTree*)fileIn->Get("digi");
    
    Tin->SetBranchAddress("A1", &A1);
    Tin->SetBranchAddress("A2", &A2);
    Tin->SetBranchAddress("A3", &A3);
    Tin->SetBranchAddress("A4", &A4);
    Tin->SetBranchAddress("A5", &A5);
    Tin->SetBranchAddress("B1", &B1);
    Tin->SetBranchAddress("B2", &B2);
    Tin->SetBranchAddress("B3", &B3);
    Tin->SetBranchAddress("B4", &B4);
    Tin->SetBranchAddress("B5", &B5);
    Tin->SetBranchAddress("C1", &C1);
    Tin->SetBranchAddress("C2", &C2);
    Tin->SetBranchAddress("C3", &C3);
    Tin->SetBranchAddress("C4", &C4);
    Tin->SetBranchAddress("C5", &C5);
    Tin->SetBranchAddress("D1", &D1);
    Tin->SetBranchAddress("D2", &D2);
    Tin->SetBranchAddress("D3", &D3);
    Tin->SetBranchAddress("D4", &D4);
    Tin->SetBranchAddress("D5", &D5);
    Tin->SetBranchAddress("E1", &E1);
    Tin->SetBranchAddress("E2", &E2);
    Tin->SetBranchAddress("E3", &E3);
    Tin->SetBranchAddress("E4", &E4);
    Tin->SetBranchAddress("E5", &E5);
    Tin->SetBranchAddress("MCP1", &MCP1);
    Tin->SetBranchAddress("MCP2", &MCP2);
    Tin->SetBranchAddress("CLK", &CLK);
    Tin->SetBranchAddress("CFD", &CFD);
    Tin->SetBranchAddress("LED", &LED);
    Tin->SetBranchAddress("index", &index);
    Tin->SetBranchAddress("n_channels", &n_channels);
    Tin->SetBranchAddress("n_timetypes", &n_timetypes);
    Tin->SetBranchAddress("gain", gain);
    Tin->SetBranchAddress("pedestal", pedestal);
    Tin->SetBranchAddress("b_charge", b_charge);
    Tin->SetBranchAddress("b_slope", b_slope);
    Tin->SetBranchAddress("b_rms", b_rms);
    Tin->SetBranchAddress("time", time);
    Tin->SetBranchAddress("time_chi2", time_chi2);
    Tin->SetBranchAddress("time_error", time_error);
    Tin->SetBranchAddress("time_slope", time_slope);
    Tin->SetBranchAddress("period", period);
    Tin->SetBranchAddress("maximum", maximum);
    Tin->SetBranchAddress("time_maximum", time_maximum);
    Tin->SetBranchAddress("amp_max", amp_max);
    Tin->SetBranchAddress("time_max", time_max);
    Tin->SetBranchAddress("chi2_max", chi2_max);
    Tin->SetBranchAddress("charge_tot", charge_tot);
    Tin->SetBranchAddress("charge_sig", charge_sig);
    Tin->SetBranchAddress("fit_ampl", fit_ampl);
    Tin->SetBranchAddress("fit_time", fit_time);
    Tin->SetBranchAddress("fit_terr", fit_terr);
    Tin->SetBranchAddress("fit_chi2", fit_chi2);
    Tin->SetBranchAddress("fit_period", fit_period);
    Tin->SetBranchAddress("fit_ampl_scint", fit_ampl_scint);
    Tin->SetBranchAddress("fit_time_scint", fit_time_scint);
    Tin->SetBranchAddress("fit_ampl_spike", fit_ampl_spike);
    Tin->SetBranchAddress("fit_time_spike", fit_time_spike);
    Tin->SetBranchAddress("fit_chi2_scint_plus_spike", fit_chi2_scint_plus_spike);
    Tin->SetBranchAddress("fit_converged_scint_plus_spike", fit_converged_scint_plus_spike);
    Tin->SetBranchAddress("ampl_calib", ampl_calib);
    Tin->SetBranchAddress("time_calib", time_calib);

   
    // cout <<"Give the input file name"<<endl;;
    // cin>> rootfiles;
  
    // ifstream file_db;
    // file_db.open(rootfiles);  

    // int len=strlen(rootfiles)-4;
    // strncpy(namex, rootfiles, len);
    // namex[len]='\0';

    // cout <<"len "<< len <<" "<<rootfiles<<" "<< namex<<" "<<strlen(rootfiles)<<" "<<strlen(namex)<<endl;
    //  sprintf(name, "output_%s.root", namex);
   
    //  cout <<"Give the  number of files"<<endl;;
    //  cin>>  nfiles;

    sprintf(name, "%s_%s_%s.root",enrg[ifl].c_str(),target_crystal.c_str(),outtag.c_str());
     
    TFile* fileOut = new TFile(name, "recreate");  
    TTree* Tout = new TTree("new_tree", "new_tree");

    Tout->Branch("run", &run, "run/i");
    Tout->Branch("spill", &spill, "spill/i");
    Tout->Branch("event", &event, "event/i");

    Tout->Branch("trg", &trg, "trg/i");
    Tout->Branch("PED", &ped, "PED/I");
    Tout->Branch("laser", &laser, "laser/I");
    Tout->Branch("phys", &phys, "phys/I");
    
    Tout->Branch("nhits", &nhits, "nhits/I");  
    Tout->Branch("A1", &A1, "A1/I");
    Tout->Branch("A2", &A2, "A2/I");
    Tout->Branch("A3", &A3, "A3/I");
    Tout->Branch("A4", &A4, "A4/I");
    Tout->Branch("A5", &A5, "A5/I");
    Tout->Branch("B1", &B1, "B1/I");
    Tout->Branch("B2", &B2, "B2/I");
    Tout->Branch("B3", &B3, "B3/I");
    Tout->Branch("B4", &B4, "B4/I");
    Tout->Branch("B5", &B5, "B5/I");
    Tout->Branch("C1", &C1, "C1/I");
    Tout->Branch("C2", &C2, "C2/I");
    Tout->Branch("C3", &C3, "C3/I");
    Tout->Branch("C4", &C4, "C4/I");
    Tout->Branch("C5", &C5, "C5/I");
    Tout->Branch("D1", &D1, "D1/I");
    Tout->Branch("D2", &D2, "D2/I");
    Tout->Branch("D3", &D3, "D3/I");
    Tout->Branch("D4", &D4, "D4/I");
    Tout->Branch("D5", &D5, "D5/I");
    Tout->Branch("E1", &E1, "E1/I");
    Tout->Branch("E2", &E2, "E2/I");
    Tout->Branch("E3", &E3, "E3/I");
    Tout->Branch("E4", &E4, "E4/I");
    Tout->Branch("E5", &E5, "E5/I");
    Tout->Branch("MCP1", &MCP1, "MCP1/I");
    Tout->Branch("MCP2", &MCP2, "MCP2/I");
    Tout->Branch("CLK", &CLK, "CLK/I");
    Tout->Branch("CFD", &CFD, "CFD/I");
    Tout->Branch("LED", &LED, "LED/I");
    Tout->Branch("index", &index, "index/l");
    Tout->Branch("n_channels", &n_channels, "n_channels/i");
    Tout->Branch("n_timetypes", &n_timetypes, "n_timetypes/i");
    Tout->Branch("gain", gain, "gain[28]/F");

    Tout->Branch("pedestal", pedestal, "pedestal[28]/F");
    //Tout->Branch("b_charge", b_charge, "b_charge[28]/F");
    Tout->Branch("b_slope", b_slope, "b_slope[28]/F");
    Tout->Branch("b_rms", b_rms, "b_rms[28]/F");
    Tout->Branch("time", time, "time[56]/F");
    //Tout->Branch("time_chi2", time_chi2, "time_chi2[56]/F");
    Tout->Branch("time_error", time_error, "time_error[56]/F");
    Tout->Branch("time_slope", time_slope, "time_slope[56]/F");
    Tout->Branch("period", period, "period[28]/F");
    Tout->Branch("maximum", maximum, "maximum[28]/F");
    Tout->Branch("time_maximum", time_maximum, "time_maximum[28]/F");
    Tout->Branch("amp_max", amp_max, "amp_max[28]/F");
    Tout->Branch("time_max", time_max, "time_max[28]/F");
    //Tout->Branch("chi2_max", chi2_max, "chi2_max[28]/F");
    //Tout->Branch("charge_tot", charge_tot, "charge_tot[28]/F");
    //Tout->Branch("charge_sig", charge_sig, "charge_sig[28]/F");
    Tout->Branch("fit_ampl", fit_ampl, "fit_ampl[28]/F");
    Tout->Branch("fit_time", fit_time, "fit_time[28]/F");
    Tout->Branch("fit_terr", fit_terr, "fit_terr[28]/F");
    //Tout->Branch("fit_chi2", fit_chi2, "fit_chi2[28]/F");
    Tout->Branch("fit_period", fit_period, "fit_period[28]/F");
    //Tout->Branch("fit_ampl_scint", fit_ampl_scint, "fit_ampl_scint[28]/F");
    //Tout->Branch("fit_time_scint", fit_time_scint, "fit_time_scint[28]/F");
    //Tout->Branch("fit_ampl_spike", fit_ampl_spike, "fit_ampl_spike[28]/F");
    //Tout->Branch("fit_time_spike", fit_time_spike, "fit_time_spike[28]/F");
    //Tout->Branch("fit_chi2_scint_plus_spike", fit_chi2_scint_plus_spike, "fit_chi2_scint_plus_spike[28]/F");
    //Tout->Branch("fit_converged_scint_plus_spike", fit_converged_scint_plus_spike, "fit_converged_scint_plus_spike[28]/F");
    //Tout->Branch("ampl_calib", ampl_calib, "ampl_calib[28]/F");
    //Tout->Branch("time_calib", time_calib, "time_calib[28]/F");

    Tout->Branch("amp_max_C2",&amp_max_C2,"amp_max_C2/F");
    Tout->Branch("E1x1",&E1x1,"E1x1/F");
    Tout->Branch("E3x3",&E3x3,"E3x3/F");
    Tout->Branch("E5x5",&E5x5,"E5x5/F");
    Tout->Branch("R13",&R13,"R13/F");
    Tout->Branch("R15",&R15,"R15/F");
    Tout->Branch("R35",&R35,"R35/F");
    Tout->Branch("dt_C2_MCP1",&dt_C2_MCP1,"dt_C2_MCP1/F");
    Tout->Branch("dt_C2_MCP2",&dt_C2_MCP2,"dt_C2_MCP2/F");
    Tout->Branch("dt_MCP1_MCP2",&dt_MCP1_MCP2,"dt_MCP1_MCP2/F"); 
    Tout->Branch("dt_C2_C1",&dt_C2_C1,"dt_C2_C1/F");
    Tout->Branch("dt_C2_C3",&dt_C2_C3,"dt_C2_C3/F");
    Tout->Branch("dt_C2_B2",&dt_C2_B2,"dt_C2_B2/F");
    Tout->Branch("dt_C2_D2",&dt_C2_D2,"dt_C2_D2/F");
    Tout->Branch("dt_C2_B3",&dt_C2_B3,"dt_C2_B3/F");
    Tout->Branch("dt_C2_D3",&dt_C2_D3,"dt_C2_D3/F");
    Tout->Branch("dt_C2_B1",&dt_C2_B1,"dt_C2_B1/F");
    Tout->Branch("dt_C2_D1",&dt_C2_D1,"dt_C2_D1/F");
    Tout->Branch("t_C2_C1",&t_C2_C1,"t_C2_C1/F");
    Tout->Branch("t_C2_C3",&t_C2_C3,"t_C2_C3/F");
    Tout->Branch("gain_C2",&gain_C2,"gain_C2/F");
    Tout->Branch("gain_MCP1",&gain_MCP1,"gain_MCP1/F");
    Tout->Branch("gain_MCP2",&gain_MCP2,"gain_MCP2/F");
    Tout->Branch("dt_C2_C1_phase_correct",&dt_C2_C1_phase_correct,"dt_C2_C1_phase_correct/F");
    Tout->Branch("dt_C2_C3_phase_correct",&dt_C2_C3_phase_correct,"dt_C2_C3_phase_correct/F");
    Tout->Branch("dt_C2_B1_phase_correct",&dt_C2_B1_phase_correct,"dt_C2_B1_phase_correct/F");
    Tout->Branch("dt_C2_B2_phase_correct",&dt_C2_B2_phase_correct,"dt_C2_B2_phase_correct/F");
    Tout->Branch("dt_C2_B3_phase_correct",&dt_C2_B3_phase_correct,"dt_C2_B3_phase_correct/F");
    Tout->Branch("dt_C2_D1_phase_correct",&dt_C2_D1_phase_correct,"dt_C2_D1_phase_correct/F");
    Tout->Branch("dt_C2_D2_phase_correct",&dt_C2_D2_phase_correct,"dt_C2_D2_phase_correct/F");
    Tout->Branch("dt_C2_D3_phase_correct",&dt_C2_D3_phase_correct,"dt_C2_D3_phase_correct/F");

    Tout->Branch("n_h1X", &x1n_clusters,"n_h1X/I");   
    Tout->Branch("h1X", x1X,"h1X[n_h1X]/F");   

    Tout->Branch("n_h2X", &x2n_clusters,"n_h2X/I");   
    Tout->Branch("h2X", x2X,"h2X[n_h2X]/F");   

    Tout->Branch("n_h1Y", &y1n_clusters,"n_h1Y/I");   
    Tout->Branch("h1Y", y1Y,"h1Y[n_h1Y]/F");   

    Tout->Branch("n_h2Y", &y2n_clusters,"n_h2Y/I");   
    Tout->Branch("h2Y", y2Y,"h2Y[n_h2Y]/F");   

    Tout->Branch("n_tracks", &n_tracks,"n_tracks/I");   
    Tout->Branch("trackX", trackX,"trackX[n_tracks]/F");   
    Tout->Branch("trackY", trackY,"trackY[n_tracks]/F");   
    
    //Tout->Branch("", , "");
    Tout->Branch("amp_max_crystal", amp_max_crystal, "amp_max_crystal[25]/F");
    Tout->Branch("posX_crystal", posX_crystal, "posX_crystal[25]/F");
    Tout->Branch("posY_crystal", posY_crystal, "posY_crystal[25]/F");
    
    Tout->Branch("dt_C3_MCP1",&dt_C3_MCP1,"dt_C3_MCP1/F");
    Tout->Branch("dt_C3_MCP2",&dt_C3_MCP2,"dt_C3_MCP2/F");

    Tout->Branch("dt_C3_C2_phase_correct",&dt_C3_C2_phase_correct,"dt_C3_C2_phase_correct/F");
    Tout->Branch("dt_C3_C4_phase_correct",&dt_C3_C4_phase_correct,"dt_C3_C4_phase_correct/F");
    Tout->Branch("dt_C3_B4_phase_correct",&dt_C3_B4_phase_correct,"dt_C3_B4_phase_correct/F");
    Tout->Branch("dt_C3_B2_phase_correct",&dt_C3_B2_phase_correct,"dt_C3_B2_phase_correct/F");
    Tout->Branch("dt_C3_B3_phase_correct",&dt_C3_B3_phase_correct,"dt_C3_B3_phase_correct/F");
    Tout->Branch("dt_C3_D4_phase_correct",&dt_C3_D4_phase_correct,"dt_C3_D4_phase_correct/F");
    Tout->Branch("dt_C3_D2_phase_correct",&dt_C3_D2_phase_correct,"dt_C3_D2_phase_correct/F");
    Tout->Branch("dt_C3_D3_phase_correct",&dt_C3_D3_phase_correct,"dt_C3_D3_phase_correct/F");
    
   
    // TH2F* h_hodo1xy = new TH2F("hodo1xy", "hodo1xy", 90, -15, 15, 90, -15., 15.); 
    // TH2F* h_hodo2xy = new TH2F("hodo2xy", "hodo2xy", 90, -15, 15, 90, -15., 15.);
    // TH2F* h_track = new TH2F("trackxy", "trackxy", 90, -15, 15, 90, -15., 15.);

    // TH2F* h_crenergyxy = new TH2F("crenergyxy", "crenergyxy", 5, -.5, 4.5, 5, -.5, 4.5);
  
    // const int nreadchn=27; //total number of channels
    // const int ncrystal=25; //total number of crystal
    // const int nbin=120;
    // double xybins[nbin+1]={1.02329, 1.10154, 1.18577, 1.27644, 1.37404, 1.47911, 1.59221, 1.71396, 1.84502, 1.98609, 2.13796, 2.30144, 2.47742, 2.66686, 2.87078, 3.0903, 3.3266, 3.58096, 3.85478, 4.14954, 4.46684, 4.80839, 5.17607, 5.57186, 5.99791, 6.45654, 6.95024, 7.4817, 8.05378, 8.66962, 9.33254, 10.0462, 10.8143, 11.6413, 12.5314, 13.4896, 14.5211, 15.6315, 16.8267, 18.1134, 19.4984, 20.9894, 22.5944, 24.322, 26.1818, 28.1838, 30.3389, 32.6588, 35.156, 37.8443, 40.738, 43.8531, 47.2063, 50.8159, 54.7016, 58.8844, 63.387, 68.2339, 73.4514, 79.0679, 85.1138, 91.622, 98.6279, 106.17, 114.288, 123.027, 132.434, 142.561, 153.462, 165.196, 177.828, 191.426, 206.063, 221.82, 238.781, 257.04, 276.694, 297.852, 320.627, 345.144, 371.535, 399.945, 430.527, 463.447, 498.884, 537.032, 578.096, 622.3, 669.885, 721.107, 776.247, 835.603, 899.498, 968.278, 1042.32, 1122.02, 1207.81, 1300.17, 1399.59, 1506.61, 1621.81, 1745.82, 1879.32, 2023.02, 2177.71, 2344.23, 2523.48, 2716.44, 2924.15, 3147.75, 3388.44, 3647.54, 3926.45, 4226.69, 4549.88, 4897.79, 5272.3, 5675.45, 6109.42, 6576.58, 7079.46};
      
    // TH2F* h_timeamp2x[nreadchn]={0};
    // TH2F* h_ampchi2x[nreadchn]={0};
    // TH2F* h_ampfit_vs_max[nreadchn]={0};



    // const char* variab[nreadchn]= {"A1", "A2", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5", "C1", "C2", "C3", "C4", "C5", "D1", "D2", "D3", "D4", "D5", "E1", "E2", "E3", "E4", "E5", "MCP1", "MCP2"};

    // for (int ij=0; ij<nreadchn; ij++) {
    //   sprintf(name, "time_vs_ampl_%s", variab[ij]);
    //   h_timeamp2x[ij] = new TH2F(name, name, nbin, xybins, 100, (ij>=ncrystal) ? 0.0 :150.0, (ij>=ncrystal) ? 200.0 :450.0); //-0.5, 255.5); //128, -0.5, 4095.5);
    
    //   sprintf(name, "ampl_vs_chi2_%s", variab[ij]);
    //   h_ampchi2x[ij] = new TH2F(name, name,  nbin, xybins, /*128, -0.5, 4095.5,*/ 120, -0.5, 1199.5);

    //   sprintf(name, "ampfit_vs_max_%s", variab[ij]);
    //   h_ampfit_vs_max[ij] = new TH2F(name, name, nbin, xybins, nbin, xybins);
    // }


    // TH2F* h_crtimecorrel[ncrystal]={0};
    // TH2F* h_crenergycorrel[ncrystal]={0};
  
    // for (int ij=0; ij<ncrystal; ij++) {
    //   sprintf(name, "crtimecorrel_%i", ij);
    //   h_crtimecorrel[ij] = new TH2F (name, name, 120, ((ij==0) ? 10.0 : lowtecal), ((ij==0) ? 130.0 : hghtecal), 120, lowtecal, hghtecal);
    //   sprintf(name, "crenergycorrel_%i", ij);
    //   h_crenergycorrel[ij] = new TH2F (name, name, nbin, xybins, nbin, xybins);
    // }
  
    // TH1F* h_energy[noption][nsigma]={{0}};
    // TH1F* h_crtime[noption][nsigma]={{0}};  
  
    // TH2F* h_crtimedif_enr[noption][nsigma]={{0}};  
    // TH2F* h_mcptimedif_enr[4]={0};

    // TH2F* h_e1by9_enr[nsigma]={0}; 
    // TH2F* h_e9by25_enr[nsigma]={0}; 

    // for (int jk=0; jk<noption; jk++) {
    //   for (int kl=0; kl<nsigma; kl++) {
    // 	sprintf(name, "crenergy_op%i_nsig%i", jk, kl);
    // 	h_energy[jk][kl] = new TH1F(name, name, 1200, 0.0, totalen);

    // 	sprintf(name, "crtime_op%i_nsig%i", jk, kl);
    // 	h_crtime[jk][kl] = new TH1F(name, name, 120, lowtecal, hghtecal);

    // 	sprintf(name, "crenergy_vstime_op%i_nsig%i", jk, kl);
    // 	h_crtimedif_enr[jk][kl] = new TH2F(name, name, 120, 0.0, totalen, 120, -70.0, 20.0); //lowtecal, hghtecal);
      
    // 	if (jk==0) {
    // 	  sprintf(name, "h_e1by9_enr_nsig%i", kl);
    // 	  h_e1by9_enr[kl] = new TH2F(name, name, 120, 0.0, totalen, 90, 0.3, 1.0);
        
    // 	  sprintf(name, "h_e9by25_enr_nsig%i", kl);
    // 	  h_e9by25_enr[kl] = new TH2F(name, name, 120, 0.0, totalen, 90, 0.5, 1.0);
    // 	}
    //   }
    // }

    // h_mcptimedif_enr[0] = new TH2F("mcp_timedf_vs_mnamp", "mcp_timedf_vs_mnamp", 120, 0.0, totalen, 120, 43.5, 45.5);
    // h_mcptimedif_enr[1] = new TH2F("mcp_timedf_vs_avamp", "mcp_timedf_vs_avamp", 120, 0.0, totalen, 120, 43.5, 45.5);
    // h_mcptimedif_enr[2] = new TH2F("mcp_timeav_vs_mnamp", "mcp_timeav_vs_mnamp", 120, 0.0, totalen, 120, 10.0, 130.0);
    // h_mcptimedif_enr[3] = new TH2F("mcp_timeav_vs_avamp", "mcp_timeav_vs_avamp", 120, 0.0, totalen, 120, 10.0, 130.0);


  
    // const int nchannel=3; //Study the channels;
    // int timevar[nchannel]={0};
    // int chrgvar[nchannel]={0};
    // const char* varname[nchannel]={"MCP1", "MCP2", "C2"};
    // double tlow[nchannel]={-0.5, -0.5, -0.5};
    // double thgh[nchannel]={1023.5, 1023.5, 1023.5};
    // double clow[nchannel]={-0.5, -0.5, -0.5};
    // double chgh[nchannel]={8191.5, 8191.5, 8191.5};
  
    // TH1F* histtime1x[nchannel]={0};
    // TH1F* histchrg1x[nchannel]={0};
    // TH2F* histtimechrg2x[nchannel]={0};
  
    // TH1F* histtime1dif[nchannel*(nchannel-1)/2]={0};
    // TH1F* histchrg1dif[nchannel*(nchannel-1)/2]={0};    
    // TH2F* histtime2d[nchannel*(nchannel-1)/2]={0};
    // TH2F* histchrg2d[nchannel*(nchannel-1)/2]={0};  

    // int ix=0;
    // for (int ij=0; ij<nchannel; ij++) {
    //   sprintf(name, "time_%s", varname[ij]);
    //   histtime1x[ij] = new TH1F(name, name, 128, tlow[ij], thgh[ij]);
    
    //   sprintf(name, "charge_%s", varname[ij]);
    //   histchrg1x[ij] = new TH1F(name, name, 128, clow[ij], chgh[ij]);

    //   sprintf(name, "charge_vs_chrg_%s", varname[ij]);
    //   histtimechrg2x[ij] = new TH2F(name, name, 128, tlow[ij], thgh[ij], 128, clow[ij], chgh[ij]);

    //   for (int jk=ij+1; jk<nchannel; jk++) {
    //     sprintf(name, "timediff_%s_%s", varname[ij], varname[jk]);
    //     histtime1dif[ix]= new TH1F(name, name, 128, -200.0, 200.0);

    //     sprintf(name, "chrgdiff_%s_%s", varname[ij], varname[jk]);
    //     histchrg1dif[ix]= new TH1F(name, name, 128, -200.0, 200.0);

    //     sprintf(name, "timed2d_%s_%s", varname[ij], varname[jk]);
    //     histtime2d[ix]= new TH2F(name, name,  128, tlow[ij], thgh[ij],  128, tlow[jk], thgh[jk]);

    //     sprintf(name, "chrg2d_%s_%s", varname[ij], varname[jk]);
    //     histchrg2d[ix]= new TH2F(name, name, 128, clow[ij], chgh[ij], 128, clow[jk], chgh[jk]);

    //     ix++;
    //   }
    // }
  
    // while(!(file_db.eof())) {
    
    //   file_db >> datafile;
    //   if (strstr(datafile,"#")) continue;
    
    //   cout <<"datafile = "<<datafile<<endl;
    
    //   if(file_db.eof()) break;
  
    int nentries = Tin->GetEntries();
    cout <<"nentr "<< datafile<<" "<<nentries<<endl;;
    
    for (int iev=0; iev< nentries; iev++) {
      fileIn->cd();
      cout<<"Entries="<<iev<<"/"<<nentries<<" ";
      Th4->GetEntry(iev);
      cout<<run<<" "<<spill<<" "<<event<<endl;
      T1X->GetEntry(iev);
      T1Y->GetEntry(iev);
      // for (int ij=0; ij<x1n_clusters; ij++) {
      //   for (int jk=0; jk<y1n_clusters; jk++) {
      //     h_hodo1xy->Fill(x1X[ij], y1Y[jk]);
      //   }
      // }

      T2X->GetEntry(iev);
      T2Y->GetEntry(iev);
      Thtrg->GetEntry(iev);
      
      // for (int ij=0; ij<x2n_clusters; ij++) {
      //   for (int jk=0; jk<y2n_clusters; jk++) {
      //     h_hodo2xy->Fill(x2X[ij], y2Y[jk]);
      //   }
      // }
      
      Ttrk->GetEntry(iev);
      //      for (int ij=0; ij<n_tracks; ij++) {
      //        cout <<"ij "<< iev<<" "<<ij<<" "<<(&trkx[ij])<<" "<< (&trky[ij])<<endl;
      //        h_track->Fill((&trkx[ij]), (&trky[ij]));
      //      }

      //nhits = 100*min(x1n_clusters, y1n_clusters) + 10*min(x2n_clusters, y2n_clusters) + n_tracks;
      Tin->GetEntry(iev);
      //First the energy in each crystal
      // for (int ij=0; ij<ncrystal; ij++) {
      //   if (fit_ampl[ij]>3*pedwidth) {
      //     h_crenergyxy->Fill(ij%5, ij/5, fit_ampl[ij]);
      //   }
      // }
      
      // fileOut->cd();
      // for (int jk=0; jk<nreadchn; jk++) {
      //   double tmpfitamp = max(1.0, min(xybins[nbin], double(fit_ampl[jk])));
      //   double tmpamp = max(1.0, min(xybins[nbin], double(amp_max[jk])));
        
      //   h_timeamp2x[jk]->Fill(tmpfitamp, time_max[jk]);
      //   h_ampchi2x[jk]->Fill(tmpfitamp, chi2_max[jk]);
      //   h_ampfit_vs_max[jk]->Fill(tmpamp, tmpfitamp);
      // }
      // double mcptime = max(43.501, min(45.499, 1.0*fit_time[MCP2]-fit_time[MCP1]));
      // double mcpampl = max(xybins[0]+0.01, min(xybins[nbin]-1.0, 0.5*(fit_ampl[MCP1]+fit_ampl[MCP2])));
      
      // h_mcptimedif_enr[0]->Fill(min(fit_ampl[MCP1],fit_ampl[MCP2]),  mcptime);
      // h_mcptimedif_enr[1]->Fill(mcpampl,  mcptime);
      
      // mcptime = max(10.001, min(129.999, 0.5*fit_time[MCP2]+fit_time[MCP1]));
      // h_mcptimedif_enr[2]->Fill(min(fit_ampl[MCP1],fit_ampl[MCP2]),  mcptime);
      // h_mcptimedif_enr[3]->Fill(mcpampl,  mcptime);

      
      // //sort out signals according to pulse height;
      // int in_r[ncrystal];
      // double tmp_ampl[ncrystal];
      // for (int ij=0; ij<ncrystal; ij++) {in_r[ij] = ij; tmp_ampl[ij] = fit_ampl[ij];}
      
      // for (int ij=0; ij<ncrystal; ij++) {
      //   for (int jk=ij+1; jk<ncrystal; jk++) {
      //     if (tmp_ampl[ij]<tmp_ampl[jk]) {
      //       int tmpid = in_r[ij];
      //       in_r[ij] = in_r[jk];
      //       in_r[jk] = tmpid;

      //       double tmpp = tmp_ampl[ij];
      //       tmp_ampl[ij] = tmp_ampl[jk];
      //       tmp_ampl[jk] = tmpp; 
      //     }
      //   }
      // }

      // for (int ij=0; ij<ncrystal; ij++) {
      //   h_crenergycorrel[ij]->Fill(((ij==0) ? mcpampl : fit_ampl[in_r[ij-1]]),  fit_ampl[in_r[ij]]);
      //   if (fit_ampl[in_r[ij]]>3*pedwidth) {
      //     h_crtimecorrel[ij]->Fill(((ij==0) ? mcptime : time_max[in_r[ij-1]]),  time_max[in_r[ij]]);
      //   }
      // }
      
      // //      for (int ij=0; ij<ncrystal; ij++) {
      // //        cout <<"ij "<<ij<<" "<<in_r[ij]<<" "<<tmp_ampl[ij]<<" "<<fit_ampl[ij]<<endl;
      // //      }
      
      // //Signals of 1, 2, 4, 9 & 25 crystals
      // //Need proper code for 2, 4 & 9, such that those are surrounding seed crystal          
      // double sumenr[noption][nsigma]={{0}};
      // double sumtime[noption][nsigma]={{0}};
      
      // int ix1=-5, ix2=-5, iy1=-5, iy2=-5;
      // for (int ij=0; ij<ncrystal; ij++) {
      //   int icr=in_r[ij];

      //   //        cout <<ij<<" "<<icr<<" "<<fit_ampl[icr]<<" "<<time_max[icr]<<endl;
        
      //   if (ij==0) {
      //     ix1 = icr%5;  iy1 = icr/5;
      //   } else {
      //     ix2 = icr%5;  iy2 = icr/5;
      //   }
        
      //   for (int jk=0; jk<noption; jk++) {
      //     if (jk==0 && ij>0) continue;
      //     if (jk==1 && ij>1) continue;
      //     if (jk==2 && ij>3) continue;
      //     if (jk==3 && ij>8) continue;
          
      //     if (ij>0 && jk>=1 && jk<=3) {
      //       if (abs(ix1-ix2)>1 || abs(iy1-iy2)>1) continue;
      //     }
          
      //     for (int kl=0; kl<nsigma; kl++) {
      //       if (fit_ampl[icr]>(kl+1)*pedwidth) {
      //         sumenr[jk][kl] += fit_ampl[icr];
      //         sumtime[jk][kl] += fit_ampl[icr]*time_max[icr];
      //       }
      //     }
      //   }
      // }
      // for (int jk=0; jk<noption; jk++) {
      //   for (int kl=0; kl<nsigma; kl++) {
      //     double tmpenr  = max(1.,min(sumenr[jk][kl],totalen-0.5));
      //     h_energy[jk][kl]->Fill(tmpenr);
      //     double tmptime = sumtime[jk][kl]/max(1.,sumenr[jk][kl]);
      //     h_crtimedif_enr[jk][kl]->Fill(tmpenr, max(-69.99, min(19.99,time_max[MCP1]+time_max[MCP2]-tmptime)));
      //     h_crtime[jk][kl]->Fill(max( lowtecal, min( hghtecal, tmptime)));
          
      //     if (jk==0) {
      //       double toten = max(1.,min(sumenr[3][kl],totalen-0.5));
            
      //       //            if (kl==0) { cout <<"toten "<<sumenr[0][kl]<<" "<<sumenr[1][kl]<<" "<<sumenr[2][kl]<<" "<<sumenr[3][kl]<<" "<<sumenr[4][kl]<<" x "<<toten<<" "<<tmpenr/toten<<" "<<toten/max(1.,min(sumenr[4][kl],totalen-0.5))<<endl;}
            
            
      //       h_e1by9_enr[kl]->Fill(toten, max(0.301,min(0.999,tmpenr/toten)));
      //       h_e9by25_enr[kl]->Fill(toten, max(0.501,min(0.999,toten/max(1.,min(sumenr[4][kl],totalen-0.5)))));
      //     }
      //   }
      // }
      
      // // int ivar[nchannel]={MCP1, MCP2, C2};
      // // int ix=0;
      // // for (int ij=0; ij<nchannel; ij++) {
      // //   if (ivar[ij]>=0) {
      // //     histtime1x[ij]->Fill(fit_time[ivar[ij]]);
      // //     histchrg1x[ij]->Fill(amp_max[ivar[ij]]);
      // //     histtimechrg2x[ij]->Fill(fit_time[ivar[ij]], amp_max[ivar[ij]]);
      // //   }
      // //   for (int jk=ij+1; jk<nchannel; jk++) {
      // //     if (ivar[ij]>=0 && ivar[jk]>=0) {
      // //       histtime1dif[ix]->Fill(fit_time[ivar[ij]] - fit_time[ivar[jk]]);
      // //       histtime2d[ix]->Fill(fit_time[ivar[ij]], fit_time[ivar[jk]]);
      // //       histchrg1dif[ix]->Fill(amp_max[ivar[ij]] - amp_max[ivar[jk]]);
      // //       histchrg2d[ix]->Fill(amp_max[ivar[ij]], amp_max[ivar[jk]]);
      // //     }
      // //     ix++;
      // //   }
      // // }

      // //      time 190 to 260
      // // ampmx 2, 20 in 10 steps
      

      for (int ij=0; ij<n_tracks; ij++) {
	trackX[ij] = trkx->at(ij);
	trackY[ij]  = trky->at(ij);
      }

      for(int ij = 0; ij < 25; ij++)
	{
	  amp_max_crystal[ij] = amp_max[ij];
	  posX_crystal[ij] = -100000.0;
	  posY_crystal[ij] = -100000.0;
	}

      posX_crystal[C3] = 0;                          posY_crystal[C3] = 0;

      float crystal_size_x = 11.0;
      float crystal_size_y = crystal_size_x;

      posX_crystal[C1] = posX_crystal[C3];                     posY_crystal[C1] = posY_crystal[C3] -2*crystal_size_y;
      posX_crystal[C2] = posX_crystal[C3];                     posY_crystal[C2] = posY_crystal[C3] -crystal_size_y;
      posX_crystal[C4] = posX_crystal[C3];                     posY_crystal[C4] = posY_crystal[C3] +crystal_size_y;
      posX_crystal[C5] = posX_crystal[C3];                     posY_crystal[C5] = posY_crystal[C3] +2*crystal_size_y;

      posX_crystal[A1] = posX_crystal[C3] -2*crystal_size_x;   posY_crystal[A1] = posY_crystal[C1];
      posX_crystal[A2] = posX_crystal[C3] -2*crystal_size_x;   posY_crystal[A2] = posY_crystal[C2];
      posX_crystal[A3] = posX_crystal[C3] -2*crystal_size_x;   posY_crystal[A3] = posY_crystal[C3];
      posX_crystal[A4] = posX_crystal[C3] -2*crystal_size_x;   posY_crystal[A4] = posY_crystal[C4];
      posX_crystal[A5] = posX_crystal[C3] -2*crystal_size_x;   posY_crystal[A5] = posY_crystal[C5];

      posX_crystal[B1] = posX_crystal[C3] -crystal_size_x;     posY_crystal[B1] = posY_crystal[C1];
      posX_crystal[B2] = posX_crystal[C3] -crystal_size_x;     posY_crystal[B2] = posY_crystal[C2];
      posX_crystal[B3] = posX_crystal[C3] -crystal_size_x;     posY_crystal[B3] = posY_crystal[C3];
      posX_crystal[B4] = posX_crystal[C3] -crystal_size_x;     posY_crystal[B4] = posY_crystal[C4];
      posX_crystal[B5] = posX_crystal[C3] -crystal_size_x;     posY_crystal[B5] = posY_crystal[C5];
      
      posX_crystal[D1] = posX_crystal[C3] +crystal_size_x;     posY_crystal[D1] = posY_crystal[C1];
      posX_crystal[D2] = posX_crystal[C3] +crystal_size_x;     posY_crystal[D2] = posY_crystal[C2];
      posX_crystal[D3] = posX_crystal[C3] +crystal_size_x;     posY_crystal[D3] = posY_crystal[C3];
      posX_crystal[D4] = posX_crystal[C3] +crystal_size_x;     posY_crystal[D4] = posY_crystal[C4];
      posX_crystal[D5] = posX_crystal[C3] +crystal_size_x;     posY_crystal[D5] = posY_crystal[C5];

      posX_crystal[E1] = posX_crystal[C3] +2*crystal_size_x;   posY_crystal[E1] = posY_crystal[C1];
      posX_crystal[E2] = posX_crystal[C3] +2*crystal_size_x;   posY_crystal[E2] = posY_crystal[C2];
      posX_crystal[E3] = posX_crystal[C3] +2*crystal_size_x;   posY_crystal[E3] = posY_crystal[C3];
      posX_crystal[E4] = posX_crystal[C3] +2*crystal_size_x;   posY_crystal[E4] = posY_crystal[C4];
      posX_crystal[E5] = posX_crystal[C3] +2*crystal_size_x;   posY_crystal[E5] = posY_crystal[C5];
      
      amp_max_C2 = amp_max[C2];
      if(target_crystal.find("C2") != string::npos)
	{
	  E1x1 = amp_max[C2];
	  E3x3 = amp_max[C1]+amp_max[C2]+amp_max[C3]+amp_max[B1]+amp_max[B2]+amp_max[B3]+amp_max[D1]+amp_max[D2]+amp_max[D3];
	}
      else if (target_crystal.find("C3") != string::npos)
	{
	  E1x1 = amp_max[C3];
	  E3x3 = amp_max[C2]+amp_max[C3]+amp_max[C4]+amp_max[B2]+amp_max[B3]+amp_max[B4]+amp_max[D2]+amp_max[D3]+amp_max[D4];
	}
      
      E5x5 = amp_max[C1]+amp_max[C2]+amp_max[C3]+amp_max[C4]+amp_max[C5]+amp_max[B1]+amp_max[B2]+amp_max[B3]+amp_max[B4]+amp_max[B5]+amp_max[D1]+amp_max[D2]+amp_max[D3]+amp_max[D4]+amp_max[D5]+amp_max[A1]+amp_max[A2]+amp_max[A3]+amp_max[A4]+amp_max[A5]+amp_max[E1]+amp_max[E2]+amp_max[E3]+amp_max[E4]+amp_max[E5];
      R13 = E1x1/E3x3;
      R15 = E1x1/E5x5;
      R35 = E3x3/E5x5;
      dt_C2_MCP1 = fit_time[C2]-fit_time[MCP1]+fit_time[CLK]-(int((fit_time[C2]-fit_time[MCP1]+fit_time[CLK])/6.238)*6.238);
      dt_C2_MCP2 = fit_time[C2]-fit_time[MCP2]+fit_time[CLK]-(int((fit_time[C2]-fit_time[MCP2]+fit_time[CLK])/6.238)*6.238);
      dt_MCP1_MCP2 = fit_time[MCP1]-fit_time[MCP2]; 
      dt_C2_C1 = fit_time[C2]-fit_time[C1]-(int((fit_time[C2]-fit_time[C1])/6.238)*6.238);
      dt_C2_C3 = fit_time[C2]-fit_time[C3]-(int((fit_time[C2]-fit_time[C3])/6.238)*6.238);
      dt_C2_B2 = fit_time[C2]-fit_time[B2]-(int((fit_time[C2]-fit_time[B2])/6.238)*6.238);
      dt_C2_D2 = fit_time[C2]-fit_time[D2]-(int((fit_time[C2]-fit_time[D2])/6.238)*6.238);
      dt_C2_B3 = fit_time[C2]-fit_time[B3]-(int((fit_time[C2]-fit_time[B3])/6.238)*6.238);
      dt_C2_D3 = fit_time[C2]-fit_time[D3]-(int((fit_time[C2]-fit_time[D3])/6.238)*6.238);
      dt_C2_B1 = fit_time[C2]-fit_time[B1]-(int((fit_time[C2]-fit_time[B1])/6.238)*6.238);
      dt_C2_D1 = fit_time[C2]-fit_time[D1]-(int((fit_time[C2]-fit_time[D1])/6.238)*6.238);
      t_C2_C1 = fit_time[C2]-fit_time[C1];
      t_C2_C3 = fit_time[C2]-fit_time[C3];
      gain_C2 = gain[C2];
      gain_MCP1 = gain[MCP1];
      gain_MCP2 = gain[MCP2];
      
      dt_C2_C1_phase_correct =fit_time[C2]-fit_time[C1]-(round((fit_time[C2]-fit_time[C1])/6.238)*6.238);
      dt_C2_C3_phase_correct = fit_time[C2]-fit_time[C3]-(round((fit_time[C2]-fit_time[C3])/6.238)*6.238);
      dt_C2_B1_phase_correct = fit_time[C2]-fit_time[B1]-(round((fit_time[C2]-fit_time[B1])/6.238)*6.238);
      dt_C2_B2_phase_correct = fit_time[C2]-fit_time[B2]-(round((fit_time[C2]-fit_time[B2])/6.238)*6.238);
      dt_C2_B3_phase_correct = fit_time[C2]-fit_time[B3]-(round((fit_time[C2]-fit_time[B3])/6.238)*6.238);
      dt_C2_D1_phase_correct = fit_time[C2]-fit_time[D1]-(round((fit_time[C2]-fit_time[D1])/6.238)*6.238);
      dt_C2_D2_phase_correct = fit_time[C2]-fit_time[D2]-(round((fit_time[C2]-fit_time[D2])/6.238)*6.238);
      dt_C2_D3_phase_correct = fit_time[C2]-fit_time[D3]-(round((fit_time[C2]-fit_time[D3])/6.238)*6.238);
      
      dt_C3_MCP1 = fit_time[C3]-fit_time[MCP1]+fit_time[CLK]-(int((fit_time[C3]-fit_time[MCP1]+fit_time[CLK])/6.238)*6.238);
      dt_C3_MCP2 = fit_time[C3]-fit_time[MCP2]+fit_time[CLK]-(int((fit_time[C3]-fit_time[MCP2]+fit_time[CLK])/6.238)*6.238);
      dt_C3_C4_phase_correct = fit_time[C3]-fit_time[C4]-(round((fit_time[C3]-fit_time[C4])/6.238)*6.238);
      dt_C3_C2_phase_correct = fit_time[C3]-fit_time[C2]-(round((fit_time[C3]-fit_time[C2])/6.238)*6.238);
      dt_C3_B4_phase_correct = fit_time[C3]-fit_time[B4]-(round((fit_time[C3]-fit_time[B4])/6.238)*6.238);
      dt_C3_B2_phase_correct = fit_time[C3]-fit_time[B2]-(round((fit_time[C3]-fit_time[B2])/6.238)*6.238);
      dt_C3_B3_phase_correct = fit_time[C3]-fit_time[B3]-(round((fit_time[C3]-fit_time[B3])/6.238)*6.238);
      dt_C3_D4_phase_correct = fit_time[C3]-fit_time[D4]-(round((fit_time[C3]-fit_time[D4])/6.238)*6.238);
      dt_C3_D2_phase_correct = fit_time[C3]-fit_time[D2]-(round((fit_time[C3]-fit_time[D2])/6.238)*6.238);
      dt_C3_D3_phase_correct = fit_time[C3]-fit_time[D3]-(round((fit_time[C3]-fit_time[D3])/6.238)*6.238);

      fileOut->cd();
      
      Tout->Fill();

      
    } //end of events loop
    
    fileOut->cd();

    //gStyle->SetOptStat(10);
    // TCanvas* c1 = new TCanvas("c1", "c1", 600, 800);
    // c1->Divide(3,3,1.e-6,1.e-6);
  
    // pss.NewPage();
    // // for (int ij=0; ij<nchannel; ij++) {
    // //   c1->cd(2*ij+1); histtime1x[ij]->Draw();
    // //   c1->cd(2*ij+2); histchrg1x[ij]->Draw();
    // // }
  
    // c1->cd(1); h_hodo1xy->Draw("colz");
    // c1->cd(2); h_hodo2xy->Draw("colz");
    // c1->cd(3); h_track->Draw("colz");
    // c1->cd(4); gPad->SetLogz(1); h_crenergyxy->Draw("colz");
    // c1->cd(5);  gPad->SetLogx(1); h_mcptimedif_enr[0]->Draw("colz");
    // c1->cd(6);  gPad->SetLogx(1); h_mcptimedif_enr[1]->Draw("colz");
    // c1->cd(7);  gPad->SetLogx(1); h_mcptimedif_enr[2]->Draw("colz");
    // c1->cd(8);  gPad->SetLogx(1); h_mcptimedif_enr[3]->Draw("colz");


    // c1->Update();


    // pss.NewPage();
    // gStyle->SetOptStat(0);
    // TCanvas* c1a = new TCanvas("c1a", "c1a", 800, 1000);
    // c1a->Divide(5,5,1.e-6,1.e-6);
  
    // pss.NewPage();
    // for (int ij=0; ij<ncrystal; ij++) {
    //   c1a->cd(ij+1);
    //   h_crtimecorrel[ij]->Draw("colz");
    // }
    // c1a->Update();

    // pss.NewPage();
    // for (int ij=0; ij<ncrystal; ij++) {
    //   c1a->cd(ij+1); gPad->SetLogx(1); gPad->SetLogy(1); 
    //   h_crenergycorrel[ij]->Draw("colz");
    // }
    // c1a->Update();
  

    // pss.NewPage();
    // gStyle->SetOptStat(0);
    // gStyle->SetOptLogx(1);
    // TCanvas* c2 = new TCanvas("c2", "c2", 800, 1000);
    // c2->Divide(5,6,1.e-6,1.e-6);
  
    // for (int ij=0; ij<3; ij++) {
    //   pss.NewPage();
    //   for (int jk=0; jk<nreadchn; jk++) {
    // 	c2->cd(jk+1);
    // 	switch(ij) {
    // 	case 0 : h_timeamp2x[jk]->Draw("colz"); break;
    // 	case 1 : h_ampchi2x[jk]->Draw("colz"); break;
    // 	case 2 : gPad->SetLogy(1); h_ampfit_vs_max[jk]->Draw("colz"); break;
    // 	default : h_timeamp2x[jk]->Draw("colz"); break;
    // 	}
    //   }
    //   c2->Update();
    // }
    // gStyle->SetOptLogy(1);
    // gStyle->SetOptStat(1110);
    // pss.NewPage();
    // TCanvas* c3 = new TCanvas("c3", "c3", 800, 1000);
    // c3->Divide(4,5,1.e-6,1.e-6);
  
    // for (int ij=0; ij<noption; ij++) {
    //   pss.NewPage();
    //   for (int jk=0; jk<nsigma; jk++) {
    // 	c3->cd(2*jk+1); gPad->SetLogx(0); h_energy[ij][jk]->Draw();
    // 	c3->cd(2*jk+2); gPad->SetLogx(0); h_crtime[ij][jk]->Draw();
    //   }
    //   c3->Update();
    // }
  
    // pss.NewPage();
    // gStyle->SetOptStat(0);
    // gStyle->SetOptLogy(0);
    // gStyle->SetOptLogx(0);
    // TCanvas* c3a = new TCanvas("c3a", "c3a", 800, 1000);
    // c3a->Divide(3,4,1.e-6,1.e-6);

    // for (int ij=0; ij<2; ij++) {
    //   pss.NewPage();
    //   for (int jk=0; jk<nsigma; jk++) {
    // 	c3a->cd(jk+1);
    // 	switch(ij) {
    // 	case 0 : h_e1by9_enr[jk]->Draw("colz"); break;
    // 	case 1 : h_e9by25_enr[jk]->Draw("colz"); break;
    // 	default : h_e1by9_enr[jk]->Draw("colz"); break;
    // 	}
    //   }
    //   c3a->Update();
    // }
  
    // pss.NewPage();
    // TCanvas* c4 = new TCanvas("c4", "c4", 800, 1000);
    // c4->Divide(3,4,1.e-6,1.e-6);
  
    // for (int ij=0; ij<noption; ij++) {
    //   pss.NewPage();
    //   for (int jk=0; jk<nsigma; jk++) {
    // 	c4->cd(jk+1); h_crtimedif_enr[ij][jk]->Draw("colz");
    //   }
    //   c4->Update();
    // }
 
    // // for (int ij=0; ij<nchannel; ij++) {
    // //   c1->cd(ij+1); histtimechrg2x[ij]->Draw("colz");
    // // }
    // // c1->Update();
  
    // // pss.NewPage();
    // // for (int ij=0; ij<nchannel*(nchannel-1)/2; ij++) {
    // //   c1->cd(2*ij+1); histtime1dif[ij]->Draw();
    // //   c1->cd(2*ij+2); histtime1dif[ij]->Draw();
    // // }
    // // c1->Update();

    // // pss.NewPage();
    // // for (int ij=0; ij<nchannel*(nchannel-1)/2; ij++) {
    // //   c1->cd(2*ij+1); histtime2d[ij]->Draw("colz");
    // //   c1->cd(2*ij+2); histtime2d[ij]->Draw("colz");
    // // }
    // // c1->Update(); 
    // // pss.NewPage();
  
    // // const int nvarx=30;
    // // const int nchn=28;
    // // const char* variab[nvarx]={"gain", "pedestal", "b_charge", "b_slope", "b_rms", "time", "time_chi2",  "time_error", "time_slope", "period", "maximum", "time_maximum", "amp_max", "time_max", "chi2_max", "charge_tot", "charge_sig", "fit_ampl", "fit_time", "fit_terr", "fit_chi2", "fit_period", "fit_ampl_scint", "fit_time_scint", "fit_ampl_spike", "fit_time_spike", "fit_chi2_scint_plus_spike", "fit_converged_scint_plus_spike", "ampl_calib", "time_calib"};
  
    // // TCanvas* c2 = new TCanvas("c2", "c2", 800, 1000);
    // // c2->Divide(5,6,1.e-6,1.e-6);
  
    // // for (int ij=0; ij<nvarx; ij++) {
    // //   pss.NewPage();
    // //   for (int jk=0; jk<nchn; jk++) {
    // //     c2->cd(jk+1); sprintf(name, "%s[%i]", variab[ij], jk); Tout->Draw(name);
    // //   }
    // //   c2->Update();
    // // }
    // pss.Close();

    fileOut->Write();
    fileOut->Close();

    fileIn->cd();
    fileIn->Close();
    delete Tin;
    delete fileIn;
    
  } //end of file loop

}

/*
  Tb2018->Draw("amp_max", "amp_max>0&&amp_max<100")
  Tb2018->Draw("fit_time[MCP2]", "abs(fit_time[MCP2]-73)<7")
  Tb2018->Draw("fit_time[MCP1]", "abs(fit_time[MCP1]-30)<10")


*/
