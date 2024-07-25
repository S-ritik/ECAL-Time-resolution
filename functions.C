#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit;
#include "My_Style.C"
//#include "/home/ritik/Desktop/latest_boosted_ttbar/all_code_lastest_verison/plotscripts.h"

//#define target_C2
#define target_C3

string files = {"100_digi.root"};

string file_tag = "VFEfixe_syncfix_v2";    /// for intersrystal time res
//string file_tag = "digisync_VFE_fix_files"; 
//string file_tag = "digisync_rearranged_events";  ///  for mcp and crsytal time res
//string file_tag = "digisync_rearranged_events_v2";  ///  for mcp and crsytal time res with all vars
//string file_tag = "digisync_VFE_fix_test_files";

//string energies[] = {"25","50","75","100","125","150","175","200"}; //////////// Low puirty C2 runs
//string amp_max_cut[] = {"amp_max[C2] > 300","amp_max[C2] > 600","amp_max[C2] > 1000","amp_max[C2] > 1200","amp_max[C2] > 2000","amp_max[C2] > 2400","amp_max[C2] > 2500","amp_max[C2] > 3000"};

string energies[] = {"50","100","150","200"}; //////////// High puirty C3 U/L/D/R runs
string amp_max_cut[] = {"amp_max[C3] > 300","amp_max[C3] > 1500","amp_max[C3] > 300","amp_max[C3] > 1500"};
double max_posx[] = {5,-2.5,1,-1};  /// D runs
double max_posy[] = {-5,-2,-5,-11};   /// D runs

// string energies[] = {"100","150","200","250"};   //////////// High puirty C3 runs  
// string amp_max_cut[] = {"amp_max[C3] > 1500","amp_max[C3] > 2300","amp_max[C3] > 3000","amp_max[C3] > 3000"};
// double max_posx[] = {3.5, 3.5, 3.5, 3.5, 1.5};
// double max_posy[] = {-5.5, -4.5, -5.5, -6.5};

int nenergies = sizeof(energies)/sizeof(energies[0]);
string treename = "new_tree";
string calc_treename = "digi";


string evt_sel_cuts_global = {"amp_max[C3] > 400 && amp_max[MCP1] >100 && amp_max[MCP2] > 100"};

string gain_cut[] = {"1","10"};


double gain_switch_calriabation_factor = 1;

map<string, float> crystals_Xpos;

map<string, float> crystals_Ypos;


struct tgraphpair
{
  double x;
  double y;
  double xerror;
  double yerror;
};

bool comparestruct(tgraphpair t1, tgraphpair t2)
{
  return t1.x <t2.x;
}
string replace_all_substrings(string s, string x, string y)
{
  size_t pos = 0;
  while (pos += y.length())
    {
      pos = s.find(x, pos);
      if (pos == std::string::npos) {
	break;
      }
 
      s.replace(pos, x.length(), y);
    }
  return s;
}

void removeDuplicates(vector<double>& vec)
{
  set<double> uniqueElements(vec.begin(), vec.end());
  vec.assign(uniqueElements.begin(), uniqueElements.end());
}

double* CalculateMedian(string fileinname, string treename, string varname, string cut, vector<double>& values)
{
  TFile *filein;
  filein = new TFile(fileinname.c_str(),"read");
  if(!filein->IsOpen())
    {cout<<"Check: file "<<fileinname<<" not present"<<endl;exit(0);}
  TTree *tree = (TTree*)filein->Get(treename.c_str());
  if(tree == NULL)
    {cout<<"Check: file "<<fileinname<<" does not have tree "<<treename<<endl;exit(0);}
  // Get the number of entries in the TTree
  Long64_t nEntries = tree->GetEntries();

  // Create a variable to store the current value of the leaf variable
  double value;

  //TBranch* branch = tree->GetBranch(varname.c_str());
  TTreeFormula* formula = new TTreeFormula("Var",varname.c_str(),tree);
  //branch->SetAddress(&formula);

  //TBranch* branch_sel = tree->GetBranch(cut.c_str());
  TTreeFormula* formula2 = new TTreeFormula("Selection",cut.c_str(),tree);
  //branch_sel->SetAddress(&formula2);

  /* TTreeFormula        *select  = nullptr;
   if (selection && strlen(selection)) {
      select = new TTreeFormula("Selection",selection,fTree);
      if (!select) return -1;
      if (!select->GetNdim()) { delete select; return -1; }
      fFormulaList->Add(select);
      }*/
  // Set the branch address to read the leaf variable
   //tree->SetBranchAddress(varname.c_str(), &value);
   //tree->SetBranchAddress(cut.c_str(), &cutvar);

  // Loop over the entries of the TTree
  for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
    tree->GetEntry(iEntry); // Get the entry
    double result = formula->EvalInstance();
    bool cutvar = formula2->EvalInstance() > 0.5 ? true:false;
    if(cutvar) {
      //cout<<endl<<"var = "<<result<<" condition :"<<cut<<" = "<<cutvar;
      values.push_back(result); // Store the value of the leaf variable
      //cout<<" selected"<<endl;
    }
    //else cout<<" not selected"<<endl;
  }

  // Sort the vector in ascending order
  std::sort(values.begin(), values.end());

  int len = int(values.size());
  // Calculate the median value
  double median;
  if(int(values.size()) > 0){
    if (len % 2 == 0) {
      median = (values[len / 2 - 1] + values[len / 2]) / 2.0; // Average of middle two values
    } else {
      median = values[len / 2]; // Middle value
    }
  }
  cout<<"Median = "<<median<<endl;

  vector<double> values_sets[10];
  for (int i = 0; i < len; ++i)
    {
      int set_nottoinclude = (int)gRandom->Integer(10);
      for(int ij = 0; ij < 10; ij++)
	{
	  if(ij != set_nottoinclude)
	    values_sets[ij].push_back(values[i]);
	}
    }
  double error = 0;
  /* for(int ij = 0; ij < 10; ij++)
    {
      double median_set;

      int len_set = int(values_sets[ij].size());

      if(len_set > 0){
	if (len_set % 2 == 0) {
	  median_set = (values_sets[ij][len_set / 2 - 1] + values_sets[ij][len_set / 2]) / 2.0; // Average of middle two values
	} else {
	  median_set = values_sets[ij][len_set / 2]; // Middle value
	}
      }
      cout<<"Set "<<ij<<" has median "<<median_set<<" with len="<<len_set<<" out of "<<len<<endl;
      error = max(error,fabs(median - median_set));
    }
  */
  cout<<"Median = "<<median<<" with error="<<error<<endl;

  double median_witherror[2] = {median,error};
  cout<<"Median = "<<median_witherror[0]<<" with error="<<median_witherror[1]<<endl;
  
  return median_witherror;

}

double calculate_std_from_array( double arr[], int size) 
{
  double sum = 0;
  for (int i = 0; i < size; ++i) {
    sum += arr[i];
  }
  double mean = sum / size;
  double variance = 0;
  for (int i = 0; i < size; ++i) {
    variance += pow(arr[i] - mean, 2);
  }
  variance /= (size-1);
  return pow(variance,1/2);
}

double calculate_ithquantile_from_array( double arr[], int size, double ith) 
{
  std::sort(arr, arr + size);

  int index = size * ith;

  double min_range = 0;

  for (int i = 0; i < size - index; i++)
    {
      double range= arr[i + index] - arr[i];
      if(i==0)
	min_range = range;

      if(range < min_range)
	min_range = range;
    }
    
  return min_range;
 

}

void remove_zerobincontents(TH1 * hist)
{
  for(int ib =1; ib < hist->GetNbinsX() +1; ib++)
    {
      if(hist->GetBinContent(ib) < 1)
	{
	  hist->SetBinContent(ib,1.1);
	  hist->SetBinError(ib,0.1);
	}
    }
}

void replaceSubstring(std::string& str, const std::string& oldStr, const std::string& newStr) {
    size_t pos = 0;
    while ((pos = str.find(oldStr, pos)) != std::string::npos) {
        str.replace(pos, oldStr.length(), newStr);
        pos += newStr.length();
    }
}
void plotTH2D(TH2D* hist, string filename) {
    // Set up CMS style
    gStyle->SetOptStat(1111); // Turn off statistics box
        //gStyle->SetPalette(kBird); // Set color palette (example: kBird)

    // Create canvas
    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas->SetRightMargin(0.15); // Adjust right margin for axis labels
    //gPad->SetLogz(1);
    // Draw TH3D histogram with "colztext" option
    hist->Draw("colz");

    // Save the plot as a PDF file
    canvas->SaveAs(("plots/" + filename + ".png").c_str());

    // Clean up
    delete canvas;
}

void plot_CB_params(TGraphErrors *tg)
{
  TCanvas *cv;
  TLegend *legv;
  
  cv = tdrCanvas("canv_d",(TH1D*)tg->GetHistogram(),8,0);
  cv->cd();
  gPad->SetGrid();

  tg->GetYaxis()->CenterTitle();
  tg->GetXaxis()->CenterTitle();
  
  //gPad->SetLogy(1);
  tg->SetMarkerSize(1.5);
  tg->SetMarkerStyle(20);
  tg->Draw("P");

  double avg=0;
  int count = 0;
  for(int ip =0; ip < tg->GetN(); ip++)
    {
      cout<<endl<<tg->GetPointY(ip);
      avg += tg->GetPointY(ip);
      if(tg->GetPointY(ip) !=  0.00001) count++;
    }
  avg = avg/count;
  cout<<"Avg = "<<avg<<endl;
  
  TLine *line = new TLine(0,avg,tg->GetXaxis()->GetBinLowEdge(tg->GetN()) + tg->GetXaxis()->GetBinWidth(tg->GetN()),avg);
  line->SetLineColor(1);
  line->Draw("same");

  TPaveText* tbox = new TPaveText(0.2,0.65,0.5,0.9,"BRNDC");
  tbox->AddText(TString::Format("Mean value of paramter is %f",avg));
  tbox->SetFillStyle(0);
  tbox->SetBorderSize(0);
  tbox->Draw("same");

  string plot_name_byvar = tg->GetYaxis()->GetTitle();
  plot_name_byvar = "_" + plot_name_byvar + tg->GetXaxis()->GetTitle();
  plot_name_byvar = replace_all_substrings(plot_name_byvar," ","_"); 
  cv->SaveAs(("plots/for_CBparam_estimate_" + plot_name_byvar + ".png").c_str());
  
}

void editrooplot_up(RooPlot* frame)
{
  //frame->SetLabelSize(0.05);
  frame->GetXaxis()->SetTitleSize(0.05);
  frame->GetXaxis()->CenterTitle();
  frame->GetXaxis()->SetTitleOffset(0.85);
  frame->GetXaxis()->SetLabelSize(0.04);

  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetYaxis()->SetTitleOffset(0.7);
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetLabelSize(0.035);

}

void editrooplot_down(RooPlot* frame)
{
  //frame->SetLabelSize(0.05);
  frame->GetXaxis()->SetTitleSize(0.08);
  frame->GetXaxis()->CenterTitle();
  frame->GetXaxis()->SetTitleOffset(0.95);
  frame->GetXaxis()->SetLabelSize(0.07);

  frame->GetYaxis()->SetTitleSize(0.095);
  frame->GetYaxis()->SetTitleOffset(0.37);
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetLabelSize(0.07);

}

void plot_frames_on_canvas(TCanvas* c, RooPlot* xframe_up, RooPlot* xframe_down)
{
  c->cd(1);
  gPad->SetLogy(0);
  
  gPad->SetTopMargin(0.085);
  gPad->SetRightMargin(0.02);
  gPad->SetBottomMargin(0.04);
  gPad->SetLeftMargin(0.08);
  gPad->SetPad(0,0.35,1,1);

  editrooplot_up(xframe_up);
  xframe_up->SetMaximum( ((gPad->GetLogy()) ? 1700 : 1.9)*xframe_up->GetMaximum());  
  xframe_up->SetMinimum(max(1.0,xframe_up->GetMinimum()));
  xframe_up->Draw();
   
  c->cd(2);

  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.04);
  gPad->SetLeftMargin(0.08);
  gPad->SetPad(0,0,1,0.35);
  gPad->SetBottomMargin(0.2);

  
  editrooplot_down(xframe_down);
  xframe_down->SetMaximum(5);
  xframe_down->SetMinimum(-5);
  xframe_down->Draw();
}

void plot_hist(TH1D *hist_obs, string title,string plotname)
{
  
  TCanvas *cv;
  
  cv = tdrCanvas("canv_d",hist_obs,8,0);
  cv->cd();

  gPad->SetGrid();
  //gPad->SetLogy(1);
  gStyle->SetOptStat(11111);
  
  double max_bincontent = hist_obs->GetMaximum();
  double min_bincontent = hist_obs->GetMinimum();
  hist_obs->SetMaximum(((gPad->GetLogy()) ? 50 : 1.05)*max_bincontent);
  hist_obs->SetMinimum(min(1.15*min_bincontent,(gPad->GetLogy() ? 0.5 : 0.85)*min_bincontent));

  hist_obs->SetFillStyle(0);
  hist_obs->SetFillColor(0);
  hist_obs->SetLineColor(2);

  hist_obs->SetLineWidth(2);
  hist_obs->SetMarkerStyle(20);
  hist_obs->SetMarkerSize(1);
  hist_obs->SetMarkerColor(2);
  hist_obs->GetYaxis()->SetTitleOffset(1.5);
  
  hist_obs->GetYaxis()->CenterTitle();
  hist_obs->GetXaxis()->CenterTitle();
  hist_obs->GetXaxis()->SetTitle(title.c_str());  

  hist_obs->Draw("hist");

  //CMS_lumi( cv, 8, 0 );

  cv->SaveAs(("plots/" + plotname).c_str());
}

void plot_1dhists(int nhist, TH1D *hist_obs[nhist], string dataname[nhist],string plotname)
{
  TCanvas *cv;
  TLegend *legv;
  
  cv = tdrCanvas("canv_d",hist_obs[0],8,0);
  cv->cd();
  legv = tdrLeg(0.35,0.65,0.875,0.915);

  legv->SetTextFont(42);
  legv->SetTextSize(0.045);
  legv->SetBorderSize(0);
  gStyle->SetPaintTextFormat( "1.2f" );
  gPad->SetGrid();

  //gStyle->SetOptStat(111111);
      
  double max_bincontent = hist_obs[0]->GetMaximum();
  double min_bincontent = hist_obs[0]->GetMinimum();
  
  for(int ih = 0; ih < nhist; ih++)
    {
      //hist_obs[ih]->Scale(1.0/max(0.00001,hist_obs[ih]->Integral()));
      hist_obs[ih]->SetFillStyle(0);
      hist_obs[ih]->SetFillColor(0);
      hist_obs[ih]->SetLineColor(ih+2);
      //     hist_obs[ih]->SetLineStyle(ih+2);
      hist_obs[ih]->SetLineWidth(2);
      hist_obs[ih]->SetMarkerStyle(20);
      hist_obs[ih]->SetMarkerSize(1);
      hist_obs[ih]->SetMarkerColor(ih+2);
      hist_obs[ih]->GetYaxis()->SetTitleOffset(1.5);

      hist_obs[ih]->GetYaxis()->CenterTitle();
      hist_obs[ih]->GetXaxis()->CenterTitle();

      legv->AddEntry(hist_obs[ih],dataname[ih].c_str(),"ep");

      max_bincontent = max(max_bincontent,double(hist_obs[ih]->GetMaximum()));
      min_bincontent = max(min_bincontent,double(hist_obs[ih]->GetMinimum()));
    }

  hist_obs[0]->SetMaximum(1.05*max_bincontent);
  hist_obs[0]->SetMinimum(0.85*min_bincontent);
  
  //gPad->SetLogy(1);

  if(nhist == 1) hist_obs[0]->Draw("hist");
  if(nhist > 1) hist_obs[0]->Draw("Pe1");

  for(int ih = 1; ih < nhist; ih++)
    hist_obs[ih]->Draw("Pe1same");

  if(nhist> 1) legv->Draw("same");
  //CMS_lumi( cv, 8, 0 );

  cv->SaveAs(("plots/" + plotname).c_str());
}

double time_aeff_fitfunc2(double *x, double *par)
{
  return (sqrt(  pow(par[0]/(2*x[0]),2) +  pow(par[1],2)/(2.0*x[0]) + pow(par[2],2)  ));  //pow(par[0]/sqrt(energy),2)
}
TF1* time_aeff_fit2()
{
  TF1 *fit1 = new TF1("fit1_time_aeff", time_aeff_fitfunc2 ,0,1000,3);
  fit1->SetParameter(0,6000);
  fit1->SetParameter(1,0);
  fit1->SetParNames("N","S","C");
  //fit1->SetParNames("N","C");
  //fit1->FixParameter(0,12870);
  fit1->SetParLimits(0,0,1e10);
  fit1->SetParLimits(1,0,1e10);
  fit1->SetParLimits(2,0,1e10);
  return fit1;
}

double time_aeff_fitfunc(double *x, double *par)
{
  return (sqrt(  /*pow(par[0]/x[0],2)/(2.0)*/ pow(par[0],2)/(2.0*x[0]) + pow(par[1],2)  )); 
}
TF1* time_aeff_fit()
{
  TF1 *fit1 = new TF1("fit1_time_aeff", time_aeff_fitfunc ,0,1000,2);
  fit1->SetParameter(0,6000);
  fit1->SetParameter(1,0);
  fit1->SetParNames("S","C");
  //fit1->SetParNames("N","C");
  //fit1->FixParameter(0,12870);
  fit1->SetParLimits(0,0,1e10);
  fit1->SetParLimits(1,0,1e10);
  return fit1;
}

double linear_fitfunc_tgraph(double *x, double *par)
{
  return (x[0]*par[0] + par[1] );  //pow(par[0]/sqrt(energy),2)
}
TF1* linear_fit_tgraph()
{
  TF1 *fit1 = new TF1("fit1_time_aeff", linear_fitfunc_tgraph ,0,1000,2);
  fit1->SetParameter(0,6000);
  fit1->SetParameter(1,45);
  fit1->SetParNames("Slope","constant");
  return fit1;
}

void plot_tgraph(TGraphErrors *hist_obs,string plotname, bool plotfit = true, string fitfuncname = "")
{
  
  TCanvas *cv;
  TLegend *legv;
  
  cv = tdrCanvas("canv_d",(TH1D*)hist_obs->GetHistogram(),9,0);
  cv->cd();
  legv = tdrLeg(0.2,0.65,0.875,0.915);

  legv->SetTextFont(42);
  legv->SetTextSize(0.045);
  legv->SetBorderSize(0);
  gStyle->SetPaintTextFormat( "1.2f" );
  gPad->SetGrid();


  //hist_obs->SetMaximum(1.5*hist_obs->GetMaximum());
  //hist_obs->SetMinimum(0.97*hist_obs->GetMinimum());
  hist_obs->GetYaxis()->CenterTitle();
  hist_obs->GetYaxis()->SetTitleOffset(0.8);
  //hist_obs->GetYaxis()->SetTitleSize(0.4);
  hist_obs->GetXaxis()->CenterTitle();
  hist_obs->SetMarkerStyle(20);
  hist_obs->SetMarkerSize(1);
    
  //gPad->SetLogy(1);
  gStyle->SetOptFit(0);
  
  hist_obs->Draw("AP");

  int numParams =0;
  const double *params,*paramErrors;
  TF1 *fit1;
  if(plotfit)
    {
      //TF1 *fit1 = new TF1("fit1", fit_function,0,220,2);
      //fit1->SetParameter(0,1);
      //fit1->SetParameter(1,2);
      if(fitfuncname == "linear")
	{
	  fit1 = linear_fit_tgraph();
	  hist_obs->Fit(fit1->GetName());
	  numParams = fit1->GetNpar();
	  params = fit1->GetParameters();
	  paramErrors = fit1->GetParErrors();
	}
      else
	{
	  fit1 = time_aeff_fit2();
	  hist_obs->Fit(fit1->GetName());
	  numParams = fit1->GetNpar();
	  params = fit1->GetParameters();
	  paramErrors = fit1->GetParErrors();
	}
    }
  TLatex *latex  = new TLatex();
  latex->SetNDC();
  latex->SetTextAngle(0);
  latex->SetTextColor(kBlack);    

  latex->SetTextFont(42);
  latex->SetTextAlign(31); 
  //latex.SetTextSize(lumiTextSize*t);
  
  TString iText = " test";

  if(fit1 != NULL)
    {

      if(numParams == 2)
	iText = TString::Format("#splitline{%s = (%1.1f #pm %1.1f) ps}{%s = (%1.1f #pm %1.1f) ps}",fit1->GetParName(0),params[0],paramErrors[0],fit1->GetParName(1),params[1],paramErrors[1]);

      if(numParams == 3)
	iText = TString::Format("#splitline{%s = (%1.1f #pm %1.1f) ps}{#splitline{%s = (%1.1f #pm %1.1f) ps}{%s = (%1.1f #pm %1.1f) ps}}",fit1->GetParName(0),params[0],paramErrors[0],fit1->GetParName(1),params[1],paramErrors[1],fit1->GetParName(2),params[2],paramErrors[2]);

      if(fitfuncname == "linear")
	iText = TString::Format("#splitline{%s = (%1.1f #pm %1.1f) }{%s = (%1.1f #pm %1.1f)}",fit1->GetParName(0),params[0],paramErrors[0],fit1->GetParName(1),params[1],paramErrors[1]);
    }
  if(fit1 != NULL)
    latex->DrawLatex(0.8,0.7,iText);

  if(fit1 != NULL)
    latex->Draw("same");
  CMS_lumi( cv, 9, 0 );
  cv->SaveAs(("plots/" + plotname).c_str());
  replaceSubstring(plotname,"png","root");
  cv->SaveAs(("plots/" + plotname).c_str());
  replaceSubstring(plotname,"root","C");
  cv->SaveAs(("plots/" + plotname).c_str());
  TFile f("plots/allfile.root","RECREATE");
  f.cd();
  hist_obs->Write();
  f.Close();
}


double energy_res_fitfunc(double *x, double *par)
{
  double energy = 25*x[0] + 12.5;
  return (sqrt(   pow(par[0]/energy,2) +  pow(par[1],2)/energy + pow(par[2],2)  ));  //pow(par[0]/sqrt(energy),2)
}
TF1* energy_res_fit()
{
  TF1 *fit1 = new TF1("fit_energy_res", energy_res_fitfunc ,0,1000,3);
  fit1->SetParameter(0,1);
  fit1->SetParameter(1,2);
  fit1->SetParNames("S","N","C");
  return fit1;
}
double time_res_fitfunc(double *x, double *par)
{
  double energy = 25*x[0] + 12.5;
  return (sqrt(   pow(par[0]/energy,2) + pow(par[1],2)  ));  //pow(par[0]/sqrt(energy),2)
}
TF1* time_res_fit()
{
  TF1 *fit1 = new TF1("fit_time_res", time_res_fitfunc ,0,1000,2);
  fit1->SetParameter(0,1);
  fit1->SetParameter(1,2);
  fit1->SetParNames("N","C");
  return fit1;
}

double linear_fitfunc(double *x, double *par)
{
  double energy = 25*x[0] + 12.5;
  return par[0]*energy + par[1];
}
TF1* linear_fit()
{
  TF1 *fit1 = new TF1("fit_linear", linear_fitfunc,0,1000,2);
  fit1->SetParameter(0,1);
  fit1->SetParameter(1,2);
  fit1->SetParNames("Slope","Intercept");
  return fit1;
}

double time_res_vsaff_fitfunc(double *x, double *par)
{
  return (sqrt(   pow(par[0]/x[0],2) + pow(par[1],2)  ));  //pow(par[0]/sqrt(energy),2)
}
TF1* time_res_vsaff_fit()
{
  TF1 *fit1 = new TF1("fit_time_res", time_res_fitfunc ,0,1000,2);
  fit1->SetParameter(0,1);
  fit1->SetParameter(1,2);
  fit1->SetParNames("N","C");
  return fit1;
}

void plot_1dhists_withfit(int nhist, TH1D *hist_obs[nhist], string dataname[nhist],string plotname)
{
  
  TCanvas *cv;
  TLegend *legv;
  
  cv = tdrCanvas("canv_d",hist_obs[0],8,0);
  cv->cd();
  legv = tdrLeg(0.35,0.65,0.875,0.915);

  legv->SetTextFont(42);
  legv->SetTextSize(0.045);
  legv->SetBorderSize(0);
  gStyle->SetPaintTextFormat( "1.2f" );
  gPad->SetGrid();

  //int nlast = sizeof(range)/sizeof(range[0]);

  cout<<"fit range = "<<hist_obs[0]->GetBinCenter(0)<<"  "<<hist_obs[0]->GetBinCenter(1);
  //TF1 *fit1 = new TF1("fit1", linear_fitfunc,0,1000,2);
  //fit1->SetParameter(0,1);
  //fit1->SetParameter(1,2);
  //fit1->SetParameter(2,3);

  TF1 *fit1 = time_res_vsaff_fit();
  //TF1 *fit1 = energy_res_fit();

  double max_bincontent = hist_obs[0]->GetMaximum();
  double min_bincontent = hist_obs[0]->GetMinimum();

  for(int ih = 0; ih < nhist; ih++)
    {
      //hist_obs[ih]->Scale(1.0/max(0.00001,hist_obs[ih]->Integral()));
      hist_obs[ih]->SetFillStyle(0);
      hist_obs[ih]->SetFillColor(0);
      hist_obs[ih]->SetLineColor(ih+2);
      //     hist_obs[ih]->SetLineStyle(ih+2);
      hist_obs[ih]->SetLineWidth(2);
      hist_obs[ih]->SetMarkerStyle(20);
      hist_obs[ih]->SetMarkerColor(ih+2);
      hist_obs[ih]->SetMarkerSize(1);
      hist_obs[ih]->GetYaxis()->SetTitleOffset(1.8);
      hist_obs[ih]->GetYaxis()->CenterTitle();
      hist_obs[ih]->GetXaxis()->CenterTitle();

      legv->AddEntry(hist_obs[ih],dataname[ih].c_str(),"ep");

      max_bincontent = max(max_bincontent,double(hist_obs[ih]->GetMaximum()));
      min_bincontent = max(min_bincontent,double(hist_obs[ih]->GetMinimum()));
    }

  hist_obs[0]->SetMaximum(1.26*max_bincontent);
  hist_obs[0]->SetMinimum(0.55*min_bincontent);

  TH1D * hist_fit = (TH1D*)hist_obs[0]->Clone();
  if(nhist> 1)
    {
      hist_fit->SetBinContent(7,hist_obs[1]->GetBinContent(7));
      hist_fit->SetBinError(7,hist_obs[1]->GetBinError(7));
      hist_fit->SetBinContent(8,hist_obs[1]->GetBinContent(8));
      hist_fit->SetBinError(8,hist_obs[1]->GetBinError(8));
      hist_fit->SetLineColor(1);
      hist_fit->SetLineWidth(2);
    }
  

  //gPad->SetLogy(1);

  hist_fit->Fit(fit1->GetName());

  hist_obs[0]->Draw("Pe1same");
  
  for(int ih = 1; ih < nhist; ih++)
    hist_obs[ih]->Draw("Pe1same");

  if(nhist> 1) legv->Draw("same");
  //CMS_lumi( cv, 8, 0 );

  cv->SaveAs(("plots/" + plotname).c_str());
}

void savehist2d(TH2D *hist, string plotname, bool plotwitherror = true, bool plotlinesx = false,float x1 = 0 , float x2 = 0, bool plotlinesy = false,float y1 = 0 , float y2 = 0)
{
  TCanvas *c1 = tdrCanvas("canv_d", hist->ProjectionX(),8,0);

  c1->cd();
  gStyle->SetOptTitle(0);

  gStyle->SetPaintTextFormat( "1.f" );
  gPad->SetLogx(0);
  gPad->SetLogy(0);

  //  c1->Update();
  hist->GetXaxis()->SetNdivisions(215);
  hist->GetYaxis()->SetNdivisions(215);
  hist->GetXaxis()->SetTickSize(0.03);
  hist->GetYaxis()->SetTickSize(0.03);
  hist->GetXaxis()->SetMoreLogLabels();
  hist->GetYaxis()->SetMoreLogLabels();

  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.06);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetLabelOffset(0.001);

  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.13);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->SetMarkerSize(1.4);
  //hist->GetYaxis()->SetLabelOffset(1.06);
  if(plotwitherror) hist->Draw("colztext");
  else hist->Draw("colz");
  //CMS_lumi( c1, 8, 0 );

  if(plotlinesx)
    {
      TLine *linex1 = new TLine(x1,hist->GetYaxis()->GetXmin(),x1,hist->GetYaxis()->GetXmax());
      linex1->SetLineColor(kRed);
      linex1->Draw("same");

      TLine *linex2 = new TLine(x2,hist->GetYaxis()->GetXmin(),x2,hist->GetYaxis()->GetXmax());
      linex2->SetLineColor(kRed);
      linex2->Draw("same");
    }

    if(plotlinesy)
    {
      TLine *liney1 = new TLine(hist->GetXaxis()->GetXmin(),y1,hist->GetXaxis()->GetXmax(),y1);
      liney1->SetLineColor(kRed);
      liney1->Draw("same");

      TLine *liney2 = new TLine(hist->GetXaxis()->GetXmin(),y2,hist->GetXaxis()->GetXmax(),y2);
      liney2->SetLineColor(kRed);
      liney2->Draw("same");
    }

  c1->SaveAs(("./plots/" + plotname + ".png").c_str());
}

double* get_resolution_value_gauss_const(TH1D *hist, string engy, float xlow, float xhigh)
{
  cout<<"Start calculation of hist:"<<hist->GetName()<<" with xlow = "<<xlow<<" and xhigh = "<<xhigh<<endl;
  
  //float xlow = hist->GetBinLowEdge(1);
  int nbins = hist->GetNbinsX();
  int nfreeparams = 3;

  string title = hist->GetTitle();
  
  RooRealVar x("x",title.c_str(),xlow,xhigh);
  RooRealVar mean("mean","mean of gaussian",(xlow+xhigh)/2,xlow,xhigh);
  RooRealVar sigma("sigma","width of gaussian",(xhigh-xlow)/4,(xhigh-xlow)/40,xhigh-xlow);
  
  RooGaussian sig("gauss","gaussian PDF",x,mean,sigma);  
  
  RooPlot* xframe = x.frame(Title("Gaussian p.d.f."));    
  
  mean.setVal((xlow+xhigh)/2);
  sigma.setVal((xhigh-xlow)/2);
  sig.plotOn(xframe);
  
  RooRealVar sig_yield("sig_yield","signal yield",0.5,0.,1/*2*hist->Integral()*/);
  
  RooPolynomial bkg("pol0","pol0",x,RooArgList());
  
  RooAddPdf model1("model","pol0+gauss",RooArgList(sig,bkg),RooArgList(sig_yield));
  
  model1.plotOn(xframe,LineColor(kRed));
  
  RooDataHist *data = new RooDataHist("data","dataset with (x)",RooArgList(x),hist);
  //data->plotOn(xframe);
  
  RooPlot* xframe2 = x.frame(Title(("Fit of " + title).c_str()));
  data->plotOn(xframe2, DataError(RooAbsData::SumW2));

  /// Change model here

  model1.fitTo(*data);
  model1.plotOn(xframe2,Components(bkg), LineColor(kRed));
  model1.plotOn(xframe2,Components(sig), LineColor(kGreen));
  model1.plotOn(xframe2);
    //model1.paramOn(xframe2,data);
  model1.paramOn(xframe2);

  Double_t skew = data->skewness(x); // 3rd standardized moment
  Double_t kurt = data->kurtosis(x); // 4th standardized moment

  TPaveText* tbox = new TPaveText(0.1,0.65,0.5,0.9,"BRNDC");
  tbox->AddText(TString::Format("#splitline{#splitline{#splitline{At energy: %s,}{Fit Function is %s}}{#frac{#chi^2}{ndof} = %.2f & Prob. = %.2f}}{Skewness = %.2f Kurtosis = %.2f}",engy.c_str(),"Gaussian + Constant Background",xframe2->chiSquare(nfreeparams),TMath::Prob(xframe2->chiSquare(nfreeparams)*(nbins-nfreeparams),(nbins-nfreeparams)),skew,kurt));
  tbox->SetFillStyle(0);
  tbox->SetBorderSize(0);
  xframe2->addObject(tbox);

  //mean.Print();
  //sigma.Print();
  //sig_yield.Print();
  
  //cout << "chi^2 = " << xframe2->chiSquare() << endl;
  
  RooHist *hresid = xframe2->residHist();
  
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist *hpull = xframe2->pullHist();
 
  RooPlot* xframe3 = x.frame(Title(" "));
  xframe3->addPlotable(hpull, "P");
  xframe3->GetYaxis()->SetTitle("Pull");
  
  TCanvas* c = new TCanvas("fit_canv","fit_canv",800,800);   //tdrDiCanvas("c",hist,hist,8,0);
  c->Divide(1,2);
  plot_frames_on_canvas(c,xframe2,xframe3);


  //delete c;
  
  cout<<"End calculation"<<endl;
  engy = replace_all_substrings(engy," ","_");  c->SaveAs(("plots/energy_" + engy + "_" + title + ".png").c_str());
  
  double return_values[6] ={double(sigma.getVal()),double(sigma.getError()),mean.getVal(),double(mean.getError()),sig_yield.getVal(),double(sig_yield.getError())};
  return return_values;

}

double* get_resolution_value_modgauss(TH1D *hist, string engy, float xlow, float xhigh)
{
  cout<<"Start calculation of hist:"<<hist->GetName()<<" with xlow = "<<xlow<<" and xhigh = "<<xhigh<<endl;
  
  //float xlow = hist->GetBinLowEdge(1);
  int nbins = hist->GetNbinsX();
  int nfreeparams = 3;

  string title = hist->GetTitle();
  
  RooRealVar x("x",title.c_str(),xlow,xhigh);
  RooRealVar mean("mean","mean of gaussian",(xlow+xhigh)/2,xlow,xhigh);
  RooRealVar sigma("sigma","width of gaussian",(xhigh-xlow)/4,(xhigh-xlow)/40,xhigh-xlow);
  RooRealVar N("N","Power of gaussian",2,0.01,10);
  
  RooGenericPdf model1("gauss","exp(-1.0*pow((abs(x-mean)/sigma),N)/N)",RooArgSet(x,mean,sigma,N));  
  
  RooPlot* xframe = x.frame(Title("Gaussian p.d.f."));
  
  mean.setVal((xlow+xhigh)/2);
  sigma.setVal((xhigh-xlow)/4);
  N.setVal(2.5);
  model1.plotOn(xframe);
  
  RooRealVar sig_yield("sig_yield","signal yield",0.5,0.,1/*2*hist->Integral()*/);
  
  RooPolynomial bkg("pol0","pol0",x,RooArgList());
  
  //RooAddPdf model1("model","pol0+gauss",RooArgList(sig,bkg),RooArgList(sig_yield));
  
  model1.plotOn(xframe,LineColor(kRed));
  
  RooDataHist *data = new RooDataHist("data","dataset with (x)",RooArgList(x),hist);
  //data->plotOn(xframe);
  
  RooPlot* xframe2 = x.frame(Title(("Fit of " + title).c_str()));
  data->plotOn(xframe2, DataError(RooAbsData::SumW2));

  /// Change model here

  
  model1.fitTo(*data);
  model1.plotOn(xframe2);
  model1.paramOn(xframe2,Layout(0.51,0.96,0.918));

  Double_t skew = data->skewness(x); // 3rd standardized moment
  Double_t kurt = data->kurtosis(x); // 4th standardized moment

  TPaveText* tbox = new TPaveText(0.1,0.65,0.5,0.9,"BRNDC");
  tbox->AddText(TString::Format("#splitline{#splitline{#splitline{At energy: %s,}{Fit Function is %s}}{#frac{#chi^2}{ndof} = %.2f & Prob. = %.2f}}{Skewness = %.2f Kurtosis = %.2f}",engy.c_str(),"Modified Gaussian",xframe2->chiSquare(nfreeparams),TMath::Prob(xframe2->chiSquare(nfreeparams)*(nbins-nfreeparams),(nbins-nfreeparams)),skew,kurt));
  tbox->SetFillStyle(0);
  tbox->SetBorderSize(0);
  xframe2->addObject(tbox) ;

  //mean.Print();
  //sigma.Print();
  //sig_yield.Print();
  
  //cout << "chi^2 = " << xframe2->chiSquare() << endl;
  
  RooHist *hresid = xframe2->residHist();
  
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist *hpull = xframe2->pullHist();
 
  RooPlot* xframe3 = x.frame(Title(" "));
  xframe3->addPlotable(hpull, "P");
  xframe3->GetYaxis()->SetTitle("Pull");
  
  TCanvas* c = new TCanvas("fit_canv","fit_canv",800,800);   //tdrDiCanvas("c",hist,hist,8,0);
  c->Divide(1,2);
  plot_frames_on_canvas(c,xframe2,xframe3);

  //delete c;
  
  cout<<"End calculation"<<endl;
  engy = replace_all_substrings(engy," ","_");  c->SaveAs(("plots/energy_" + engy + "_" + title + ".png").c_str());  
  double reutn_values[6] ={double(sigma.getVal()),double(sigma.getError()),mean.getVal(),double(mean.getError()),double(N.getVal()),double(N.getError())};
  return reutn_values;
}

double* get_resolution_value_gauss(TH1D *hist, string engy, float xlow, float xhigh)
{
  cout<<"Start calculation of hist:"<<hist->GetName()<<" with xlow = "<<xlow<<" and xhigh = "<<xhigh<<endl;
  
  //float xlow = hist->GetBinLowEdge(1);
  int nbins = hist->GetNbinsX();
  int nfreeparams = 2;

  string title = hist->GetTitle();
  
  RooRealVar x("x",title.c_str(),xlow,xhigh);
  RooRealVar mean("mean","mean of gaussian",(xlow+xhigh)/2,xlow,xhigh);
  RooRealVar sigma("sigma","width of gaussian",(xhigh-xlow)/4,(xhigh-xlow)/40,xhigh-xlow);
  
  RooGaussian model1("gauss","gaussian PDF",x,mean,sigma);  
  
  RooPlot* xframe = x.frame(Title("Gaussian p.d.f."));    
  
  mean.setVal((xlow+xhigh)/2);
  sigma.setVal((xhigh-xlow)/4);
  model1.plotOn(xframe);
  
  RooRealVar sig_yield("sig_yield","signal yield",0.5,0.,1/*2*hist->Integral()*/);
  
  RooPolynomial bkg("pol0","pol0",x,RooArgList());
  
  //RooAddPdf model1("model","pol0+gauss",RooArgList(sig,bkg),RooArgList(sig_yield));
  
  model1.plotOn(xframe,LineColor(kRed));
  
  RooDataHist *data = new RooDataHist("data","dataset with (x)",RooArgList(x),hist);
  //data->plotOn(xframe);
  
  RooPlot* xframe2 = x.frame(Title(("Fit of " + title).c_str()));
  data->plotOn(xframe2, DataError(RooAbsData::SumW2));

  /// Change model here

  
  model1.fitTo(*data);
  model1.plotOn(xframe2);
  model1.paramOn(xframe2,Layout(0.51,0.96,0.918));

  Double_t skew = data->skewness(x); // 3rd standardized moment
  Double_t kurt = data->kurtosis(x); // 4th standardized moment

  TPaveText* tbox = new TPaveText(0.1,0.65,0.5,0.9,"BRNDC");
  tbox->AddText(TString::Format("#splitline{#splitline{#splitline{At energy: %s,}{Fit Function is %s}}{#frac{#chi^2}{ndof} = %.2f & Prob. = %.2f}}{Skewness = %.2f Kurtosis = %.2f}",engy.c_str(),"Gaussian",xframe2->chiSquare(nfreeparams),TMath::Prob(xframe2->chiSquare(nfreeparams)*(nbins-nfreeparams),(nbins-nfreeparams)),skew,kurt));
  tbox->SetFillStyle(0);
  tbox->SetBorderSize(0);
  xframe2->addObject(tbox) ;

  //mean.Print();
  //sigma.Print();
  //sig_yield.Print();
  
  //cout << "chi^2 = " << xframe2->chiSquare() << endl;
  
  RooHist *hresid = xframe2->residHist();
  
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist *hpull = xframe2->pullHist();
 
  RooPlot* xframe3 = x.frame(Title(" "));
  xframe3->addPlotable(hpull, "P");
  xframe3->GetYaxis()->SetTitle("Pull");
  
  TCanvas* c = new TCanvas("fit_canv","fit_canv",800,800);   //tdrDiCanvas("c",hist,hist,8,0);
  c->Divide(1,2);
  plot_frames_on_canvas(c,xframe2,xframe3);

  //delete c;
  
  cout<<"End calculation"<<endl;
  engy = replace_all_substrings(engy," ","_");  c->SaveAs(("plots/energy_" + engy + "_" + title + ".png").c_str());  
  double reutn_values[4] ={double(sigma.getVal()),double(sigma.getError()),mean.getVal(),double(mean.getError())};
  return reutn_values;
}

double* get_resolution_value_modgauss_const(TH1D *hist, string engy, float xlow, float xhigh)
{
  cout<<"Start calculation of hist:"<<hist->GetName()<<" with xlow = "<<xlow<<" and xhigh = "<<xhigh<<endl;
  
  //float xlow = hist->GetBinLowEdge(1);
  int nbins = hist->GetNbinsX();
  int nfreeparams = 4;

  string title = hist->GetTitle();
  
  RooRealVar x("x",title.c_str(),xlow,xhigh);
  RooRealVar mean("mean","mean of gaussian",(xlow+xhigh)/2,xlow,xhigh);
  RooRealVar sigma("sigma","width of gaussian",(xhigh-xlow)/4,(xhigh-xlow)/40,xhigh-xlow);
  RooRealVar N("N","Power of gaussian",2,0.01,10);
  
  RooGenericPdf sig("gauss","exp(-1.0*pow((abs(x-mean)/sigma),N)/N)",RooArgSet(x,mean,sigma,N));  
  
  RooPlot* xframe = x.frame(Title("Gaussian p.d.f."));
  
  mean.setVal((xlow+xhigh)/2);
  sigma.setVal((xhigh-xlow)/4);
  N.setVal(2.5);
  sig.plotOn(xframe);
  
  RooRealVar sig_yield("sig_yield","signal yield",0.5,0.,1/*2*hist->Integral()*/);
  
  RooPolynomial bkg("pol0","pol0",x,RooArgList());
  
  RooAddPdf model1("model","pol0+gauss",RooArgList(sig,bkg),RooArgList(sig_yield));
  
  model1.plotOn(xframe,LineColor(kRed));
  
  RooDataHist *data = new RooDataHist("data","dataset with (x)",RooArgList(x),hist);
  //data->plotOn(xframe);
  
  RooPlot* xframe2 = x.frame(Title(("Fit of " + title).c_str()));
  data->plotOn(xframe2, DataError(RooAbsData::SumW2));

  /// Change model here

  
  model1.fitTo(*data);
  model1.plotOn(xframe2,Components(bkg), LineColor(kRed));
  model1.plotOn(xframe2,Components(sig), LineColor(kGreen));
  model1.plotOn(xframe2);
  model1.paramOn(xframe2,Layout(0.51,0.96,0.918));

  Double_t skew = data->skewness(x); // 3rd standardized moment
  Double_t kurt = data->kurtosis(x); // 4th standardized moment

  TPaveText* tbox = new TPaveText(0.1,0.65,0.5,0.9,"BRNDC");
  tbox->AddText(TString::Format("#splitline{#splitline{#splitline{At energy: %s,}{Fit Function is %s}}{#frac{#chi^2}{ndof} = %.2f & Prob. = %.2f}}{Skewness = %.2f Kurtosis = %.2f}",engy.c_str(),"Modified Gaussian + Const background",xframe2->chiSquare(nfreeparams),TMath::Prob(xframe2->chiSquare(nfreeparams)*(nbins-nfreeparams),(nbins-nfreeparams)),skew,kurt));
  tbox->SetFillStyle(0);
  tbox->SetBorderSize(0);
  xframe2->addObject(tbox) ;

  //mean.Print();
  //sigma.Print();
  //sig_yield.Print();
  
  //cout << "chi^2 = " << xframe2->chiSquare() << endl;
  
  RooHist *hresid = xframe2->residHist();
  
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist *hpull = xframe2->pullHist();
 
  RooPlot* xframe3 = x.frame(Title(" "));
  xframe3->addPlotable(hpull, "P");
  xframe3->GetYaxis()->SetTitle("Pull");
  
  TCanvas* c = new TCanvas("fit_canv","fit_canv",800,800);   //tdrDiCanvas("c",hist,hist,8,0);
  c->Divide(1,2);
  plot_frames_on_canvas(c,xframe2,xframe3);

  //delete c;
  
  cout<<"End calculation"<<endl;
  engy = replace_all_substrings(engy," ","_");  c->SaveAs(("plots/energy_" + engy + "_" + title + ".png").c_str());  
  double reutn_values[8] ={double(sigma.getVal()),double(sigma.getError()),mean.getVal(),double(mean.getError()),double(N.getVal()),double(N.getError()),double(sig_yield.getVal()),double(sig_yield.getError())};
  return reutn_values;
}

double* get_resolution_value_2gauss_diffmean(TH1D *hist, string engy, float xlow, float xhigh)
{
  cout<<"Start calculation of hist:"<<hist->GetName()<<" with xlow = "<<xlow<<" and xhigh = "<<xhigh<<endl;
  
  //float xlow = hist->GetBinLowEdge(1);
  int nbins = hist->GetNbinsX();
  int nfreeparams = 5;

  //float xhigh = hist->GetBinLowEdge(nbins+1);

  string title = hist->GetTitle();
  
  RooRealVar x("x",title.c_str(),xlow,xhigh);
  RooRealVar mean1("mean1","mean of gaussian",(xlow+xhigh)/2,xlow,xhigh);
  RooRealVar sigma1("sigma1","width of gaussian",(xhigh-xlow)/4,(xhigh-xlow)/40,xhigh-xlow);

  RooRealVar mean2("mean2","mean of gaussian",(xlow+xhigh)/2,xlow,xhigh);
  RooRealVar sigma2("sigma2","width of gaussian",(xhigh-xlow)/4,(xhigh-xlow)/40,xhigh-xlow);
  
  RooGaussian gau1("gauss1","gaussian PDF",x,mean1,sigma1);  
  RooGaussian gau2("gauss2","gaussian PDF",x,mean2,sigma2);  
  
  RooPlot* xframe = x.frame(Title("Gaussian p.d.f."));    
  
  mean1.setVal((xlow+xhigh)/2);
  sigma1.setVal((xhigh-xlow)/4);
  mean2.setVal((xlow+xhigh)/3);
  sigma2.setVal((xhigh-xlow)/5);
  // model1.plotOn(xframe);
  
  RooRealVar sig_yield("sig_yield","signal yield",0.5,0.,1/*2*hist->Integral()*/);
  
  RooPolynomial bkg("pol0","pol0",x,RooArgList());
  
  RooAddPdf model1("model","gauss+gauss",RooArgList(gau1,gau2),RooArgList(sig_yield));
  
  model1.plotOn(xframe,LineColor(kRed));
  
  RooDataHist *data = new RooDataHist("data","dataset with (x)",RooArgList(x),hist);
  //data->plotOn(xframe);
  
  RooPlot* xframe2 = x.frame(Title(("Fit of " + title).c_str()));
  data->plotOn(xframe2, DataError(RooAbsData::SumW2));

  /// Change model here

  
  model1.fitTo(*data);
  model1.plotOn(xframe2,Components(gau1), LineColor(kRed));
  model1.plotOn(xframe2,Components(gau2), LineColor(kGreen));
  model1.plotOn(xframe2);
  model1.paramOn(xframe2,Layout(0.55,0.96,0.91));

  Double_t skew = data->skewness(x) ; // 3rd standardized moment
  Double_t kurt = data->kurtosis(x) ; // 4th standardized moment

  TPaveText* tbox = new TPaveText(0.1,0.65,0.5,0.9,"BRNDC");
  tbox->AddText(TString::Format("#splitline{#splitline{#splitline{At energy: %s,}{Fit Function is %s}}{#frac{#chi^2}{ndof} = %.2f & Prob. = %.2f}}{Skewness = %.2f Kurtosis = %.2f}",engy.c_str(),"Sum of 2 Gaussian(with different mean)",xframe2->chiSquare(nfreeparams),TMath::Prob(xframe2->chiSquare(nfreeparams)*(nbins-nfreeparams),(nbins-nfreeparams)),skew,kurt));
  tbox->SetFillStyle(0);
  tbox->SetBorderSize(0);
  xframe2->addObject(tbox) ;

  //mean.Print();
  //sigma.Print();
  //sig_yield.Print();
  
  //cout << "chi^2 = " << xframe2->chiSquare() << endl;
  
  RooHist *hresid = xframe2->residHist();
  
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist *hpull = xframe2->pullHist();
 
  RooPlot* xframe3 = x.frame(Title(" "));
  xframe3->addPlotable(hpull, "P");
  xframe3->GetYaxis()->SetTitle("Pull");
  
  TCanvas* c = new TCanvas("fit_canv","fit_canv",800,800);   //tdrDiCanvas("c",hist,hist,8,0);
  c->Divide(1,2);
  plot_frames_on_canvas(c,xframe2,xframe3);

  //delete c;
  
  cout<<"End calculation"<<endl;
  engy = replace_all_substrings(engy," ","_");  c->SaveAs(("plots/energy_" + engy + "_" + title + ".png").c_str());  
  double return_values[10] ={double(sigma1.getVal()),double(sigma1.getError()),mean1.getVal(),double(mean1.getError()),double(sigma2.getVal()),double(sigma2.getError()),sig_yield.getVal(),double(sig_yield.getError()),mean2.getVal(),double(mean2.getError())}; 
  return return_values;
}

double* get_resolution_value_2gauss_samemean(TH1D *hist, string engy, float xlow, float xhigh)
{
  cout<<"Start calculation of hist:"<<hist->GetName()<<" with xlow = "<<xlow<<" and xhigh = "<<xhigh<<endl;
  
  //float xlow = hist->GetBinLowEdge(1);
  int nbins = hist->GetNbinsX();
  int nfreeparams = 4;

  //float xhigh = hist->GetBinLowEdge(nbins+1);

  string title = hist->GetTitle();
  
  RooRealVar x("x",title.c_str(),xlow,xhigh);
  RooRealVar mean1("mean1","mean of gaussian",(xlow+xhigh)/2,xlow,xhigh);
  RooRealVar sigma1("sigma1","width of gaussian",(xhigh-xlow)/4,(xhigh-xlow)/40,xhigh-xlow);

  RooRealVar mean2("mean2","mean of gaussian",(xlow+xhigh)/2,xlow,xhigh);
  RooRealVar sigma2("sigma2","width of gaussian",(xhigh-xlow)/4,(xhigh-xlow)/40,xhigh-xlow);
  
  RooGaussian gau1("gauss1","gaussian PDF",x,mean1,sigma1);  
  RooGaussian gau2("gauss2","gaussian PDF",x,mean1,sigma2);  
  
  RooPlot* xframe = x.frame(Title("Gaussian p.d.f."));    
  
  mean1.setVal((xlow+xhigh)/2);
  sigma1.setVal((xhigh-xlow)/4);
  mean2.setVal((xlow+xhigh)/3);
  sigma2.setVal((xhigh-xlow)/5);
  // model1.plotOn(xframe);
  
  RooRealVar sig_yield("sig_yield","signal yield",0.5,0.,1/*2*hist->Integral()*/);
  
  RooPolynomial bkg("pol0","pol0",x,RooArgList());
  
  RooAddPdf model1("model","gauss+gauss",RooArgList(gau1,gau2),RooArgList(sig_yield));
  
  model1.plotOn(xframe,LineColor(kRed));
  
  RooDataHist *data = new RooDataHist("data","dataset with (x)",RooArgList(x),hist);
  //data->plotOn(xframe);
  
  RooPlot* xframe2 = x.frame(Title(("Fit of " + title).c_str()));
  data->plotOn(xframe2, DataError(RooAbsData::SumW2));

  /// Change model here

  
  model1.fitTo(*data);
  model1.plotOn(xframe2,Components(gau1), LineColor(kRed));
  model1.plotOn(xframe2,Components(gau2), LineColor(kGreen));
  model1.plotOn(xframe2);
  model1.paramOn(xframe2,Layout(0.55,0.96,0.91));

  Double_t skew = data->skewness(x) ; // 3rd standardized moment
  Double_t kurt = data->kurtosis(x) ; // 4th standardized moment

  TPaveText* tbox = new TPaveText(0.1,0.65,0.5,0.9,"BRNDC");
  tbox->AddText(TString::Format("#splitline{#splitline{#splitline{At energy: %s,}{Fit Function is %s}}{#frac{#chi^2}{ndof} = %.2f & Prob. = %.2f}}{Skewness = %.2f Kurtosis = %.2f}",engy.c_str(),"Sum of 2 Gaussian (with same mean)",xframe2->chiSquare(nfreeparams),TMath::Prob(xframe2->chiSquare(nfreeparams)*(nbins-nfreeparams),(nbins-nfreeparams)),skew,kurt));
  tbox->SetFillStyle(0);
  tbox->SetBorderSize(0);
  xframe2->addObject(tbox) ;

  //mean.Print();
  //sigma.Print();
  //sig_yield.Print();
  
  //cout << "chi^2 = " << xframe2->chiSquare() << endl;
  
  RooHist *hresid = xframe2->residHist();
  
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist *hpull = xframe2->pullHist();
 
  RooPlot* xframe3 = x.frame(Title(" "));
  xframe3->addPlotable(hpull, "P");
  xframe3->GetYaxis()->SetTitle("Pull");
  
  TCanvas* c = new TCanvas("fit_canv","fit_canv",800,800);   //tdrDiCanvas("c",hist,hist,8,0);
  c->Divide(1,2);
  plot_frames_on_canvas(c,xframe2,xframe3);

  //delete c;
  
  cout<<"End calculation"<<endl;
  engy = replace_all_substrings(engy," ","_");  c->SaveAs(("plots/energy_" + engy + "_" + title + ".png").c_str());  
  double return_values[8] ={double(sigma1.getVal()),double(sigma1.getError()),mean1.getVal(),double(mean1.getError()),double(sigma2.getVal()),double(sigma2.getError()),sig_yield.getVal(),double(sig_yield.getError())}; 
  return return_values;
}

double* get_resolution_value_crystalball_const(TH1D *hist, string engy, float xlow, float xhigh)
{
  cout<<"Start calculation of hist:"<<hist->GetName()<<" with xlow = "<<xlow<<" and xhigh = "<<xhigh<<endl;
  
  //float xlow = hist->GetBinLowEdge(1);
  int nbins = hist->GetNbinsX();
  //float xhigh = hist->GetBinLowEdge(nbins+1);

  string title = hist->GetTitle();
  
  RooRealVar x("x",title.c_str(),10000,xlow,xhigh);
  RooRealVar mean("mean","mean of gaussian",(xlow+xhigh)/2,xlow,xhigh);
  RooRealVar sigma("sigma","width of gaussian",(xhigh-xlow)/10,(xhigh - xlow)/40,xhigh - xlow);
  RooRealVar alpha("alpha","alpha of gaussian",0.91,0.001,10);
  RooRealVar n_cb("N","n of gaussian",3,0.001,25);
  
  RooCrystalBall sig("CB","Crystal Ball PDF",x,mean,sigma,alpha,n_cb);  
  
  RooPlot* xframe = x.frame(Title("Gaussian p.d.f."));    
  
  /*mean.setVal((xlow+xhigh)/2);
  sigma.setVal((xhigh-xlow)/10);
  n_cb.setVal(10);
  alpha.setVal(0.91);
  sig.plotOn(xframe);
  */
  //n_cb.setConstant(kTRUE);
  //alpha.setConstant(kTRUE);
  
  RooRealVar sig_yield("sig_yield","signal yield",0.5,0.,1/*2*hist->Integral()*/);
  
  RooPolynomial bkg("pol0","pol0",x,RooArgList());
  
  RooAddPdf model1("model","pol0+gauss",RooArgList(sig,bkg),RooArgList(sig_yield));
  
  model1.plotOn(xframe,LineColor(kRed));
  
  RooDataHist *data = new RooDataHist("data","dataset with (x)",RooArgList(x),hist);
  //data->plotOn(xframe);
  
  RooPlot* xframe2 = x.frame(Title(("Fit of " + title).c_str()));
  data->plotOn(xframe2, DataError(RooAbsData::SumW2));

  /// Change model here

  model1.fitTo(*data,PrintLevel(0));
  model1.plotOn(xframe2,Components(bkg), LineColor(kRed));
  model1.plotOn(xframe2);
    //model1.paramOn(xframe2,data);
  model1.paramOn(xframe2);

  Double_t skew = data->skewness(x) ; // 3rd standardized moment
  Double_t kurt = data->kurtosis(x) ; // 4th standardized moment

  TPaveText* tbox = new TPaveText(0.11,0.55,0.65,0.91,"BRNDC");
  tbox->AddText(TString::Format("#splitline{#splitline{#splitline{At energy: %s,}{Fit Function is %s}}{#frac{#chi^2}{ndof} = %.2f}}{Skewness = %.2f Kurtosis = %.2f}",engy.c_str(),"Single sided Crystal Ball + const background",xframe2->chiSquare(6),skew,kurt));
  tbox->SetFillStyle(0);
  tbox->SetBorderSize(0);
  xframe2->addObject(tbox);

  //mean.Print();
  //sigma.Print();
  //sig_yield.Print();
  
  //cout << "chi^2 = " << xframe2->chiSquare() << endl;
  
  RooHist *hresid = xframe2->residHist();
  
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist *hpull = xframe2->pullHist();
 
  RooPlot* xframe3 = x.frame(Title(" "));
  xframe3->addPlotable(hpull, "P");
  xframe3->GetYaxis()->SetTitle("Pull");
  
  TCanvas* c = new TCanvas("fit_canv","fit_canv",800,800);   //tdrDiCanvas("c",hist,hist,8,0);
  c->Divide(1,2);
  plot_frames_on_canvas(c,xframe2,xframe3);
  
  //delete c;
  
  cout<<"End calculation"<<endl;
  engy = replace_all_substrings(engy," ","_");  c->SaveAs(("plots/energy_" + engy + "_" + title + ".png").c_str());  

  double return_values[10] = {sigma.getVal(), sigma.getError(), mean.getVal(), mean.getError(), alpha.getVal(), alpha.getError(), n_cb.getVal(), n_cb.getError(), sig_yield.getVal(), sig_yield.getError()}; 
  return return_values ;

}

double* get_resolution_value_crystalball(TH1D *hist, string engy, float xlow, float xhigh)
{
  cout<<"Start calculation of hist:"<<hist->GetName()<<" with xlow = "<<xlow<<" and xhigh = "<<xhigh<<endl;
  
  //float xlow = hist->GetBinLowEdge(1);
  int nbins = hist->GetNbinsX();
  int nfreeparams = 4;
  //float xhigh = hist->GetBinLowEdge(nbins+1);

  string title = hist->GetTitle();
  
  RooRealVar x("x",title.c_str(),xlow,xhigh);
  RooRealVar mean("mean","mean of gaussian",(xlow+xhigh)/2,xlow,xhigh);
  RooRealVar sigma("sigma","width of gaussian",(xhigh-xlow)/10,(xhigh - xlow)/nbins,xhigh - xlow);
  RooRealVar alpha("alpha","alpha of gaussian",0.91,0.01,10);
  RooRealVar n_cb("N","n of gaussian",3,0.01,25);
  
  RooCrystalBall sig("CB","Crystal Ball PDF",x,mean,sigma,alpha,n_cb);  
  
  RooPlot* xframe = x.frame(Title("Gaussian p.d.f."));    
  
  mean.setVal((xlow+xhigh)/2);
  sigma.setVal((xhigh-xlow)/10);
  n_cb.setVal(3);
  alpha.setVal(0.91);
  sig.plotOn(xframe);

  //n_cb.setConstant(kTRUE);
  //alpha.setConstant(kTRUE);
  
  RooRealVar sig_yield("sig_yield","signal yield",0.5,0.,1);
  
  RooPolynomial bkg("pol0","pol0",x,RooArgList());
  
  RooAddPdf model1("model","cb",RooArgList(sig),RooArgList());
  
  model1.plotOn(xframe,LineColor(kRed));
  
  RooDataHist *data = new RooDataHist("data","dataset with (x)",RooArgList(x),hist);
  //data->plotOn(xframe);
 
  RooPlot* xframe2 = x.frame(Title(("Fit of " + title).c_str()));
  data->plotOn(xframe2, DataError(RooAbsData::SumW2));

  /// Change model here

  model1.fitTo(*data);
  model1.plotOn(xframe2);
    //model1.paramOn(xframe2,data);
  model1.paramOn(xframe2);

  TPaveText* tbox = new TPaveText(0.1,0.65,0.5,0.9,"BRNDC");
  tbox->AddText(TString::Format("#splitline{#splitline{At energy: %s,}{Fit Function is %s}}{#frac{#chi^2}{ndof} = %.2f & Prob = %.2f}",engy.c_str(),"One sided Crystal Ball",xframe2->chiSquare(nfreeparams),TMath::Prob(xframe2->chiSquare(nfreeparams)*(nbins-nfreeparams),(nbins-nfreeparams))));
  tbox->SetFillStyle(0);
  tbox->SetBorderSize(0);
  xframe2->addObject(tbox);

  //mean.Print();
  //sigma.Print();
  //sig_yield.Print();
  
  //cout << "chi^2 = " << xframe2->chiSquare() << endl;
  
  RooHist *hresid = xframe2->residHist();
  
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist *hpull = xframe2->pullHist();
 
  RooPlot* xframe3 = x.frame(Title(" "));
  xframe3->addPlotable(hpull, "P");
  xframe3->GetYaxis()->SetTitle("Pull");
  
  TCanvas* c = new TCanvas("fit_canv","fit_canv",800,800);   //tdrDiCanvas("c",hist,hist,8,0);
  c->Divide(1,2);
  plot_frames_on_canvas(c,xframe2,xframe3);
  
  //delete c;
  
  cout<<"End calculation"<<endl;
  engy = replace_all_substrings(engy," ","_");  c->SaveAs(("plots/energy_" + engy + "_" + title + ".png").c_str());  

  double return_values[8] = {sigma.getVal(), sigma.getError(), mean.getVal(), mean.getError(), alpha.getVal(), alpha.getError(), n_cb.getVal(), n_cb.getError()}; 

  //unique_ptr<double[]> rv = unique_ptr<double[]>(return_values);
  return return_values ;

  //return rv
}

double* get_resolution_value_crystalball_gauss(TH1D *hist, string engy, float xlow, float xhigh)
{
  cout<<"Start calculation of hist:"<<hist->GetName()<<" with xlow = "<<xlow<<" and xhigh = "<<xhigh<<endl;
  
  //float xlow = hist->GetBinLowEdge(1);
  int nbins = hist->GetNbinsX();
  int nfreeparams = 7;
  //float xhigh = hist->GetBinLowEdge(nbins+1);

  string title = hist->GetTitle();
  
  RooRealVar x("x",title.c_str(),xlow,xhigh);
  RooRealVar mean("mean","mean of gaussian",(xlow+xhigh)/2,xlow,xhigh);
  RooRealVar sigma("sigma","width of gaussian",(xhigh-xlow)/10,(xhigh - xlow)/nbins,xhigh - xlow);
  RooRealVar alpha("alpha","alpha of gaussian",0.91,0.001,10);
  RooRealVar n_cb("N","n of gaussian",3,0.001,20);

  RooRealVar mean_gauss("mean_gauss","mean of gaussian",(xlow+xhigh)/2,xlow,xhigh);
  RooRealVar sigma_gauss("sigma_gauss","width of gaussian",(xhigh-xlow)/10,(xhigh - xlow)/50,xhigh - xlow);
  
  RooGaussian gau("gauss","gaussian PDF",x,mean_gauss,sigma_gauss);  

  RooCrystalBall sig("CB","Crystal Ball PDF",x,mean,sigma,alpha,n_cb);  
  
  RooPlot* xframe = x.frame(Title("Gaussian p.d.f."));    
  
  mean.setVal((xlow+xhigh)/2);
  sigma.setVal((xhigh-xlow)/10);
  n_cb.setVal(3);
  alpha.setVal(0.91);
  sig.plotOn(xframe);

  //n_cb.setConstant(kTRUE);
  //alpha.setConstant(kTRUE);
  
  RooRealVar sig_yield("sig_yield","signal yield",0.5,0.,1);
  
  RooPolynomial bkg("pol0","pol0",x,RooArgList());
  
  RooAddPdf model1("model","cb",RooArgList(sig,gau),RooArgList(sig_yield));
  
  model1.plotOn(xframe,LineColor(kRed));
  
  RooDataHist *data = new RooDataHist("data","dataset with (x)",RooArgList(x),hist);
  //data->plotOn(xframe);
  
  RooPlot* xframe2 = x.frame(Title(("Fit of " + title).c_str()));
  data->plotOn(xframe2, DataError(RooAbsData::SumW2));

  /// Change model here

  model1.fitTo(*data);
  model1.plotOn(xframe2,Components(gau), LineColor(kRed));
  model1.plotOn(xframe2,Components(sig), LineColor(kGreen));
  model1.plotOn(xframe2);
    //model1.paramOn(xframe2,data);
  model1.paramOn(xframe2);

  TPaveText* tbox = new TPaveText(0.1,0.65,0.5,0.9,"BRNDC");
  tbox->AddText(TString::Format("#splitline{#splitline{At energy: %s,}{Fit Function is %s}}{#frac{#chi^2}{ndof} = %.2f & Prob = %.2f}",engy.c_str(),"One sided Crystal Ball + Gaussian",xframe2->chiSquare(nfreeparams),TMath::Prob(xframe2->chiSquare(nfreeparams)*(nbins-nfreeparams),(nbins-nfreeparams))));
  tbox->SetFillStyle(0);
  tbox->SetBorderSize(0);
  xframe2->addObject(tbox);

  //mean.Print();
  //sigma.Print();
  //sig_yield.Print();
  
  //cout << "chi^2 = " << xframe2->chiSquare() << endl;
  
  RooHist *hresid = xframe2->residHist();
  
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist *hpull = xframe2->pullHist();
 
  RooPlot* xframe3 = x.frame(Title(" "));
  xframe3->addPlotable(hpull, "P");
  xframe3->GetYaxis()->SetTitle("Pull");
  
  TCanvas* c = new TCanvas("fit_canv","fit_canv",800,800);   //tdrDiCanvas("c",hist,hist,8,0);
  c->Divide(1,2);
  plot_frames_on_canvas(c,xframe2,xframe3);
  
  //delete c;
  
  cout<<"End calculation"<<endl;
  engy = replace_all_substrings(engy," ","_");  c->SaveAs(("plots/energy_" + engy + "_" + title + ".png").c_str());  

  double return_values[] = {sigma.getVal(), sigma.getError(), mean.getVal(), mean.getError(), alpha.getVal(), alpha.getError(), n_cb.getVal(), n_cb.getError(), sig_yield.getVal(), sig_yield.getError(), mean_gauss.getVal(), mean_gauss.getError(), sigma_gauss.getVal(), sigma_gauss.getError()}; 
  return return_values ;

}

double* get_resolution_value_crystalball_modgauss(TH1D *hist, string engy, float xlow, float xhigh)
{
  cout<<"Start calculation of hist:"<<hist->GetName()<<" with xlow = "<<xlow<<" and xhigh = "<<xhigh<<endl;
  
  //float xlow = hist->GetBinLowEdge(1);
  int nbins = hist->GetNbinsX();
  int nfreeparams = 8;
  //float xhigh = hist->GetBinLowEdge(nbins+1);

  string title = hist->GetTitle();
  
  RooRealVar x("x",title.c_str(),xlow,xhigh);
  RooRealVar mean("mean","mean of gaussian",(xlow+xhigh)/2,xlow,xhigh);
  RooRealVar sigma("sigma","width of gaussian",(xhigh-xlow)/10,(xhigh - xlow)/50,xhigh - xlow);
  RooRealVar alpha("alpha","alpha of gaussian",0.91,0.001,7);
  RooRealVar n_cb("N","n of gaussian",3,0.001,20);

  RooRealVar mean_gauss("mean_gauss","mean of gaussian",(xlow+xhigh)/2,xlow,xhigh);
  RooRealVar sigma_gauss("sigma_gauss","width of gaussian",(xhigh-xlow)/10,(xhigh - xlow)/50,xhigh - xlow);
  RooRealVar N("N_gauss","Power of gaussian",2,0.01,10);
  
  RooGenericPdf gau("gauss","exp(-1.0*pow((abs(x-mean_gauss)/sigma_gauss),N_gauss)/N_gauss)",RooArgSet(x,mean_gauss,sigma_gauss,N));  

  RooCrystalBall sig("CB","Crystal Ball PDF",x,mean,sigma,alpha,n_cb);  
  
  RooPlot* xframe = x.frame(Title("Gaussian p.d.f."));    
  
  mean.setVal((xlow+xhigh)/2);
  sigma.setVal((xhigh-xlow)/10);
  n_cb.setVal(3);
  alpha.setVal(0.91);
  sig.plotOn(xframe);

  //n_cb.setConstant(kTRUE);
  //alpha.setConstant(kTRUE);
  
  RooRealVar sig_yield("sig_yield","signal yield",0.5,0.,1);
  
  RooPolynomial bkg("pol0","pol0",x,RooArgList());
  
  RooAddPdf model1("model","cb",RooArgList(sig,gau),RooArgList(sig_yield));
  
  model1.plotOn(xframe,LineColor(kRed));
  
  RooDataHist *data = new RooDataHist("data","dataset with (x)",RooArgList(x),hist);
  //data->plotOn(xframe);
  
  RooPlot* xframe2 = x.frame(Title(("Fit of " + title).c_str()));
  data->plotOn(xframe2, DataError(RooAbsData::SumW2));

  /// Change model here

  model1.fitTo(*data);
  model1.plotOn(xframe2,Components(gau), LineColor(kRed));
  model1.plotOn(xframe2);
  //model1.paramOn(xframe2,data);
  model1.paramOn(xframe2);
  
  TPaveText* tbox = new TPaveText(0.1,0.65,0.5,0.9,"BRNDC");
  tbox->AddText(TString::Format("#splitline{#splitline{At energy: %s,}{Fit Function is %s}}{#frac{#chi^2}{ndof} = %.2f & Prob = %.2f}",engy.c_str(),"One sided Crystal Ball + Modified Gaussian",xframe2->chiSquare(nfreeparams),TMath::Prob(xframe2->chiSquare(nfreeparams)*(nbins-nfreeparams),(nbins-nfreeparams))));
  tbox->SetFillStyle(0);
  tbox->SetBorderSize(0);
  xframe2->addObject(tbox);

  //mean.Print();
  //sigma.Print();
  //sig_yield.Print();
  
  //cout << "chi^2 = " << xframe2->chiSquare() << endl;
  
  RooHist *hresid = xframe2->residHist();
  
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist *hpull = xframe2->pullHist();
 
  RooPlot* xframe3 = x.frame(Title(" "));
  xframe3->addPlotable(hpull, "P");
  xframe3->GetYaxis()->SetTitle("Pull");
  
  TCanvas* c = new TCanvas("fit_canv","fit_canv",800,800);   //tdrDiCanvas("c",hist,hist,8,0);
  c->Divide(1,2);
  plot_frames_on_canvas(c,xframe2,xframe3);
  
  //delete c;
  
  cout<<"End calculation"<<endl;
  engy = replace_all_substrings(engy," ","_");  c->SaveAs(("plots/energy_" + engy + "_" + title + ".png").c_str());  

  double return_values[] = {sigma.getVal(), sigma.getError(), mean.getVal(), mean.getError(), alpha.getVal(), alpha.getError(), n_cb.getVal(), n_cb.getError(), sig_yield.getVal(), sig_yield.getError(), mean_gauss.getVal(), mean_gauss.getError(), sigma_gauss.getVal(), sigma_gauss.getError()}; 
  return return_values ;

}

double* get_resolution_value_doublecrystalball(TH1D *hist, string engy, float xlow, float xhigh)
{
  cout<<"Start calculation of hist:"<<hist->GetName()<<" with xlow = "<<xlow<<" and xhigh = "<<xhigh<<endl;
  
  //float xlow = hist->GetBinLowEdge(1);
  int nbins = hist->GetNbinsX();
  int nfreeparams = 7;
  //float xhigh = hist->GetBinLowEdge(nbins+1);

  string title = hist->GetTitle();
  
  RooRealVar x("x",title.c_str(),xlow,xlow,xhigh);
  RooRealVar mean("mean","mean of gaussian",(xlow+xhigh)/2,xlow,xhigh);
  RooRealVar sigmaL("sigmaL","width of gaussian",(xhigh - xlow)/10,(xhigh - xlow)/50,xhigh - xlow);
  RooRealVar alphaL("alphaL","alpha of gaussian",0.91,0.001,7);
  RooRealVar nL_cb("NL","n of gaussian",3,0.001,20);
  RooRealVar sigmaR("sigmaR","width of gaussian",(xhigh - xlow)/5,(xhigh - xlow)/50,xhigh - xlow);
  RooRealVar alphaR("alphaR","alpha of gaussian",0.91,0.001,7);
  RooRealVar nR_cb("NR","n of gaussian",3,0.001,20);
  
  RooCrystalBall sig("CB","Crystal Ball PDF",x,mean,sigmaL,sigmaR,alphaL,nL_cb,alphaR,nR_cb);  
  
  RooPlot* xframe = x.frame(Title("Gaussian p.d.f."));    
  
  mean.setVal((xlow+xhigh)/2);
  sigmaL.setVal((xhigh-xlow)/10);
  nL_cb.setVal(3);
  alphaL.setVal(0.91);
  sigmaR.setVal((xhigh-xlow)/10);
  nR_cb.setVal(3);
  alphaR.setVal(0.91);
  sig.plotOn(xframe);
  
  RooRealVar sig_yield("sig_yield","signal yield",0.5,0.,1/*2*hist->Integral()*/);
  
  RooPolynomial bkg("pol0","pol0",x,RooArgList());
  
  RooAddPdf model1("model","cb",RooArgList(sig),RooArgList());
  
  model1.plotOn(xframe,LineColor(kRed));
  
  RooDataHist *data = new RooDataHist("data","dataset with (x)",RooArgList(x),hist);
  //data->plotOn(xframe);
  
  RooPlot* xframe2 = x.frame(Title(("Fit of " + title).c_str()));
  data->plotOn(xframe2, DataError(RooAbsData::SumW2));

  /// Change model here

  model1.fitTo(*data,PrintLevel(0));
  model1.plotOn(xframe2);
  model1.paramOn(xframe2,Layout(0.66,0.96,0.918));
  //model1.paramOn(xframe2,data,"",2,Format("NEU",AutoPrecision(1)));

  TPaveText* tbox = new TPaveText(0.1,0.65,0.5,0.9,"BRNDC");
  tbox->AddText(TString::Format("#splitline{#splitline{At energy: %s,}{Fit Function is %s}}{#frac{#chi^2}{ndof} = %.2f & Prob. = %.2f}",engy.c_str(),"Double sided Crystal Ball",xframe2->chiSquare(nfreeparams),TMath::Prob(xframe2->chiSquare(nfreeparams)*(nbins-nfreeparams),(nbins-nfreeparams))));
  tbox->SetFillStyle(0);
  tbox->SetBorderSize(0);
  xframe2->addObject(tbox) ;

  //mean.Print();
  //sigma.Print();
  //sig_yield.Print();
  
  //cout << "chi^2 = " << xframe2->chiSquare() << endl;
  
  RooHist *hresid = xframe2->residHist();
  
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist *hpull = xframe2->pullHist();
 
  RooPlot* xframe3 = x.frame(Title(" "));
  xframe3->addPlotable(hpull, "P");
  xframe3->GetYaxis()->SetTitle("Pull");
  
  TCanvas* c = new TCanvas("fit_canv","fit_canv",800,800);   //tdrDiCanvas("c",hist,hist,8,0);
  c->Divide(1,2);

  plot_frames_on_canvas(c,xframe2,xframe3);
  
  //delete c;
  
  cout<<"End calculation"<<endl;
  engy = replace_all_substrings(engy," ","_");  c->SaveAs(("plots/energy_" + engy + "_" + title + ".png").c_str());  

  double return_values[14] = {double(sigmaR.getVal()), sigmaR.getError(), mean.getVal(), mean.getError(), sigmaL.getVal(), sigmaL.getError(), alphaL.getVal(), alphaL.getError(), nL_cb.getVal(), nL_cb.getError(), alphaR.getVal(), alphaR.getError(), nR_cb.getVal(), nR_cb.getError()}; 
  return return_values ;
}


double* get_resolution_value_doublecrystalball_const(TH1D *hist, string engy, float xlow, float xhigh)
{
  cout<<"Start calculation of hist:"<<hist->GetName()<<" with xlow = "<<xlow<<" and xhigh = "<<xhigh<<endl;
  
  //float xlow = hist->GetBinLowEdge(1);
  int nbins = hist->GetNbinsX();
  //float xhigh = hist->GetBinLowEdge(nbins+1);

  string title = hist->GetTitle();
  
  RooRealVar x("x",title.c_str(),xlow,xlow,xhigh);
  RooRealVar mean("mean","mean of gaussian",100,xlow,xhigh);
  RooRealVar sigmaL("sigmaL","width of gaussian",100,0.00001,xhigh - xlow);
  RooRealVar alphaL("alphaL","alpha of gaussian",10,0.00001,20.);
  RooRealVar nL_cb("NL","n of gaussian",10,0.0001,18);
  RooRealVar sigmaR("sigmaR","width of gaussian",10,0.00001,xhigh - xlow);
  RooRealVar alphaR("alphaR","alpha of gaussian",10,0.00001,20.);
  RooRealVar nR_cb("NR","n of gaussian",1,0.00001,18);
  
  RooCrystalBall sig("CB","Crystal Ball PDF",x,mean,sigmaL,sigmaR,alphaL,nL_cb,alphaR,nR_cb);  
  
  RooPlot* xframe = x.frame(Title("Gaussian p.d.f."));    
  
  mean.setVal((xlow+xhigh)/2);
  sigmaL.setVal((xhigh-xlow)/10);
  nL_cb.setVal(3);
  alphaL.setVal(0.91);
  sigmaR.setVal((xhigh-xlow)/10);
  nR_cb.setVal(3);
  alphaR.setVal(0.91);
  sig.plotOn(xframe);
  
  RooRealVar sig_yield("sig_yield","signal yield",0.5,0.,1/*2*hist->Integral()*/); 
  sig_yield.setVal(0.5);
 
  RooRealVar p1("p1","coeff #1", 10000, -1000., 1000.);
  
  RooPolynomial bkg("pol0","pol0",x,RooArgList());
  
  RooAddPdf model1("model","pol0+gauss",RooArgList(sig,bkg),RooArgList(sig_yield));

  model1.plotOn(xframe,LineColor(kRed));
  
  RooDataHist *data = new RooDataHist("data","dataset with (x)",RooArgList(x),hist);
  //data->plotOn(xframe);
  
  RooPlot* xframe2 = x.frame(Title(("Fit of " + title).c_str()));
  data->plotOn(xframe2, DataError(RooAbsData::SumW2));

  /// Change model here

  model1.fitTo(*data,PrintLevel(0));
  model1.plotOn(xframe2,Components(bkg), LineColor(kRed));
  model1.plotOn(xframe2);
  model1.paramOn(xframe2,Layout(0.68,0.98,0.918));
  //model1.paramOn(xframe2,data,"",2,Format("NEU",AutoPrecision(1)));
  
  // Skewness and kurtosis
  Double_t skew = data->skewness(x) ; // 3rd standardized moment
  Double_t kurt = data->kurtosis(x) ; // 4th standardized moment

  TPaveText* tbox = new TPaveText(0.11,0.55,0.65,0.91,"BRNDC");
  tbox->AddText(TString::Format("#splitline{#splitline{#splitline{At energy: %s,}{Fit Function is %s}}{#frac{#chi^2}{ndof} = %.2f}}{Skewness = %f Kurtosis = %f}",engy.c_str(),"Double sided Crystal Ball + constant background",xframe2->chiSquare(9),skew,kurt));
  tbox->SetFillStyle(0);
  tbox->SetBorderSize(0);
  xframe2->addObject(tbox) ;

  //mean.Print();
  //sigma.Print();
  //sig_yield.Print();
  
  //cout << "chi^2 = " << xframe2->chiSquare() << endl;
  
  RooHist *hresid = xframe2->residHist();
  
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist *hpull = xframe2->pullHist();
 
  RooPlot* xframe3 = x.frame(Title(" "));
  xframe3->addPlotable(hpull, "P");
  xframe3->GetYaxis()->SetTitle("Pull");
  
  TCanvas* c = new TCanvas("fit_canv","fit_canv",800,800);   //tdrDiCanvas("c",hist,hist,8,0);
  c->Divide(1,2);

  plot_frames_on_canvas(c,xframe2,xframe3);
  
  //delete c;
  
  cout<<"End calculation"<<endl;
  engy = replace_all_substrings(engy," ","_");  c->SaveAs(("plots/energy_" + engy + "_" + title + ".png").c_str()); 

  double return_values[16] = {double(sigmaR.getVal()), sigmaR.getError(), mean.getVal(), mean.getError(), sigmaL.getVal(), sigmaL.getError(), alphaL.getVal(), alphaL.getError(), nL_cb.getVal(), nL_cb.getError(), alphaR.getVal(), alphaR.getError(), nR_cb.getVal(), nR_cb.getError(), sig_yield.getVal(), sig_yield.getError()}; 
  return return_values ;

}

double* get_resolution_simfit(TH1D *hist_MCP1,TH1D *hist_MCP2, string engy)
 {
   cout<<"Start calculation of hist:"<<hist_MCP1->GetName()<<" and hist:"<<hist_MCP2->GetName()<<endl;

   //float xlow = hist->GetBinLowEdge(1);
   int nbins1 = hist_MCP1->GetNbinsX();
   int nbins2 = hist_MCP2->GetNbinsX();
   int nbins = nbins1 + nbins2;
   int nfreeparams = 7;
   
   float xlow1 = hist_MCP1->GetXaxis()->GetBinLowEdge(1);
   float xhigh1 = hist_MCP1->GetXaxis()->GetBinLowEdge(nbins1) + hist_MCP1->GetXaxis()->GetBinWidth(nbins1);
   float xlow2 = hist_MCP2->GetXaxis()->GetBinLowEdge(1);
   float xhigh2 = hist_MCP2->GetXaxis()->GetBinLowEdge(nbins2) + hist_MCP1->GetXaxis()->GetBinWidth(nbins2);

   float xlow = min(xlow1,xlow2);
   float xhigh = max(xhigh1,xhigh2);
   //float xhigh = hist->GetBinLowEdge(nbins+1);

   string title1 = hist_MCP1->GetTitle();
   string title2 = hist_MCP2->GetTitle();
  
   RooRealVar x1("x1",title1.c_str(),xlow1,xhigh1);
   RooRealVar x2("x2",title2.c_str(),xlow2,xhigh2);

   RooRealVar x("x","total range",min(xlow1,xlow2),max(xhigh1,xhigh2));
   
   RooRealVar mean1("mean1","mean of gaussian",(xlow1+xhigh1)/2.0,xlow1,xhigh1);
   RooRealVar sigma1("sigma1","width of gaussian",(xhigh-xlow)/4.0,(xhigh-xlow)/40,xhigh-xlow);

   RooRealVar mean2("mean2","mean of gaussian",(xlow2+xhigh2)/2.0,xlow2,xhigh2);
   RooRealVar sigma2("sigma2","width of gaussian",(xhigh-xlow)/4.0,(xhigh-xlow)/40,xhigh-xlow);
  
   RooRealVar mean3("mean3","mean of gaussian",(xlow1+xhigh1)/2.0,xlow1,xhigh1);
   RooRealVar mean4("mean4","mean of gaussian",(xlow2+xhigh2)/2.0,xlow2,xhigh2);

   RooGaussian gau1("gauss1","gaussian PDF",x1,mean1,sigma1);  
   RooGaussian gau2("gauss2","gaussian PDF",x1,mean2,sigma2);  
   RooGaussian gau3("gauss3","gaussian PDF",x2,mean3,sigma1);  
   RooGaussian gau4("gauss4","gaussian PDF",x2,mean4,sigma2);  
    
   mean1.setVal((xlow1+xhigh1)/2.0);
   sigma1.setVal((xhigh-xlow)/4.0);
   mean2.setVal((xlow2+xhigh2)/2.0);
   sigma2.setVal((xhigh-xlow)/5.0);
   mean3.setVal((xlow1+xhigh1)/2.0);
   mean4.setVal((xlow2+xhigh2)/2.0);
   // model1.plotOn(xframe);
  
   RooRealVar sig_yield1("sig_yield1","signal yield",0.5,0.,1/*2*hist->Integral()*/);
   RooRealVar sig_yield2("sig_yield2","signal yield",0.5,0.,1/*2*hist->Integral()*/);
  
   RooPolynomial bkg("pol0","pol0",x1,RooArgList());
  
   RooAddPdf model1("model1","gauss+gauss",RooArgList(gau1,gau2),RooArgList(sig_yield1));

   RooAddPdf model2("model2","gauss+gauss",RooArgList(gau3,gau4),RooArgList(sig_yield1));

   RooDataHist *data_MCP1 = new RooDataHist("data_MCP1","dataset with (x1)",RooArgList(x),hist_MCP1);
   RooDataHist *data_MCP2 = new RooDataHist("data_MCP2","dataset with (x2)",RooArgList(x),hist_MCP2);

   //RooDataSet *data_MCP1 = model1.generate(RooArgSet(x1),100);
   //RooDataSet *data_MCP2 = model2.generate(RooArgSet(x2),100) ;
 
   RooCategory sample("sample","sample");
   sample.defineType("MCP1");
   sample.defineType("MCP2");

   cout<<"check1"<<endl;
   RooDataHist combData("combData","combined data",x,Index(sample),Import("MCP1",*data_MCP1),Import("MCP2",*data_MCP2));
   cout<<"check2"<<endl;
   
   RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
  
   // Associate model with the physics state and model_ctl with the control state
   simPdf.addPdf(model1,"MCP1") ;
   simPdf.addPdf(model2,"MCP2") ;

   simPdf.fitTo(combData);

   RooPlot* xframe1 = x1.frame(Title(("Fit of " + title1).c_str()));
   combData.plotOn(xframe1,Cut("sample==sample::MCP1"));
   simPdf.plotOn(xframe1,Slice(sample,"MCP1"),Components("gau1"),ProjWData(sample,combData),LineColor(kRed),LineStyle(kDashed));  
   simPdf.plotOn(xframe1,Slice(sample,"MCP1"),Components("gau2"),ProjWData(sample,combData),LineColor(kGreen),LineStyle(kDashed));  
   simPdf.plotOn(xframe1,Slice(sample,"MCP1"),ProjWData(sample,combData));
   simPdf.paramOn(xframe1,Slice(sample,"MCP1"),ProjWData(sample,combData),Layout(0.55,0.96,0.91));
  
   Double_t skew1 = data_MCP1->skewness(x1) ; // 3rd standardized moment
   Double_t kurt1 = data_MCP1->kurtosis(x1) ; // 4th standardized moment

   TPaveText* tbox = new TPaveText(0.1,0.65,0.5,0.9,"BRNDC");
   tbox->AddText(TString::Format("#splitline{#splitline{#splitline{At energy: %s,}{Simultaneius Fit Function is %s}}{#frac{#chi^2}{ndof} = %.2f & Prob. = %.2f}}{Skewness = %.2f Kurtosis = %.2f}",engy.c_str(),"Sum of 2 Gaussian(with different mean)",xframe1->chiSquare(nfreeparams),TMath::Prob(xframe1->chiSquare(nfreeparams)*(nbins-nfreeparams),(nbins-nfreeparams)),skew1,kurt1));
   tbox->SetFillStyle(0);
   tbox->SetBorderSize(0);
   xframe1->addObject(tbox) ;
  
   RooHist *hresid1 = xframe1->residHist();
   RooHist *hpull1 = xframe1->pullHist();
 
   RooPlot* xframe3 = x1.frame(Title(" "));
   xframe3->addPlotable(hpull1, "P");
   xframe3->GetYaxis()->SetTitle("Pull");
  
   TCanvas* c1 = new TCanvas("fit_canv1","fit_canv1",800,800);   //tdrDiCanvas("c",hist,hist,8,0);
   c1->Divide(1,2);
   plot_frames_on_canvas(c1,xframe1,xframe3);

   
   RooPlot* xframe2 = x2.frame(Title(("Fit of " + title1).c_str()));
   combData.plotOn(xframe2,Cut("sample==sample::MCP2"));
   simPdf.plotOn(xframe2,Slice(sample,"MCP2"),Components("gau3"),ProjWData(sample,combData),LineColor(kRed),LineStyle(kDashed));  
   simPdf.plotOn(xframe2,Slice(sample,"MCP2"),Components("gau4"),ProjWData(sample,combData),LineColor(kGreen),LineStyle(kDashed));  
   simPdf.plotOn(xframe2,Slice(sample,"MCP2"),ProjWData(sample,combData));
   simPdf.paramOn(xframe2,Slice(sample,"MCP2"),ProjWData(sample,combData),Layout(0.55,0.96,0.91));
  
   Double_t skew2 = data_MCP2->skewness(x2) ; // 3rd standardized moment
   Double_t kurt2 = data_MCP2->kurtosis(x2) ; // 4th standardized moment

   TPaveText* tbox2 = new TPaveText(0.1,0.65,0.5,0.9,"BRNDC");
   tbox2->AddText(TString::Format("#splitline{#splitline{#splitline{At energy: %s,}{Simultaneius Fit Function is %s}}{#frac{#chi^2}{ndof} = %.2f & Prob. = %.2f}}{Skewness = %.2f Kurtosis = %.2f}",engy.c_str(),"Sum of 2 Gaussian(with different mean)",xframe2->chiSquare(nfreeparams),TMath::Prob(xframe2->chiSquare(nfreeparams)*(nbins-nfreeparams),(nbins-nfreeparams)),skew2,kurt2));
   tbox2->SetFillStyle(0);
   tbox2->SetBorderSize(0);
   xframe2->addObject(tbox2) ;
  
   RooHist *hresid2 = xframe2->residHist();
   RooHist *hpull2 = xframe2->pullHist();
 
   RooPlot* xframe4 = x2.frame(Title(" "));
   xframe4->addPlotable(hpull2, "P");
   xframe4->GetYaxis()->SetTitle("Pull");
  
   TCanvas* c2 = new TCanvas("fit_canv2","fit_canv2",800,800);   //tdrDiCanvas("c",hist,hist,8,0);
   c2->Divide(1,2);
   plot_frames_on_canvas(c2,xframe2,xframe4);

  
   cout<<"End calculation"<<endl;
   engy = replace_all_substrings(engy," ","_");
   c1->SaveAs(("plots/energy_" + engy + "_simfit_" + title1 + ".png").c_str());  
   c2->SaveAs(("plots/energy_" + engy + "_simfit_" + title2 + ".png").c_str());  

   double return_values[10] ={double(sigma1.getVal()),double(sigma1.getError()),mean1.getVal(),double(mean1.getError()),double(sigma2.getVal()),double(sigma2.getError()),sig_yield1.getVal(),double(sig_yield1.getError()),mean2.getVal(),double(mean2.getError())}; 
   return return_values;
 }


double get_mean_value_fromfit(TH1D *hist, string engy, float xlow, float xhigh)
   {
     cout<<"Start calculation of hist:"<<hist->GetName()<<" with xlow = "<<xlow<<" and xhigh = "<<xhigh<<endl;
  
     //float xlow = hist->GetBinLowEdge(1);
     int nbins = hist->GetNbinsX();
     //float xhigh = hist->GetBinLowEdge(nbins+1);

     string title = hist->GetTitle();
  
     RooRealVar x("x",title.c_str(),xlow,xlow,xhigh);
     RooRealVar mean("mean","mean of gaussian",(xlow+xhigh)/2,xlow,xhigh);
     RooRealVar sigma("sigma","width of gaussian",1,(xhigh - xlow)/50,xhigh - xlow);
     RooRealVar alpha("alpha","alpha of gaussian",10,0.00001,30.);
     RooRealVar n_cb("N","n of gaussian",10,0.0001,18);
  
     RooCrystalBall sig("CB","Crystal Ball PDF",x,mean,sigma,alpha,n_cb);  
  
     RooPlot* xframe = x.frame(Title("Gaussian p.d.f."));    
  
     mean.setVal((xlow+xhigh)/2);
     sigma.setVal((xhigh-xlow)/10);
     n_cb.setVal(3);
     alpha.setVal(0.91);
     sig.plotOn(xframe);

     //n_cb.setConstant(kTRUE);
     //alpha.setConstant(kTRUE);
  
     RooRealVar sig_yield("sig_yield","signal yield",0.5,0.,1);
  
     RooPolynomial bkg("pol0","pol0",x,RooArgList());
  
     RooAddPdf model1("model","cb",RooArgList(sig),RooArgList());
  
     model1.plotOn(xframe,LineColor(kRed));
  
     RooDataHist *data = new RooDataHist("data","dataset with (x)",RooArgList(x),hist);
     //data->plotOn(xframe);
  
     RooPlot* xframe2 = x.frame(Title(("Fit of " + title).c_str()));
     data->plotOn(xframe2, DataError(RooAbsData::SumW2));

     /// Change model here

     model1.fitTo(*data,PrintLevel(0));
     model1.plotOn(xframe2);
       //model1.paramOn(xframe2,data);
  model1.paramOn(xframe2);

     TPaveText* tbox = new TPaveText(0.1,0.65,0.5,0.9,"BRNDC");
     tbox->AddText(TString::Format("#splitline{#splitline{At energy: %s,}{Fit Function is %s}}{#frac{#chi^2}{ndof} = %.2f}",engy.c_str(),"One sided Crystal Ball",xframe2->chiSquare(4)));
     tbox->SetFillStyle(0);
     tbox->SetBorderSize(0);
     xframe2->addObject(tbox);

     //mean.Print();
     //sigma.Print();
     //sig_yield.Print();
  
     //cout << "chi^2 = " << xframe2->chiSquare() << endl;
  
     RooHist *hresid = xframe2->residHist();
  
     // Construct a histogram with the pulls of the data w.r.t the curve
     RooHist *hpull = xframe2->pullHist();
 
     RooPlot* xframe3 = x.frame(Title(" "));
     xframe3->addPlotable(hpull, "P");
     xframe3->GetYaxis()->SetTitle("Pull");
  
     TCanvas* c = new TCanvas("fit_canv","fit_canv",800,800);   //tdrDiCanvas("c",hist,hist,8,0);
     c->Divide(1,2);
     plot_frames_on_canvas(c,xframe2,xframe3);
  
     //delete c;
  
     cout<<"End calculation"<<endl;
     engy = replace_all_substrings(engy," ","_");  c->SaveAs(("plots/energy_" + engy + "_" + title + ".png").c_str());  
     return double(mean.getVal());
}


  void gethist(string fileinname, string treename, string fillhist, string weight_sel, float weight, TH1D *hist_out)
{
    TFile *filein;
    filein = new TFile(fileinname.c_str(),"read");
    fillhist = fillhist + ">>+" + "tmp";
    if(!filein->IsOpen())
    {cout<<"Check: file "<<fileinname<<" not present"<<endl;exit(0);}
    TTree *tree = (TTree*)filein->Get(treename.c_str());
    if(tree == NULL)
    {cout<<"Check: file "<<fileinname<<" does not have tree "<<treename<<endl;exit(0);}
    TH1D *tmp = (TH1D*)hist_out->Clone("tmp");
    tmp->Reset("ICES");
    cout<<fillhist<<" , "<<weight_sel<<endl;
    tree->Draw(fillhist.c_str(),weight_sel.c_str());
    //tree->Draw(fillhist.c_str());
        //cout<<fileinname<<" "<<treename<<" "<<fillhist<<" , "<<hist_basic->GetNbinsX()<<"  "<<hist_basic->GetEntries();
    //TH1D *tmp = (TH1D*)gDirectory->Get(hist_out->GetName());
    //for(int i =1; i<tmp->GetNbinsX()+1; i++)
    //cout<<"bin "<<i<<" "<<tmp->GetBinLowEdge(i)<<" to "<<tmp->GetBinLowEdge(i)+tmp->GetBinWidth(i)<<" "<<tmp->GetBinContent(i)<<endl;
    if(tmp != NULL)
    {
        //tmp->Sumw2();
        tmp->Scale(weight);
        hist_out->Add(tmp);
        //cout<<fileinname<<" "<<treename<<" "<<fillhist<<" , "<<hist_out->GetName()<<" "<<tmp->GetNbinsX()<<"  "<<tmp->GetEntries();
    }
    filein->Close();
}

void get2dhist(string fileinname, string treename, string fillhist, string weight_sel, float weight, TH2D *hist_out)
{
  TFile *filein2;
  filein2 = new TFile(fileinname.c_str(),"read");
  if(!filein2->IsOpen())
    {cout<<"Check: file "<<fileinname<<" not present"<<endl;exit(0);}
  filein2->cd();
  TTree *tree = (TTree*)filein2->Get(treename.c_str());
  fillhist = fillhist + ">>+" + "tmp";
   if(tree == NULL)
    {cout<<"Check: file "<<fileinname<<" does not have tree "<<treename<<endl;exit(0);}
  //tree->Draw(fillhist.c_str(),weight_sel.c_str());
  //cout<<fileinname<<" "<<treename<<" "<<fillhist<<" , "<<hist_basic->GetNbinsX()<<"  "<<hist_basic->GetEntries();
  //TH2D *tmp = (TH2D*)gDirectory->Get(hist_out->GetName());
  TH2D *tmp = (TH2D*)hist_out->Clone("tmp");
  tmp->Reset("ICES");
  cout<<fillhist<<" , "<<weight_sel<<" "<<filein2<<endl;
  tree->Draw(fillhist.c_str(),weight_sel.c_str());
  if(tmp != NULL)
    {
      tmp->Scale(weight);
      hist_out->Add(tmp);
      //cout<<fileinname<<" "<<treename<<" "<<fillhist<<" , "<<hist_out->GetName()<<" "<<tmp->GetNbinsX()<<"  "<<tmp->GetEntries();
      cout<<fillhist<<" , "<<weight_sel<<" "<<filein2<<endl;
    }
  if(filein2!= NULL && filein2->IsOpen()) filein2->Close();
}

void plot_time_resolution()
{

  string central_crystal;
#ifdef target_C2
  central_crystal = "C2";
#elif defined(target_C3)
  central_crystal = "C3";
#endif

  ///////////////// Plot basic histogram to study
  
  string A_eff_C2_MCP1 = "1.0 / sqrt( pow(b_rms[" + central_crystal  + "]/amp_max[" + central_crystal + "],2) + pow(b_rms[MCP1]/amp_max[MCP1],2))"; 
  string A_eff_C2_MCP2 = "1.0 / sqrt( pow(b_rms[" + central_crystal + "]/amp_max[" + central_crystal + "],2) + pow(b_rms[MCP2]/amp_max[MCP2],2))"; 
  string A_eff_MCP1_MCP2 = "1.0 / sqrt( pow(b_rms[MCP2]/amp_max[MCP2],2) + pow(b_rms[MCP1]/amp_max[MCP1],2))"; 
  string A_eff_C2_C3 = "1.0 / sqrt( pow(b_rms[C2]/amp_max[C2],2) + pow(b_rms[C3]/amp_max[C3],2))"; 

  //string basicvar[8] = {"amp_max["+central_crystal+"]/b_rms["+central_crystal+"]","amp_max[MCP1]/b_rms[MCP1]","amp_max[MCP2]/b_rms[MCP2]","amp_max[C3]/b_rms[b_3]",A_eff_C2_MCP1,A_eff_C2_MCP2,A_eff_MCP1_MCP2,A_eff_C2_C3};
  //string basicvarname[8] = {"#frac{A}{#sigma_{n}} of " + central_crystal,"#frac{A}{#sigma_{n}} of MCP1","#frac{A}{#sigma_{n}} of MCP2","#frac{A}{#sigma_{n}} of C3","#frac{A_{eff}}{#sigma_{n}} of "+central_crystal+" and MCP1","#frac{A_{eff}}{#sigma_{n}} of "+central_crystal+" and MCP2","#frac{A_{eff}}{#sigma_{n}} of MCP1 and MCP2","#frac{A_{eff}}{#sigma_{n}} of C2 and C3"};

  string basicvar[8] = {"amp_max["+central_crystal+"]/b_rms["+central_crystal+"]","amp_max[MCP1]/b_rms[MCP1]","amp_max[MCP2]/b_rms[MCP2]","amp_max[C3]/b_rms[b_3]","amp_max["+central_crystal+"]/b_rms["+central_crystal+"]","amp_max["+central_crystal+"]/b_rms["+central_crystal+"]",A_eff_MCP1_MCP2,A_eff_C2_C3};
  string basicvarname[8] = {"#frac{A}{#sigma_{n}} of " + central_crystal,"#frac{A}{#sigma_{n}} of MCP1","#frac{A}{#sigma_{n}} of MCP2","#frac{A}{#sigma_{n}} of C3","#frac{A_{eff}}{#sigma_{n}} of "+central_crystal,"#frac{A_{eff}}{#sigma_{n}} of "+central_crystal,"#frac{A_{eff}}{#sigma_{n}} of MCP1 and MCP2","#frac{A_{eff}}{#sigma_{n}} of C2 and C3"};

  string basicfitvar[8] = {"dt_"+central_crystal+"_MCP1","dt_"+central_crystal+"_MCP2","dt_MCP1_MCP2","fit_time_C4","dt_"+central_crystal+"_MCP1","dt_"+central_crystal+"_MCP2","dt_MCP1_MCP2","dt_C2_C3_phase_corrected"};
  string basicfitvarname[8] = {"dt_"+central_crystal+"_MCP1","dt_"+central_crystal+"_MCP2","dt_MCP1_MCP2","t_C4","dt_"+central_crystal+"_MCP1","dt_"+central_crystal+"_MCP2","dt_MCP1_MCP2","dt_C2_C3"};

  string basic_cut[8] = {"amp_max["+central_crystal+"] > 400","amp_max[MCP1]> 100","amp_max[MCP2]>100","amp_max[C4] > 20","amp_max["+central_crystal+"] > 400 && amp_max[MCP1]> 100","amp_max["+central_crystal+"] > 400 && amp_max[MCP2]> 100","amp_max[MCP1] > 100 && amp_max[MCP2]>100","amp_max[C2] > 400 && amp_max[C3] > 20"}; //// Always check

  const int nbins_aeff = 7;
  //float Aeff_bin_edges[nbins_aeff+1] = {0., 100., 150., 200., 300., 400., 600.0, 10000.};
  //float Aeff_bin_edges[nbins_aeff+1] = {0., 100.0, 160., 300., 700};
  //float Aeff_bin_edges[nbins_aeff+1] = {0., 100.0, 130, 200, 300, 450, 900};
  float Aeff_bin_edges[nbins_aeff+1] = { 0., 100.0,  200, 300, 400, 550, 1000, 1500};
  double bin_low[8] = {0,0,0,0,0,0,0,0};
  double bin_up[8] = {3000,700,700,500,1000,1000,700,1000};
  int nbin[8] = {20,20,20,20,20,20,20,20};

  string aeff_gain_cutC2 = " && ( (amp_max["+central_crystal+"]/b_rms["+central_crystal+"] > 1500 && gain[C2] == 10 ) ||  (amp_max["+central_crystal+"]/b_rms[C2] < 1500 && gain["+central_crystal+"] == 1 ) )";

  string time_res_title[4] = {"Combined time resolution of C3 and MCP1 (in ps)","Combined time resolution of C3 and MCP2 (in ps)","Time resolution of MCP (in ps)","Time resolution of crsytal (in ps)"};


  /////////// Compute 1d hist
  string fillvar_1d[] = {"amp_max[MCP1]","amp_max[MCP2]","amp_max["+central_crystal+"]",A_eff_MCP1_MCP2,"amp_max["+central_crystal+"]/b_rms["+central_crystal+"]","b_rms["+central_crystal+"]"};
  string fillvar_1d_name[] = {"amp_max[MCP1]","amp_max[MCP2]","amp_max["+central_crystal+"]","A_eff_MCP1_MCP2","A_eff_"+central_crystal,"b_rms["+central_crystal+"]"};

  int nvars_1d = sizeof(fillvar_1d)/sizeof(fillvar_1d[0]);

  double bin_low_1d[] = {0,0,0,0,0,0};
  double bin_up_1d[] = {1000,1000,6000,1000,1500,8};
  int nbin_1d[] = {30,30,30,30,30,30};

  if(nvars_1d != sizeof(bin_low_1d)/sizeof(bin_low_1d[0]) || nvars_1d != sizeof(bin_up_1d)/sizeof(bin_up_1d[0]) || nvars_1d != sizeof(nbin_1d)/sizeof(nbin_1d[0]))
    {cout<<"No of variable and bin information dont match for 1d hist for time resolution"<<endl;exit(0);}
  for(int ig = 10; ig <2; ig++)
    {
      for(int ie =0; ie < nenergies; ie++)
	{
	  for(int iv =0; iv < nvars_1d; iv++)
	    {
	      TH1D* hist_1d;
	      hist_1d = new TH1D(("hist_nocuts_" + fillvar_1d_name[iv] + energies[ie]).c_str(),("hist" + fillvar_1d_name[iv]).c_str(),nbin_1d[iv],bin_low_1d[iv],bin_up_1d[iv]);
	      gethist(energies[ie] + "_" + central_crystal + "2x2_"+file_tag+".root",treename,fillvar_1d[iv],"",1,hist_1d);
	  
	      string str = fillvar_1d_name[iv];
	      plot_hist(hist_1d,str,"energy_" + energies[ie] + central_crystal + "_" + "_no_cuts_" + fillvar_1d_name[iv] + ".png");

	      string evt_sel_cut = "n_h1X > 0 && n_h1Y>0 && h1X > " + to_string(max_posx[ie]-3)  + " && h1X < " + to_string(max_posx[ie]+3) + " && h1Y > " + to_string(max_posy[ie]-3)  + " && h1Y < " + to_string(max_posy[ie]+3) + " &&  gain["+central_crystal+"] == " + gain_cut[ig] + " && amp_max[MCP1] > 200 && amp_max[MCP2] > 200 && " + amp_max_cut[ie];

	      TH1D* hist_1d_cuts;
	      hist_1d_cuts = new TH1D(("hist_nocuts_" + fillvar_1d_name[iv] + energies[ie]).c_str(),("hist" + fillvar_1d_name[iv]).c_str(),nbin_1d[iv],bin_low_1d[iv],bin_up_1d[iv]);
	      gethist(energies[ie] + "_" + central_crystal + "2x2_"+file_tag+".root",treename,fillvar_1d[iv],evt_sel_cut,1,hist_1d_cuts);
	  
	      str = fillvar_1d_name[iv];
	      plot_hist(hist_1d_cuts,str,"energy_" + energies[ie] + central_crystal + "_" + "_with_eventsel_gain" + gain_cut[ig] + "_" + fillvar_1d_name[iv] + ".png");


	    }
	}
    }
 
  //////////// Compute 2D hist

  string energy_weighted_posX = "Sum$((max(8,amp_max_crystal[])-8)*posX_crystal[])/Sum$((max(8,amp_max_crystal[])-8))";
  string energy_weighted_posY = "Sum$((max(8,amp_max_crystal[])-8)*posY_crystal[])/Sum$((max(8,amp_max_crystal[])-8))";
  string logenergy_weighted_posX = "Sum$(log((max(8,amp_max_crystal[])-7))*posX_crystal[])/Sum$(log((max(8,amp_max_crystal[])-7)))";
  string logenergy_weighted_posY = "Sum$(log((max(8,amp_max_crystal[])-7))*posY_crystal[])/Sum$(log((max(8,amp_max_crystal[])-7)))";

  string fillvar1_2d[] = {"h1X","h2X","h1X","h1Y","fit_time["+central_crystal+"]",A_eff_C2_MCP1,A_eff_C2_MCP2,A_eff_MCP1_MCP2,"h1X","h1Y","h1X","h1Y"};
  string fillvar1_2d_name[] = {"h1X","h2X","h1X","h1Y","fit_time["+central_crystal+"]","A_eff_C2_MCP1","A_eff_C2_MCP2","A_eff_MCP1_MCP2","h1X","h1Y","h1X","h1Y"};
  string fillvar2_2d[] = {"h1Y","h2Y","amp_max[" + central_crystal +"]","amp_max[" + central_crystal +"]","dt_"+central_crystal+"_MCP1","dt_"+central_crystal+"_MCP1","dt_"+central_crystal+"_MCP2","dt_MCP1_MCP2",energy_weighted_posX,energy_weighted_posY,logenergy_weighted_posX,logenergy_weighted_posY};
  string fillvar2_2d_name[] = {"h1Y","h2Y","amp_max[" + central_crystal +"]","amp_max[" + central_crystal +"]","dt_"+central_crystal+"_MCP1","dt_"+central_crystal+"_MCP1","dt_"+central_crystal+"_MCP2","dt_MCP1_MCP2","energy weighted X","energy weighted X","log(energy weighted X),log(energy weighted Y)"};

  int nvars1_2d = sizeof(fillvar1_2d)/sizeof(fillvar1_2d[0]);
  int nvars2_2d = sizeof(fillvar2_2d)/sizeof(fillvar2_2d[0]);

  double bin_low1[] = {-15.5,-15.5,-15.5,-15.5,0,0,0,-15.5,-15.5,-15.5,-15.5};
  double bin_up1[] = {15.5,15.5,15.5,15.5,628,600,600,200,15.5,15.5,15.5,15.5};
  int nbin1[] = {62,62,62,62,30,30,30,30,62,62,62,62};

  double bin_low2[] = {-15.5,-15.5,0,0,0,0,0,-45,-15.5,-15.5,-15.5,-15.5};
  double bin_up2[] = {15.5,15.5,6000,6000,6.28,6.28,6.28,-44,15.5,15.5,15.5,15.5};
  int nbin2[] = {62,62,30,30,30,30,30,30,62,62,62,62};

  if(nvars1_2d != nvars2_2d || nvars1_2d != sizeof(bin_up1)/sizeof(bin_up1[0]) || nvars1_2d != sizeof(nbin1)/sizeof(nbin1[0]) ||  nvars1_2d != sizeof(nbin2)/sizeof(nbin2[0]) ||  nvars1_2d != sizeof(bin_up2)/sizeof(bin_up2[0]))
    {cout<<"No of variable and bin information dont match for 2d histograms"<<endl;exit(0);}

  for(int ig = 10; ig <1; ig++)
    {
      for(int ie =0; ie < nenergies; ie++)
	{
	  for(int iv =0; iv < nvars1_2d; iv++)
	    {
	      TH2D* hist;
	      hist = new TH2D(("hist_nocuts_" + fillvar1_2d_name[iv] + "_" + fillvar1_2d_name[iv] + energies[ie]).c_str(),("hist_" + fillvar1_2d_name[iv] + "_"+ fillvar2_2d_name[iv]).c_str(),nbin1[iv],bin_low1[iv],bin_up1[iv],nbin2[iv],bin_low2[iv],bin_up2[iv]);
	      hist->GetXaxis()->SetTitle(fillvar1_2d_name[iv].c_str());  
	      hist->GetYaxis()->SetTitle(fillvar2_2d_name[iv].c_str());  

	      string evt_sel_cuts = "n_h1X > 0 && n_h1Y>0 && h1X > " + to_string(max_posx[ie]-3)  + " && h1X < " + to_string(max_posx[ie]+3) + " && h1Y > " + to_string(max_posy[ie]-3)  + " && h1Y < " + to_string(max_posy[ie]+3) + " &&  gain["+central_crystal+"] == " + gain_cut[ig] + " && amp_max[MCP1] > 200 && amp_max[MCP2] > 200 && " + amp_max_cut[ie];
		       
	      get2dhist(energies[ie] + "_" + central_crystal + "2x2_"+file_tag+".root",treename,fillvar2_2d[iv] + ":" + fillvar1_2d[iv],"",1,hist);

	      if(fillvar1_2d_name[iv] == "h1X" && fillvar2_2d_name[iv] == "h1Y")
		savehist2d(hist,"for_time_resmcp_energy_" + energies[ie] + "_" + central_crystal  + "_no_cuts_" + fillvar1_2d_name[iv] + "_" +  fillvar2_2d_name[iv],false,true,max_posx[ie]-3,max_posx[ie]+3,true,max_posy[ie]-3,max_posy[ie]+3);
	      else if(fillvar1_2d_name[iv] == "h1Y")
		savehist2d(hist,"for_time_resmcp_energy_" + energies[ie] + "_" + central_crystal  + "_no_cuts_" + fillvar1_2d_name[iv] + "_" +  fillvar2_2d_name[iv],false,true,max_posy[ie]-3,max_posy[ie]+3);
	      
	      else if(fillvar1_2d_name[iv] == "h1X")
		savehist2d(hist,"for_time_resmcp_energy_" + energies[ie] + "_" + central_crystal  + "_no_cuts_" + fillvar1_2d_name[iv] + "_" +  fillvar2_2d_name[iv],false,true,max_posx[ie]-3,max_posx[ie]+3);


	      else
		savehist2d(hist,"for_time_resmcp_energy_" + energies[ie] + "_" + central_crystal  + "_no_cuts_" + fillvar1_2d_name[iv] + "_" +  fillvar2_2d_name[iv],false);

	      if(iv < 8)
		{
		  TH2D* hist;
		  hist = new TH2D(("hist_nocuts_" + fillvar1_2d[iv] + "_" + fillvar1_2d[iv] + energies[ie]).c_str(),("hist_" + basicvar[iv] + "_"+ basicfitvar[iv]).c_str(),30,0,1000,30,0,6.28);
		  hist->GetXaxis()->SetTitle(basicvarname[iv].c_str());  
		  hist->GetYaxis()->SetTitle(basicfitvarname[iv].c_str());  

		  string evt_sel_cuts = "n_h1X > 0 && n_h1Y>0 && h1X > " + to_string(max_posx[ie]-3)  + " && h1X < " + to_string(max_posx[ie]+3) + " && h1Y > " + to_string(max_posy[ie]-3)  + " && h1Y < " + to_string(max_posy[ie]+3) + " &&  gain["+central_crystal+"] == " + gain_cut[ig] + " && amp_max[MCP1] > 200 && amp_max[MCP2] > 200 && " + amp_max_cut[ie];
		       
		  get2dhist(energies[ie] + "_" + central_crystal + "2x2_"+file_tag+".root",treename,basicfitvar[iv] + ":" + basicvar[iv],evt_sel_cuts,1,hist);
	  
		  savehist2d(hist,"for_time_resmcp_energy_" + energies[ie] + "_" + central_crystal  + "_no_cuts_" + basicvarname[iv] + "_" +  basicfitvarname[iv],false);

		}
	    }
	}
    }

  ///////// calculte time resolution in bins of Aeff
  
  for(int ig = 0; ig <1; ig++)
    {
      TH1D *time_res[4];
      time_res[0] = new TH1D("time_Res_0","Time resolution",nbins_aeff,Aeff_bin_edges);
      time_res[1] = new TH1D("time_Res_1","Time resolution",nbins_aeff,Aeff_bin_edges);
      time_res[2] = new TH1D("time_Res_2","Time resolution",nbins_aeff,Aeff_bin_edges);
      time_res[3] = new TH1D("time_Res_3","Time resolution",nbins_aeff,Aeff_bin_edges);

      // TGraphErrors *time_res[4];
      // time_res[0] = new TGraphErrors(nbins_aeff);
      // time_res[1] = new TGraphErrors(nbins_aeff);
      // time_res[2] = new TGraphErrors(nbins_aeff);
      // time_res[3] = new TGraphErrors(nbins_aeff);

      TH2D *h2d_time_res_vs_aeff[4];
      h2d_time_res_vs_aeff[0] = new TH2D("h2d_time_res_vsaeff_0","Time resolution",30,0,1500,30,0,6.238);
      h2d_time_res_vs_aeff[1] = new TH2D("time_Res_vsaeff_1","Time resolution",30,0,1500,30,0,6.238);
      h2d_time_res_vs_aeff[2] = new TH2D("time_Res_vsaeff_2","Time resolution",30,0,1500,30,-43.5,-45.5);
      h2d_time_res_vs_aeff[3] = new TH2D("time_Res_vsaeff_3","Time resolution",30,0,1500,30,0,6.238);

      for(int iv =0; iv < 7; iv++)  //////// Plot Aeff variable
	{
	  string str = basicvarname[iv];

	  for(int ie =0; ie < nenergies; ie++)
	    {
	
	      TH1D* hist1 = new TH1D(("hist_nocuts1_" + basicvarname[iv]).c_str(),("hist" + basicvarname[iv]).c_str(),nbin[iv],bin_low[iv],bin_up[iv]);
	      gethist(energies[ie]+"_"+central_crystal+"2x2_"+file_tag+".root",treename,basicvar[iv],"gain["+central_crystal+"] == 1 && " + basic_cut[iv],1,hist1);
	      if(hist1->Integral() >  100) plot_hist(hist1,str,energies[ie]+"_energy_time_res_gain1_" + basicvarname[iv] + ".png");
	      
	      TH1D* hist10 = new TH1D(("hist_nocuts10_" + basicvarname[iv]).c_str(),("hist" + basicvarname[iv]).c_str(),nbin[iv],bin_low[iv],bin_up[iv]);
	      gethist(energies[ie]+"_"+central_crystal+"2x2_"+file_tag+".root",treename,basicvar[iv],"gain["+central_crystal+"] == 10 && " + basic_cut[iv],1,hist10);
	      if(hist10->Integral() >  100) plot_hist(hist10,str,energies[ie]+"_energy_time_res_gain10_" + basicvarname[iv] + ".png");
	      
	      TH1D* histall = new TH1D(("hist_nocutsall_" + basicvarname[iv]).c_str(),("hist" + basicvarname[iv]).c_str(),nbin[iv],bin_low[iv],bin_up[iv]);
	      gethist(energies[ie]+"_"+central_crystal+"2x2_"+file_tag+".root",treename,basicvar[iv],basic_cut[iv],1,histall);
	      if(histall->Integral() >  100) plot_hist(histall,str,energies[ie]+"_energy_time_res_gainall_" + basicvarname[iv] + ".png");

	      if(iv==0)
		{
		  str = "amp max of C3";
		  TH1D* hist_amp = new TH1D("hist_amp_maxC3","amp_max[C3]",60,0,6000);
		  gethist(energies[ie]+"_"+central_crystal+"2x2_"+file_tag+".root",treename,"amp_max[C3]","h1X > " + to_string(max_posx[ie]-3)  + " && h1X < " + to_string(max_posx[ie]+3) + " && h1Y > " + to_string(max_posy[ie]-3)  + " && h1Y < " + to_string(max_posy[ie]+3) + " &&  gain["+central_crystal+"] == " + gain_cut[ig] + " && amp_max[MCP1] > 0 && amp_max[MCP2] > 0 && " + amp_max_cut[ie],1,hist_amp);
		   plot_hist(hist_amp,str,energies[ie]+"_energy_time_res_gain10_amp_maxC3.png");
		   //exit(0);
		}
	    }
	  
	  if(iv > 3)     /////////// time resolution
	    {
	      time_res[iv-4]->GetXaxis()->SetTitle(basicvarname[iv].c_str());
	      time_res[iv-4]->GetYaxis()->SetTitle(time_res_title[iv-4].c_str());

	      h2d_time_res_vs_aeff[iv-4]->GetXaxis()->SetTitle(basicvarname[iv].c_str());
	      h2d_time_res_vs_aeff[iv-4]->GetYaxis()->SetTitle(basicfitvarname[iv].c_str());

	      for(int ib = 0; ib < nbins_aeff; ib++)
		{
		  string evt_sel_cuts[nenergies];
		  string cut[nenergies];
		  TH1D *hist;
		  
		  for(int ie =0; ie < nenergies; ie++)
		    {
		      // evt_sel_cuts[ie] = "amp_max[MCP1] > 150 && amp_max[MCP2] > 150 && amp_max[C2] > 150";
		      
		      // evt_sel_cuts[ie] = "amp_max[MCP1] > 100 && amp_max[MCP2] > 100 && amp_max["+central_crystal+"] > 400 &&  gain["+central_crystal+"] == " + gain_cut[ig];

		      evt_sel_cuts[ie] = "n_h1X > 0 && n_h1Y>0 && h1X > " + to_string(max_posx[ie]-3)  + " && h1X < " + to_string(max_posx[ie]+3) + " && h1Y > " + to_string(max_posy[ie]-3)  + " && h1Y < " + to_string(max_posy[ie]+3) + " &&  gain["+central_crystal+"] == " + gain_cut[ig] + " && amp_max[MCP1] > 200 && amp_max[MCP2] > 200 && " + amp_max_cut[ie];
		      cut[ie] = basicvar[iv] + " > " + to_string(int(Aeff_bin_edges[ib])) + " && " + basicvar[iv] + " < " + to_string(int(Aeff_bin_edges[ib+1])) + " && " + evt_sel_cuts[ie];
		      //if(iv == 4 || iv == 5 || iv == 7) cut = cut + aeff_gain_cutC2;
		      
		      if(iv == 6)
			{
			  evt_sel_cuts[ie] = "amp_max[MCP1] > 0 && amp_max[MCP2] > 0 ";
			  hist = new TH1D(("hist_" + to_string(iv) + basicfitvarname[iv]).c_str(),"dt_MCP1_MCP2",40,-500,500);
			}
		      else
			hist = new TH1D(("hist_" + to_string(iv) + basicfitvarname[iv]).c_str(),("dt_"+central_crystal+"_MCP1").c_str(),12,0,6.238);
		      
		      gethist(energies[ie]+"_"+central_crystal+"2x2_"+file_tag+".root",treename,basicfitvarname[iv],cut[ie],1,hist);
		      get2dhist(energies[ie] + "_" + central_crystal + "2x2_"+file_tag+".root",treename,basicfitvar[iv] + ":" + basicvar[iv] ,evt_sel_cuts[ie],1,h2d_time_res_vs_aeff[iv-4]);

		    }

		  string cut_descp = "All" + basicvarname[iv] + " > " + to_string(int(Aeff_bin_edges[ib])) + " && " + basicvarname[iv] + " < " + to_string(int(Aeff_bin_edges[ib+1])) + "_gain["+central_crystal+"]_" + gain_cut[ig] ;
		  
		  int binmax = hist->GetMaximumBin();
		  double xmean = hist->GetXaxis()->GetBinCenter(binmax);
		  
		  //double std = hist->GetStdDev();	  
		  double std;
		  std = abs(xmean - hist->GetXaxis()->GetBinCenter(hist->FindFirstBinAbove(hist->GetBinContent(binmax)/2.5)));
		  //std = max(std, 0.3);
		  if(iv == 7)  std = hist->GetStdDev();

		  
		  double xmin,xmax;
		  //xmin = double(int((xmean - 4.8*max(0.01,std))*100))/100; xmax = double(int((xmean + 5*max(0.01,std))*100))/100;		  
		  xmin = 0.8*xmean; xmax = 1.2*xmean;
		  //xmin = 4.; xmax = 5;
		  if(xmin < 0.1 && xmean > 0.1) {xmin = 0.0; xmax = 0.4;}
		  if(xmean < 6.238 && xmax > 6.238) { xmax = 6.238; xmin = 2*xmean - 1.017*xmax;}
		  //xmin = 0.0; xmax = 6.28;
		  if (iv == 6) { xmin = -44.55; xmax = -44.33;}
		  
		  int nbins = 30;
		  //if(iv == 6)
		    //  nbins = 20;
		  if(iv==6 && ib == nbins_aeff-1) nbins = 20;
		  TH1D *hist_fit = new TH1D((basicfitvarname[iv] + to_string(iv)).c_str(),basicfitvarname[iv].c_str(),nbins,xmin,xmax);

		  //hist_fit->Draw();
		  
		  hist_fit->Sumw2();

		  for(int ie =0; ie < nenergies; ie++)
		    {
		      gethist(energies[ie]+"_"+central_crystal+"2x2_"+file_tag+".root",treename,basicfitvar[iv],cut[ie],1,hist_fit);
		    }
		  //cout<<"At energy "<<energies[ie]<<", fit of "<<fillvar[iv]<<" range choosen = "<<xmin<<" - "<<xmax<<" with max at "<<xmean<<" std dev = "<<std<<endl;
		  
		  /*   double *r_values = get_resolution_value_2gauss_diffmean(hist_fit,cut_descp,xmin,xmax);
		  double resolution_values = (r_values[0]*r_values[6] + r_values[4]*(1-r_values[6]) );
		  double resolution_error = (r_values[1]*r_values[6] + r_values[5]*(1-r_values[6]) );*/

		  double *r_values = get_resolution_value_gauss_const(hist_fit,cut_descp,xmin,xmax);
		  double resolution_values = r_values[0];
		  double resolution_error = r_values[1];

		  
		  if (iv == 6) { resolution_error = resolution_error/sqrt(2.0); resolution_values = resolution_values/sqrt(2.0);}
		  
		  time_res[iv-4]->SetBinContent(ib+1,resolution_values*1000.0);
		  time_res[iv-4]->SetBinError(ib+1,resolution_error*1000.0);
		}
	      
	      TH1D * hists[1] = {time_res[iv-4]};
	      string str[1] = {""};
	      plot_1dhists_withfit(1,hists,str,"all_energy_"+central_crystal+"_gain" + gain_cut[ig]+ "_timeres_vs" + basicfitvarname[iv] + ".png");

	      savehist2d(h2d_time_res_vs_aeff[iv-4],"all_energy_"+central_crystal+"_gain" + gain_cut[ig] + "_timeres_vs" + basicfitvarname[iv-4]+"_hist2d",false);

	    }
      if(iv == 6)
	{
	  TH1D * time_res_withoutmcp = (TH1D*)time_res[1]->Clone();

	  for(int ib = 0; ib < 7; ib++)
	    {
	      double dt_c2_mcp = time_res[0]->GetBinContent(ib+1);
	      double dt_mcps = time_res[2]->GetBinContent(ib+1);
	      double dt_c2 = sqrt(dt_c2_mcp*dt_c2_mcp - dt_mcps*dt_mcps);

	      double dt_c2_mcp_error = time_res[0]->GetBinError(ib+1);
	      double dt_mcps_error = time_res[2]->GetBinError(ib+1);
	      double dt_error = (pow(dt_mcps*dt_mcps_error,1) + pow(dt_c2_mcp*dt_c2_mcp_error,1)) / dt_c2  ;

	      double dt_c2_mcp2 = time_res[0]->GetBinContent(ib+1);
	      double dt_c2_2 = sqrt(dt_c2_mcp2*dt_c2_mcp2 - dt_mcps*dt_mcps);

	      double dt_c2_mcp2_error = time_res[0]->GetBinError(ib+1);
	      double dt_error_2 = (pow(dt_mcps*dt_mcps_error,1) + pow(dt_c2_mcp2*dt_c2_mcp2_error,1)) / dt_c2  ;

	      time_res_withoutmcp->SetBinContent(ib+1,dt_c2);
	      time_res_withoutmcp->SetBinError(ib+1,dt_error);
	      time_res_withoutmcp->GetYaxis()->SetTitle(("#frac{A_{eff}}{#sigma_{n}} of "+central_crystal+" and MCP1").c_str());
	      time_res_withoutmcp->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");
	    
	    }
	  TH1D * hists[1] = {time_res_withoutmcp};
	  string str[1] = {""};
	  plot_1dhists_withfit(1,hists,str,"all_energy_"+central_crystal+"_gain" + gain_cut[ig] + "_timeres_vs" + basicfitvarname[4] + "_withoutmcp.png");

	}
      	}
    }

  ///////////////// Time resolution
  string fillvar[] = {"dt_"+central_crystal+"_MCP1","dt_"+central_crystal+"_MCP2","dt_MCP1_MCP2","amp_max[MCP1]","amp_max[MCP2]","amp_max["+central_crystal+"]"};
  int nvars = sizeof(fillvar)/sizeof(fillvar[0]);
  

  TH1D * hist = new TH1D(("dt_"+central_crystal+"_MCP1").c_str(),"dt_C2_MCP1",200,0,6.28);
  hist->Sumw2();
  TH1D * hist_dt_mcps = new TH1D("dt_MCP1_MCP2_all","dt_MCP1_MCP2",100,-44.6,-44.3);
  hist_dt_mcps->Sumw2();
  TH1D * hist_res = new TH1D(("dt_"+central_crystal+"_MCP1").c_str(),"dt_C2_MCP1",60,0,6.238);

  TH1D * hist_resvs_energy[2];
  hist_resvs_energy[0] = new TH1D("time_res_energy_gain1","Time resolution vs energy",nenergies,0,nenergies);
  hist_resvs_energy[0]->GetXaxis()->SetTitle("Energy of beam (in GeV)");
  hist_resvs_energy[0]->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");

  hist_resvs_energy[1] = new TH1D("time_res_energy_gain10","Time resolution vs energy",nenergies,0,nenergies);
  hist_resvs_energy[1]->GetXaxis()->SetTitle("Energy of beam (in GeV)");
  hist_resvs_energy[1]->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");

  TH1D * hist_resvs_energy_withmcp[2];
  hist_resvs_energy_withmcp[0] = new TH1D("time_res_energy_gain1","Time resolution vs energy",nenergies,0,nenergies);
  hist_resvs_energy_withmcp[0]->GetXaxis()->SetTitle("Energy of beam (in GeV)");
  hist_resvs_energy_withmcp[0]->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");

  hist_resvs_energy_withmcp[1] = new TH1D("time_res_energy_gain10","Time resolution vs energy",nenergies,0,nenergies);
  hist_resvs_energy_withmcp[1]->GetXaxis()->SetTitle("Energy of beam (in GeV)");
  hist_resvs_energy_withmcp[1]->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");

  TH1D * hist_sigvs_mcpcut[nenergies];
  for(int ie =0; ie < nenergies; ie++)
    {
      hist_sigvs_mcpcut[ie] = new TH1D(("time_res" + energies[ie]).c_str(),"Time resolution vs energy",21,0,50*21);
      hist_sigvs_mcpcut[ie]->GetXaxis()->SetTitle("Cut on MCP amp");
      hist_sigvs_mcpcut[ie]->GetYaxis()->SetTitle("Time resolution of MCP (in ps)");
    }

  TGraphErrors * hist_sigvs_ampmaxmcp[2];
  hist_sigvs_ampmaxmcp[0] = new TGraphErrors(nenergies);
  hist_sigvs_ampmaxmcp[0]->GetXaxis()->SetTitle("A_{eff}of MCP");
  hist_sigvs_ampmaxmcp[0]->GetYaxis()->SetTitle("Time resolution of MCP (in ps)");

  hist_sigvs_ampmaxmcp[1] = new TGraphErrors(nenergies);
  hist_sigvs_ampmaxmcp[1]->GetXaxis()->SetTitle("Beam energy (in GeV) ");
  hist_sigvs_ampmaxmcp[1]->GetYaxis()->SetTitle("Time resolution of MCP (in ps)");

  TGraphErrors * hist_timeresvs_ampeff;
  hist_timeresvs_ampeff = new TGraphErrors(nenergies);
  hist_timeresvs_ampeff->GetXaxis()->SetTitle("A_{eff}of crystal");
  hist_timeresvs_ampeff->GetYaxis()->SetTitle("Time resolution of crystal (in ps)");

  vector<string> func_name = {"mod_gauss","2gauss_diffmean","2gauss_samemean","DCB","gauss","gauss_constbkg","mod_gauss_constbkg"};
  for(int ifunc = 15; ifunc < 6; ifunc++)
    {
      hist_resvs_energy[0]->Reset("ICESM");
      hist_resvs_energy[1]->Reset("ICESM");
      hist_resvs_energy_withmcp[0]->Reset("ICESM");
      hist_resvs_energy_withmcp[1]->Reset("ICESM");

      TGraphErrors *time_res[5];
      //time_res[0] = new TGraphErrors("time_Res_0","Time resolution",nbins_aeff);
      //time_res[1] = new TGraphErrors("time_Res_1","Time resolution",nbins_aeff);
      //time_res[2] = new TGraphErrors("time_Res_2","Time resolution",nbins_aeff);

      time_res[0] = new TGraphErrors(nenergies);
      time_res[1] = new TGraphErrors(nenergies);
      time_res[2] = new TGraphErrors(nenergies);
      time_res[3] = new TGraphErrors(nenergies);
      time_res[4] = new TGraphErrors(nenergies);

      time_res[3]->GetXaxis()->SetTitle((basicvarname[0+4]).c_str());
      time_res[3]->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");

      time_res[4]->GetXaxis()->SetTitle((basicvarname[1+4]).c_str());
      time_res[4]->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");

      vector<tgraphpair>mcp_reso_array;
      vector<tgraphpair>time_reso_array;

      for(int ig =0; ig < 2; ig++)
	{	  
	  for(int ie =0; ie < nenergies; ie++)
	    {
	      string evt_sel_cuts;
	      string beam_energy = energies[ie] + " GeV and gain of "+central_crystal+" = "  + gain_cut[ig];
	      
	      cout<<"Starting for energy "<<beam_energy<<endl;
	      TH1D * hist_fit[nvars];
	      double resolution_values[3] = {0,0,0};
	      double resolution_error[3] = {0,0,0};

	      double sigma_amp[3] = {0,0,0};
	      double mean_amp[3] = {0,0,0};
	      double mean_amp_error[3] = {0,0,0};
	      
	      TH1D * hists_list[3] = {NULL,NULL,NULL};
	      string str_list[3] = {" "};
	      
	      for(int iv =0; iv < 3; iv++)
		{

		  evt_sel_cuts = "n_h1X > 0 && n_h1Y>0 && h1X > " + to_string(max_posx[ie]-3)  + " && h1X < " + to_string(max_posx[ie]+3) + " && h1Y > " + to_string(max_posy[ie]-3)  + " && h1Y < " + to_string(max_posy[ie]+3) + " &&  gain["+central_crystal+"] == " + gain_cut[ig] + " && amp_max[MCP1] > 200 && amp_max[MCP2] > 200 && " + amp_max_cut[ie];
	      
		  //if(ie!=6 || iv !=5)
		  //continue;
		  time_res[iv]->GetXaxis()->SetTitle((basicvarname[iv+4]).c_str());
		  time_res[iv]->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");
		  
		 
		  hist_dt_mcps->Reset("ICESM");
		  
		  if(iv == 2)
		    hist = new TH1D(("hist_" + fillvar[iv] + energies[ie] + gain_cut[ig] ).c_str(),"dt_MCP1_MCP2",80,-500,500);
		  else if(iv > 2)
		    hist = new TH1D(("hist_" + fillvar[iv] + energies[ie] + gain_cut[ig] ).c_str(),"amp_max",40,100,2500);
		  else
		    hist = new TH1D(("hist_" + fillvar[iv] + energies[ie] + gain_cut[ig]).c_str(),("dt_"+central_crystal+"_MCP1").c_str(),50,0,6.238);
		  
		
		  double xbinmin=200+ 100*(ie),xbinmax = 1300 + 200*(ie+1);
		  int nbinx_binvar = 30;
		  if(iv == 2) {xbinmin=200; xbinmax = 600; nbinx_binvar = 30;}

		  
		  TH1D * hist_bin = new TH1D(("hist_" + to_string(iv)).c_str(),basicvarname[iv+4].c_str(),nbinx_binvar,xbinmin,xbinmax);
		  gethist(energies[ie] + "_"+central_crystal+"2x2_"+file_tag+".root",treename,basicvar[iv+4],evt_sel_cuts,1,hist_bin);

		  double* r_values;
		  r_values = get_resolution_value_doublecrystalball(hist_bin,beam_energy,xbinmin,xbinmax);
		  sigma_amp[iv] = r_values[0];
		  mean_amp[iv] = r_values[2];
		  double mean_aeff_error_L = r_values[4];
		  double mean_aeff_error_R = r_values[0];
		  mean_amp_error[iv] = 0; //(mean_aeff_error_R + mean_aeff_error_L ) / 2.0;

		  double multiplier = 1;
		  evt_sel_cuts = evt_sel_cuts +  " && " + basicvar[iv+4] + " > " + to_string(mean_amp[iv] - multiplier*mean_aeff_error_L) + " && " + basicvar[iv+4] + " < " + to_string( mean_amp[iv] + multiplier*mean_aeff_error_R);
		  gethist(energies[ie] + "_"+central_crystal+"2x2_"+file_tag+".root",treename,fillvar[iv],evt_sel_cuts,1,hist);
		  if(hist->Integral()<10)
		  continue;

		  hists_list[iv] = (TH1D*)hist->Clone();
		  
		  int binmax = hist->GetMaximumBin();
		  double xmean = hist->GetXaxis()->GetBinCenter(binmax);
		  
		  //double std = hist->GetStdDev();	  
		  double std;
		  std = abs(xmean - hist->GetXaxis()->GetBinCenter(hist->FindFirstBinAbove(hist->GetBinContent(binmax)/2.5)));
		  if(iv == 2)  std = hist->GetStdDev();
		  //double std = xmean/10;
		  
		  double xmin,xmax;
		  xmin = double(int((xmean - 4.*max(0.01,std))*100))/100; xmax = double(int((xmean + 4*max(0.01,std))*100))/100;
		  
		  if(xmin < 0 && xmean > 0) {xmin = 0.0; xmax = 2*xmean;}
		  if(xmean < 6.238 && xmax > 6.238) { xmax = 6.238; xmin = 2*xmean - 1.017*xmax;}
	      
		  //if(iv == 3){xmin = double(int((xmean - 2*max(0.01,std))*100))/100; xmax = double(int((xmean + 2*max(0.01,std))*100))/100;}
		  if (iv == 2) { xmin = -44.6; xmax = -44.3;}


		  if(iv > 2) hist_fit[iv] = new TH1D((fillvar[iv] + energies[ie] + gain_cut[ig] + "_CBgauss").c_str(),(fillvar[iv] + "_CB_modgauss").c_str(),20,xmin,xmax);
		  else hist_fit[iv] = new TH1D((fillvar[iv] + energies[ie] + gain_cut[ig] + func_name[ifunc]).c_str(),(fillvar[iv] + "_"  + func_name[ifunc]).c_str(),40,xmin,xmax);
		  //   }
		  hist_fit[iv]->Draw();
	  
		  hist_fit[iv]->Sumw2();
		  gethist(energies[ie] + "_" + central_crystal + "2x2_"+file_tag+".root",treename,fillvar[iv],evt_sel_cuts,1,hist_fit[iv]);
	  
		  cout<<endl<<" gain = "<<gain_cut[ig]<<" Integral = "<<hist_fit[iv]->Integral()<<endl;
		  if(iv == 2)
		    hist_dt_mcps->Add(hist_fit[iv]);

		  cout<<"At energy "<<energies[ie]<<", fit of "<<fillvar[iv]<<" range choosen = "<<xmin<<" - "<<xmax<<" with max at "<<xmean<<" std dev = "<<std<<endl;
	  
		  if(iv==2)
		    { //"mod_gauss","2gauss_diffmean","2gauss_samemean","DCB","gauss","gauss_constbkg","mod_gauss_constbkg"
		      if(ifunc == 0)
			{
			  r_values = get_resolution_value_modgauss(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[0] / sqrt(2.0);
			  resolution_error[iv] = r_values[1]/sqrt(2.0);
			}
		      else if(ifunc == 1)
			{
			  r_values = get_resolution_value_2gauss_diffmean(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = (r_values[0]*r_values[6] + r_values[4]*(1-r_values[6]) ) /sqrt(2.0);
			  resolution_error[iv] = (r_values[1]*r_values[6] + r_values[5]*(1-r_values[6]) ) /sqrt(2.0);
			}
		      else if(ifunc == 2)
			{
			  r_values = get_resolution_value_2gauss_samemean(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = (r_values[0]*r_values[6] + r_values[4]*(1-r_values[6]) )/sqrt(2.0);
			  resolution_error[iv] = (r_values[1]*r_values[6] + r_values[5]*(1-r_values[6]) )/sqrt(2.0);
			}
		      else if(ifunc == 3)
			{
			  r_values = get_resolution_value_doublecrystalball(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[4]/sqrt(2.0);
			  resolution_error[iv] = r_values[5]/sqrt(2.0);
			}
		      else if(ifunc == 4)
			{
			  r_values = get_resolution_value_gauss(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[0]/sqrt(2.0);
			  resolution_error[iv] = r_values[1]/sqrt(2.0);
			}
		      else if(ifunc == 5)
			{
			  r_values = get_resolution_value_gauss_const(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[0]/sqrt(2.0);
			  resolution_error[iv] = r_values[1]/sqrt(2.0);
			}
		      else if(ifunc == 6)
			{
			  r_values = get_resolution_value_modgauss_const(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[0]/sqrt(2.0);
			  resolution_error[iv] = r_values[1]/sqrt(2.0);
			}
		    }
	       		
		  else
		    {
		      //_crystalball  /// _gauss_const  doublecrystalball
		      if(ifunc == 0)
			{
			  r_values = get_resolution_value_modgauss(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[0];
			  resolution_error[iv] = r_values[1];
			}
		      else if(ifunc == 1)
			{
			  r_values = get_resolution_value_2gauss_diffmean(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[0]*r_values[6] + r_values[4]*(1-r_values[6]);
			  resolution_error[iv] = r_values[1]*r_values[6] + r_values[5]*(1-r_values[6]);
			}
		      else if(ifunc == 2)
			{
			  r_values = get_resolution_value_2gauss_samemean(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[0]*r_values[6] + r_values[4]*(1-r_values[6]);
			  resolution_error[iv] = r_values[1]*r_values[6] + r_values[5]*(1-r_values[6]);
			}
		      else if(ifunc == 3)
			{
			  r_values = get_resolution_value_doublecrystalball(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[4];
			  resolution_error[iv] = r_values[5];
			}
		      else if(ifunc == 4)
			{
			  r_values = get_resolution_value_gauss(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[0];
			  resolution_error[iv] = r_values[1];
			}
		      else if(ifunc == 5)
			{
			  r_values = get_resolution_value_gauss_const(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[0];
			  resolution_error[iv] = r_values[1];
			}
		      else if(ifunc == 6)
			{
			  r_values = get_resolution_value_modgauss_const(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[0];
			  resolution_error[iv] = r_values[1];
			}
		    }
		  //if( ( iv <3 && (resolution_values[iv] == 0 || resolution_error[iv] == 0) )|| ( iv > 2 && (sigma_amp[iv-3] == 0 || mean_amp[iv-3]  == 0)) ) { cout<<endl<<"check for varibale no. "<<iv<<" "<<sigma_amp[min(0,iv-3)]<<" "<<resolution_values[max(2,iv)]<<" "<<double(r_values[2]);exit(0);}
		  

		  if(ig==0)time_res[iv]->SetPoint(ie,mean_amp[iv],resolution_values[iv]*1000);
		  if(ig==0)time_res[iv]->SetPointError(ie,mean_amp_error[iv],resolution_error[iv]*1000);


		  cout<<"Resolution = "<<resolution_values[iv]<<" with error = "<<resolution_error[iv]<<endl;
		}
	      double resolution1 = sqrt(max(0.000000001,resolution_values[0]*resolution_values[0] - resolution_values[2]*resolution_values[2]))  *1000.0;
	      double resolution2 = sqrt(max(0.000000001,resolution_values[1]*resolution_values[1] - resolution_values[2]*resolution_values[2]))  *1000.0;

	      double reso_erro1 = (pow(resolution_values[0]*resolution_error[0],1) + pow(resolution_values[2]*resolution_error[2],1)) *1000000.0/ resolution1  ;
	      double reso_erro2 = (pow(resolution_values[1]*resolution_error[1],1) + pow(resolution_values[2]*resolution_error[2],1)) *1000000.0/ resolution2  ;

	      double resolution = (resolution1/(reso_erro1*reso_erro1) + resolution2/(reso_erro2*reso_erro2)) / (1.0/(reso_erro1*reso_erro1) + 1.0/(reso_erro2*reso_erro2));
	      double resolution_toterror = sqrt(1.0 / (1.0/(reso_erro1*reso_erro1) + 1.0/(reso_erro2*reso_erro2)));

	      double resolution_withmcp = (resolution_values[0]/(resolution_error[0]*resolution_error[0]) + resolution_values[1]/(resolution_error[1]*resolution_error[1])) / (1.0/(resolution_error[0]*resolution_error[0]) + 1.0/(resolution_error[1]*resolution_error[1])) *1000.0;
	      double resolution_toterror_withmcp = sqrt(1.0 / (1.0/(resolution_error[0]*resolution_error[0]) + 1.0/(resolution_error[1]*resolution_error[1]))) *1000.0;



	      if(ig==0)
		{
		  time_res[3]->SetPoint(ie,mean_amp[0],resolution1);
		  time_res[3]->SetPointError(ie,0,reso_erro1);

		  time_res[4]->SetPoint(ie,mean_amp[1],resolution2);
		  time_res[4]->SetPointError(ie,0,reso_erro2);
		}


	      hist_resvs_energy[ig]->SetBinContent(ie+1, resolution);
	      hist_resvs_energy[ig]->SetBinError(ie+1,resolution_toterror);
	      hist_resvs_energy[ig]->GetXaxis()->SetBinLabel(ie+1,energies[ie].c_str());

	      if(resolution_values[0] == 0)
		{
		  hist_resvs_energy[ig]->SetBinContent(ie+1, -1);
		  hist_resvs_energy[ig]->SetBinError(ie+1,0);	       
		}

	      hist_resvs_energy_withmcp[ig]->SetBinContent(ie+1, resolution_withmcp);
	      hist_resvs_energy_withmcp[ig]->SetBinError(ie+1,resolution_toterror_withmcp);
	      hist_resvs_energy_withmcp[ig]->GetXaxis()->SetBinLabel(ie+1,energies[ie].c_str());

	      if(resolution_values[0] == 0)
		{
		  hist_resvs_energy_withmcp[ig]->SetBinContent(ie+1, -1);
		  hist_resvs_energy_withmcp[ig]->SetBinError(ie+1,0);	       
		}

	   
	      //delete hist_fit[0];delete hist_fit[1];delete hist_fit[2];
	      //plot_1dhists(3,hists_list,fillvar,energies[ie] + "energy_dt_crystalMCP.png" );
	      cout<<" Crystal - MCP1 = "<<resolution_values[0]*1000<<" with error = "<<resolution_error[0]*1000<<" Crystal - MCP2 = "<<resolution_values[1]*1000<<" with error = "<<resolution_error[1]*1000<<" MCP1 - MCP2 = "<<resolution_values[2]*1000<<" with error = "<<resolution_error[2]*1000<<endl;

	      cout<<"Time resolution of crystal is "<<resolution1<<" or "<<resolution2<<" with error "<< reso_erro1<< " and "<< reso_erro2 <<" overall = "<<resolution<<" with error "<<resolution_toterror<<" ie ig = "<<" "<<ie<<" "<<ig<<endl;

	      hist_sigvs_mcpcut[ie]->SetBinContent(ig+1,resolution_values[2] *1000.0/sqrt(2.0));
	      hist_sigvs_mcpcut[ie]->SetBinError(ig+1, resolution_error[2]*1000.0/sqrt(2.0));
	      //hist_sigvs_mcpcut[ie]->GetXaxis()->SetBinLabel(ig+1,(to_string(50*ig)).c_str());

	      double amp_eff_mcp = sqrt(2.0/ ( pow((sigma_amp[0]/mean_amp[0]),2) + pow((sigma_amp[1]/mean_amp[1]),2)) );
	      double amp_eff_crystal = mean_amp[2]/sigma_amp[2];

	      //double amp_eff_mcp = sqrt(2.0/ ( pow((1.0/mean_amp[0]),2) + pow((1.0/mean_amp[1]),2)) );
	      //double amp_eff_crystal = mean_amp[2]/1.0;

	      if(ig == 0){
		tgraphpair t;
		t.x = amp_eff_mcp;
		t.y = resolution_values[2]*1000.0;
		t.xerror = 0;
		t.yerror = resolution_error[2]*1000.0;
		mcp_reso_array.push_back(t); 
	      }
	
	      if(ig== 0)
		{
		  hist_sigvs_ampmaxmcp[1]->SetPoint(ie, 25.0*(ie+1) , resolution_values[2]*1000);
		  hist_sigvs_ampmaxmcp[1]->SetPointError(ie, 25.0/2.0, resolution_error[2]*1000);
		}
	      cout<<"Sigma and mean for MCPs and crystal are "<<sigma_amp[0]<<" "<<sigma_amp[1]<<" "<<sigma_amp[2]<<" and mean = "<<mean_amp[0]<<" "<<mean_amp[1]<<" "<<mean_amp[2]<<endl;
	      cout<<"For MCP resolution : "<<resolution_values[2]<<"  , amp eff is "<<amp_eff_mcp<<endl;

	      /*if( (ig == 0 && ie < 8) || (ig == 1 && ie > 7) ){
		if(isnan(amp_eff_crystal)) { cout<<"Check varibale amp_eff_crystal for ig = "<<ig<<" ie = "<<ie<<" ("<<energies[ie]<<")"<<ifunc<<endl;exit(0);}
	     	tgraphpair t;
		t.x = amp_eff_crystal;
		t.y = resolution_withmcp;
		t.xerror = 0;
		t.yerror = resolution_toterror_withmcp;
		time_reso_array.push_back(t);
		}*/
	      //if(ifunc == 0) double* resolution_simfit = get_resolution_simfit(hist_fit[0],hist_fit[1],beam_energy);
	    }
	  TH1D * hists[2] = {hist_resvs_energy[0],hist_resvs_energy[1]};
	  string str[2] = {"gain of "+central_crystal+" = " + gain_cut[0],"gain of "+central_crystal+" = " + gain_cut[1]};
	  plot_1dhists_withfit(2,hist_resvs_energy,str, gain_cut[ig] + "energyvstimeresolution_" + func_name[ifunc] + ".png");
      
	  plot_1dhists_withfit(2,hist_resvs_energy_withmcp,str, gain_cut[ig] + "energyvstimeresolution_withmcp_" + func_name[ifunc] + ".png");

	  //double* r_values = get_resolution_value_2gauss(hist_dt_mcps,"amp_max_MCP > " + to_string(50*ig),-44.6,-44.3);
	  // double resolution_values = r_values[0]*r_values[6] + r_values[4]*(1-r_values[6]);
	  // double resolution_error = r_values[1]*r_values[6] + r_values[5]*(1-r_values[6]);

	  // hist_sigvs_mcpcut->SetBinContent(ig+1,resolution_values *1000);
	  // hist_sigvs_mcpcut->SetBinError(ig+1,resolution_error *1000);
	  // hist_sigvs_mcpcut->GetXaxis()->SetBinLabel(ig+1,(to_string(50*ig)).c_str());
      
	}
      
      
      sort(mcp_reso_array.begin(),mcp_reso_array.end(), comparestruct);	
      sort(time_reso_array.begin(),time_reso_array.end(), comparestruct);	
      
      int n =0 ;
      for(tgraphpair values : mcp_reso_array)
	{
	  hist_sigvs_ampmaxmcp[0]->SetPoint(n, values.x, values.y);
	  hist_sigvs_ampmaxmcp[0]->SetPointError(n, values.xerror, values.yerror);
	  n++;
	}
      
      n =0 ;
      for(tgraphpair values : time_reso_array)
	{
	  hist_timeresvs_ampeff->SetPoint(n, values.x, values.y);
	  hist_timeresvs_ampeff->SetPointError(n, values.xerror, values.yerror);
	  n++;
	  cout<<"Crystal resolution at "<<values.x<<" is "<<values.y<<" with "<<values.xerror<<" - "<<values.yerror<<endl;
	}

      //plot_tgraph(hist_sigvs_ampmaxmcp[0],"mcp_timeresolution_vs_ampeff_" + func_name[ifunc] + ".png");
      //plot_tgraph(hist_sigvs_ampmaxmcp[1],"mcp_timeresolution_vs_energy_" + func_name[ifunc] + ".png");
      //plot_tgraph(hist_timeresvs_ampeff,"crystal_timeresolution_vs_ampeff_" + func_name[ifunc] + ".png");
      
      TH1D * hists[1] = {hist_sigvs_mcpcut[0]};
      string str[1] = {""};
      //plot_1dhists(nenergies,hist_sigvs_mcpcut ,energies, "mcp_timeresolution_vs_amp_mcp_cut_vs_energy_" + func_name[ifunc] + ".png");

      plot_tgraph(time_res[0],"all_energy_"+central_crystal+"_gain" + gain_cut[0]+ "_timeres_vs" + basicfitvarname[4] + "_energywisebins.png");
      plot_tgraph(time_res[1],"all_energy_"+central_crystal+"_gain" + gain_cut[0]+ "_timeres_vs" + basicfitvarname[5] + "_energywisebins.png");
      plot_tgraph(time_res[2],"all_energy_"+central_crystal+"_gain" + gain_cut[0]+ "_timeres_vs" + basicfitvarname[6] + "_energywisebins.png");
      plot_tgraph(time_res[0],"all_energy_"+central_crystal+"_gain" + gain_cut[0]+ "_timeres_vs" + basicfitvarname[4] + "_energywisebins_withoutmcp.png");
      plot_tgraph(time_res[1],"all_energy_"+central_crystal+"_gain" + gain_cut[0]+ "_timeres_vs" + basicfitvarname[5] + "_energywisebins_withoutmcp.png");

    }

}

void intercystal_timeresolution()
{
  string central_crystal;

  central_crystal = "C3";

  string target_crystal = "C3DL2x2";
  string ref_crystal[3] = {"B3","B3","B2"};
  string fillvar[] = {"dt_" + central_crystal + "_" + ref_crystal[0] + "_phase_correct","dt_" + central_crystal + "_" + ref_crystal[1] + "_phase_correct","dt_" + central_crystal + "_" + ref_crystal[2] + "_phase_correct","amp_max[" + ref_crystal[0] + "]","amp_max[" + ref_crystal[1] + "]","amp_max[" + ref_crystal[2] + "]","amp_max[" + central_crystal + "]"};
  int nvars = sizeof(fillvar)/sizeof(fillvar[0]);

  TH1D * hist_timeres_vs_energy_C2_ref_crystal0[2];
  hist_timeres_vs_energy_C2_ref_crystal0[0] = new TH1D(("time_res_energy_gain1_" + central_crystal + "_" + ref_crystal[0]).c_str(),"Time resolution vs energy",nenergies,0,nenergies);
  hist_timeres_vs_energy_C2_ref_crystal0[0]->GetXaxis()->SetTitle("Energy of beam (in GeV)");
  hist_timeres_vs_energy_C2_ref_crystal0[0]->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");

  hist_timeres_vs_energy_C2_ref_crystal0[1] = new TH1D(("time_res_energy_gain10_" + central_crystal + "_" + ref_crystal[0]).c_str(),"Time resolution vs energy",nenergies,0,nenergies);
  hist_timeres_vs_energy_C2_ref_crystal0[1]->GetXaxis()->SetTitle("Energy of beam (in GeV)");
  hist_timeres_vs_energy_C2_ref_crystal0[1]->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");

  TH1D * hist_timeres_vs_energy_C2_ref_crystal1[2];
  hist_timeres_vs_energy_C2_ref_crystal1[0] = new TH1D(("time_res_energy_gain1_" + central_crystal + "_" + ref_crystal[1]).c_str(),"Time resolution vs energy",nenergies,0,nenergies);
  hist_timeres_vs_energy_C2_ref_crystal1[0]->GetXaxis()->SetTitle("Energy of beam (in GeV)");
  hist_timeres_vs_energy_C2_ref_crystal1[0]->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");

  hist_timeres_vs_energy_C2_ref_crystal1[1] = new TH1D(("time_res_energy_gain10_" + central_crystal + "_" + ref_crystal[1]).c_str(),"Time resolution vs energy",nenergies,0,nenergies);
  hist_timeres_vs_energy_C2_ref_crystal1[1]->GetXaxis()->SetTitle("Energy of beam (in GeV)");
  hist_timeres_vs_energy_C2_ref_crystal1[1]->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");

  TH1D * hist_timeres_vs_energy_C2_ref_crystal2[2];
  hist_timeres_vs_energy_C2_ref_crystal2[0] = new TH1D(("time_res_energy_gain1_" + central_crystal + "_" + ref_crystal[2]).c_str(),"Time resolution vs energy",nenergies,0,nenergies);
  hist_timeres_vs_energy_C2_ref_crystal2[0]->GetXaxis()->SetTitle("Energy of beam (in GeV)");
  hist_timeres_vs_energy_C2_ref_crystal2[0]->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");

  hist_timeres_vs_energy_C2_ref_crystal2[1] = new TH1D(("time_res_energy_gain10_" + central_crystal + "_" + ref_crystal[2]).c_str(),"Time resolution vs energy",nenergies,0,nenergies);
  hist_timeres_vs_energy_C2_ref_crystal2[1]->GetXaxis()->SetTitle("Energy of beam (in GeV)");
  hist_timeres_vs_energy_C2_ref_crystal2[1]->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");

  TGraphErrors * hist_timeresvs_ampeff_C2_C3;
  hist_timeresvs_ampeff_C2_C3 = new TGraphErrors(nenergies);
  hist_timeresvs_ampeff_C2_C3->GetXaxis()->SetTitle(("A_{eff}of " + central_crystal + " and " + ref_crystal[0]).c_str());
  hist_timeresvs_ampeff_C2_C3->GetYaxis()->SetTitle("Time resolution of crystal (in ps)");

  TGraphErrors * hist_timeresvs_ampeff_C2_C1;
  hist_timeresvs_ampeff_C2_C1 = new TGraphErrors(nenergies);
  hist_timeresvs_ampeff_C2_C1->GetXaxis()->SetTitle(("A_{eff}of " + central_crystal + " and " + ref_crystal[1]).c_str());
  hist_timeresvs_ampeff_C2_C1->GetYaxis()->SetTitle("Time resolution of crystal (in ps)");

  vector<string> func_name = {"mod_gauss","2gauss_diffmean","2gauss_samemean","DCB","gauss","gauss_constbkg","mod_gauss_constbkg"};

  int ifunc = 0;
  
  const int nbins_aeff = 6;
  //float Aeff_bin_edges[nbins_aeff+1] = {0., 100., 150., 200., 300., 400., 600.0, 10000.};
  //float Aeff_bin_edges[nbins_aeff+1] = {0., 50, 100, 150, 250, 400, 600, 900};
  //float Aeff_bin_edges[nbins_aeff+1] = {0., 80.0,  150, 300, 550 ,900,1500};
  //float Aeff_bin_edges[nbins_aeff+1] = {0.,  100, 250, 550 ,1500};
  //float Aeff_bin_edges[nbins_aeff+1] = {0., 100.0,  200, 300, 400, 550, 900, 1500};
  float Aeff_bin_edges[nbins_aeff+1] = {0., 100.0,  200, 300, 400, 550, 1000}; /// currently finialised bins
  //float Aeff_bin_edges[nbins_aeff+1] = {0., 100.0,  200, 250, 320, 400, 550, 1000};
  //float Aeff_bin_edges[nbins_aeff+1] = {0., 80.0,  150, 250, 400, 600, 1000};

  string A_eff_C2_refcrystal1 = "1.0 / sqrt( pow(b_rms[" + central_crystal  + "]/amp_max[" + central_crystal + "],2) + pow(b_rms["+ref_crystal[0]+"]/amp_max["+ref_crystal[0]+"],2))"; 
  //string A_eff_C2_MCP1 = "1.0 / sqrt( pow(2/amp_max[" + central_crystal + "],2) + pow(2.0/amp_max["+ref_crystal[0]+"],2))"; 
  string A_eff_C2_refcrystal2 = "1.0 / sqrt( pow(b_rms[" + central_crystal + "]/amp_max[" + central_crystal + "],2) + pow(b_rms["+ref_crystal[1]+"]/amp_max["+ref_crystal[1]+"],2))"; 
  string A_eff_C2_refcrystal3 = "1.0 / sqrt( pow(b_rms[" + central_crystal + "]/amp_max[" + central_crystal + "],2) + pow(b_rms["+ref_crystal[2]+"]/amp_max["+ref_crystal[2]+"],2))"; 
  
  string binvar[3] = {A_eff_C2_refcrystal1,A_eff_C2_refcrystal2,A_eff_C2_refcrystal3};
  string basicvarname[3] = {"Aeff_" + central_crystal  + "_" + ref_crystal[0],"Aeff_" + central_crystal  + "_"+ref_crystal[1],"Aeff_" + central_crystal  + "_"+ref_crystal[2]};

  string basicfitvar[3] = {"dt_" + central_crystal + "_" + ref_crystal[0] + "_phase_correct","dt_" + central_crystal + "_" + ref_crystal[1] + "_phase_correct","dt_" + central_crystal + "_" + ref_crystal[2] + "_phase_correct"};
  string basicfitvarname[3] = {"dt_" + central_crystal + "_" + ref_crystal[0],"dt_" + central_crystal + "_" + ref_crystal[1],"dt_" + central_crystal + "_" + ref_crystal[2]};

  string basic_cut[8] = {"amp_max["+central_crystal+"] > 400","amp_max[MCP1]> 100","amp_max[MCP2]>100","amp_max[C4] > 20","amp_max["+central_crystal+"] > 400 && amp_max[MCP1]> 100","amp_max["+central_crystal+"] > 400 && amp_max[MCP2]> 100","amp_max[MCP1] > 100 && amp_max[MCP2]>100","amp_max[C2] > 400 && amp_max[C3] > 20"}; //// Always check

  //// to define cut values
  
  ofstream outputFile1("plots/file_amp_diff_selection_cuts_energy_all.txt",std::ofstream::trunc);
  outputFile1 << "diff amp max cut          std of dt_C2_C3 distribution" << endl;

  ofstream outputFile2("plots/file_h1X_selection_cuts_energy_all.txt",std::ofstream::trunc);
  outputFile2 << "hodoscope h1X cut range          std of dt_C2_C3 distribution" << endl;
  
  ofstream outputFile3("plots/file_h1Y_selection_cuts_energy_all.txt",std::ofstream::trunc);
  outputFile3 << "hodoscope h1Y cut range          std of dt_C2_C3 distribution" << endl;

  // For bins of energy
  //double amp_diff_cut[4] = {0.4,0.4,0.4,0.4};
  //double h1Xcut_min_cut[4] = {2,4,4,-2};
  //double h1Xcut_max_cut[4] = {6,8,8,2};
  //double h1Ycut_min_cut[4] = {-4,0,-4,-2};
  //double h1Ycut_max_cut[4] = {-2,4,0,2};

  // For bins of Aeff  /// For C3D runs
  // double amp_diff_cut[4] = {0.25,0.3,0.3,0.3};
  // double h1Xcut_min_cut[4] = {4,1,0,-3};
  // double h1Xcut_max_cut[4] = {11,7,5,4};
  // double h1Ycut_min_cut[4] = {-6,-2,-6,-2};
  // double h1Ycut_max_cut[4] = {-1,3,-4,1};

  ///trail
  // double h1Xcut_min_cut[4] = {4,-4,-1,-5};
  // double h1Xcut_max_cut[4] = {10,2,5,1};
  // double h1Ycut_min_cut[4] = {-8,-7,-8,-6};
  // double h1Ycut_max_cut[4] = {-2,-1,-2,0};

  // for bins of Aeff  /// For C3DL runs and C2
  // double amp_diff_cut[4] = {0.25,0.3,0.3,0.35};
  // double h1Xcut_min_cut[4] = {5,-4,2,3};
  // double h1Xcut_max_cut[4] = {10,0,7,9};
  // double h1Ycut_min_cut[4] = {-6,-5,-9,-8};
  // double h1Ycut_max_cut[4] = {-2,1,-4,-2};

  // For bins of Aeff  /// For C3DL runs and B3
  double amp_diff_cut[4] = {0.25,0.3,0.3,0.35};
  double h1Xcut_min_cut[4] = {-4,-5,-2,-6};
  double h1Xcut_max_cut[4] = {1,-1,3,1};
  double h1Ycut_min_cut[4] = {-5,-5,-9,-8};
  double h1Ycut_max_cut[4] = {-1,1,-4,-2};

  // For bins of Aeff  /// For C3DL runs and B2
  // double amp_diff_cut[4] = {0.25,0.3,0.3,0.35};
  // double h1Xcut_min_cut[4] = {-4,-5,-2,-6};
  // double h1Xcut_max_cut[4] = {1,-1,3,1};
  // double h1Ycut_min_cut[4] = {-5,-5,-9,-8};
  // double h1Ycut_max_cut[4] = {-1,1,-4,-2};

  if(nenergies != sizeof(amp_diff_cut)/sizeof(amp_diff_cut[0]))
    {cout<<"amp diff cut values are not equal to no of beams."<<endl;exit(0);}

  string reference_crystal = ref_crystal[0];

  ///////// To decide cut values for hodoscope and relative amp difference for each energy
  for(int ie =10; ie < nenergies; ie++) 
    {
      outputFile1 << energies[ie] <<"GeV energy beam" << endl;
      outputFile2 << energies[ie] <<"GeV energy beam" << endl;
      outputFile3 << energies[ie] <<"GeV energy beam" << endl;

      TGraphErrors *relampmax_cuts= new TGraphErrors(20);
      relampmax_cuts->GetXaxis()->SetTitle(("Relative amp_max difference of " + reference_crystal +" and "+central_crystal+" (in fraction)").c_str());
      relampmax_cuts->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");

      TGraphErrors *hodoX_cuts= new TGraphErrors(20);
      hodoX_cuts->GetXaxis()->SetTitle("Cut on hodoscope X");
      hodoX_cuts->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");

      TGraphErrors *hodoY_cuts= new TGraphErrors(20);
      hodoY_cuts->GetXaxis()->SetTitle("Cut on hodoscope Y");
      hodoY_cuts->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");

      double xmin = -2;
      double xmax = -2;
      
      /*if(ie == 3)
	{xmin = -1.7; xmax = -1;}
      else
      {xmin = -0.95; xmax = 0.25;}*/
      
      xmin = -1.25; xmax = 0;
      
      int ip = 0;
      int iph = 0;
      for(float diff_amp = 10.05 ; diff_amp < 1.01; diff_amp+=0.05)
	{
	  TH1D* hist;
	  
	  hist = new TH1D(("hist_dt_"+reference_crystal+"_"+reference_crystal + energies[ie]).c_str(),"hist",60,xmin,xmax);
	  
	  gethist(energies[ie] + "_" + target_crystal + "_" + file_tag +".root",treename,"dt_"+central_crystal+"_" + reference_crystal + "_phase_correct","amp_max[" + reference_crystal +"] > 0 && amp_max["+central_crystal+"] > 0 && fabs(amp_max["+reference_crystal+"]-amp_max["+central_crystal+"])/max(amp_max["+reference_crystal+"],amp_max["+central_crystal+"]) < " + to_string(diff_amp),1,hist);
	  
	  string str = "dt_"+central_crystal+"_"+reference_crystal+"_phase_correct";
	  plot_hist(hist,str,"energy_" + energies[ie] + target_crystal + "_" + "_decide_amp_max_cut_value" + to_string(diff_amp) + ".png");
	  
	  double * r_vs = get_resolution_value_doublecrystalball(hist,energies[ie] + " GeV and relative amp max cut " + to_string(diff_amp),xmin,xmax);  //doublecrystalball
	  double resol = max(r_vs[0],r_vs[4])/sqrt(2.0);
	  
	  relampmax_cuts->SetPoint(ip,diff_amp,resol*1000.0);
	  relampmax_cuts->SetPointError(ip,0,0);
	  ip++;
	  
	  outputFile1 << diff_amp << setw(20) << resol*1000<< setw(20) << hist->GetEntries()<< endl;
	  
	}
      
      plot_tgraph(relampmax_cuts,"selection_cuts_amp_max_diff_energy" + energies[ie] + ".png",false);
      
      for(int istart_hodocut = -15; istart_hodocut<14; istart_hodocut++)
	{
	  for(int inc = 6; inc < 7; inc++)
	    {
	      for(int hodocut = istart_hodocut ; hodocut < 16; hodocut+=inc)
		{
		  TH1D* hist;
		  hist = new TH1D(("hist_dt_"+reference_crystal+"_"+central_crystal+"_h1"  + energies[ie]).c_str(),"hist",60,xmin,xmax);
		  gethist(energies[ie] + "_" + target_crystal + "_" + file_tag +".root",treename,"dt_"+central_crystal+"_"+reference_crystal+"_phase_correct","n_h1X>0 && h1X >" + to_string(hodocut) + " &&  h1X < " + to_string(hodocut+inc) + " && amp_max["+reference_crystal+"] > 0 && amp_max["+central_crystal+"] > 0 && fabs(amp_max["+reference_crystal+"]-amp_max["+central_crystal+"])/max(amp_max["+reference_crystal+"],amp_max["+central_crystal+"]) < " + to_string(amp_diff_cut[ie]),1,hist);
		  
		  string str = "dt_"+central_crystal+"_"+reference_crystal+"_phase_correct";
		  plot_hist(hist,str,"energy_" + energies[ie] + target_crystal + "_" + "_decide_hodoscopecut_h1X_value" + to_string(hodocut) + "_to_" + to_string(hodocut+inc) + ".png");
		  
		  double * r_vs = get_resolution_value_doublecrystalball(hist,energies[ie] + " GeV and hodocope X cut " + to_string(hodocut),xmin,xmax);  //doublecrystalball
		  double resol = max(r_vs[0],r_vs[4])/sqrt(2.0);
		  
		  hodoX_cuts->SetPoint(iph,iph ,resol*1000.0);
		  hodoX_cuts->SetPointError(iph,0,0);
		  
		  outputFile2 << "cut no " << iph << " "<<hodocut<<" to "<<hodocut+inc << setw(20) << resol*1000<<  setw(20) << hist->GetEntries()<<endl;
		  
		  TH1D* hist2;
		  hist2 = new TH1D(("hist_dt_"+reference_crystal+"_"+central_crystal+"_h2"  + energies[ie]).c_str(),"hist2",60,xmin,xmax);
		  gethist(energies[ie] + "_" + target_crystal + "_" + file_tag +".root",treename,"dt_"+central_crystal+"_"+reference_crystal+"_phase_correct","n_h1Y>0 && h1Y >" + to_string(hodocut) + " &&  h1Y < " + to_string(hodocut+inc) + "&& amp_max["+reference_crystal+"] > 0 && amp_max["+central_crystal+"] > 0 && fabs(amp_max["+reference_crystal+"]-amp_max["+central_crystal+"])/max(amp_max["+reference_crystal+"],amp_max["+central_crystal+"]) < " + to_string(amp_diff_cut[ie]),1,hist2);
	      
		  plot_hist(hist2,str,"energy_" + energies[ie] + target_crystal + "_" + "_decide_hodoscopecut_h1Y_value" + to_string(hodocut) + "_to_" + to_string(hodocut+inc) + ".png");

		  r_vs = get_resolution_value_doublecrystalball(hist,energies[ie] + " GeV and hodocope X cut " + to_string(hodocut),xmin,xmax);  //doublecrystalball
		  resol = max(r_vs[0],r_vs[4])/sqrt(2.0);
		  
		  hodoY_cuts->SetPoint(iph,iph,resol*1000.0);
		  hodoY_cuts->SetPointError(iph,0,0);

		  iph++;

		  outputFile3 << "cut no " << iph <<" "<< hodocut<<" to "<<hodocut+inc << setw(20) << resol*1000<<  setw(20) << hist2->GetEntries()<<endl;

		}
	    }
	  plot_tgraph(hodoX_cuts,"selection_hodoX_diff_energy" + energies[ie] + ".png",false);
	  plot_tgraph(hodoY_cuts,"selection_hodoY_diff_energy" + energies[ie] + ".png",false);

	}
    }
    ///

  // Close the file
  outputFile1.close();
  outputFile2.close();
  outputFile3.close();

  //int max_posx = 5;
  //int max_posy = 1;

  string fillvar1_2d[] = {"h1X","h2X","h1X","h1Y","h1X","h1Y","amp_max[" + central_crystal +"]"};
  string fillvar2_2d[] = {"h1Y","h2Y","fabs(amp_max["+reference_crystal+"]-amp_max["+central_crystal+"])/max(amp_max["+reference_crystal+"],amp_max["+central_crystal+"])","fabs(amp_max["+reference_crystal+"]-amp_max["+central_crystal+"])/max(amp_max["+reference_crystal+"],amp_max["+central_crystal+"])","dt_"+central_crystal+"_"+reference_crystal+"_phase_correct","dt_"+central_crystal+"_"+reference_crystal+"_phase_correct","dt_"+central_crystal+"_"+reference_crystal+"_phase_correct"};

  int nvars1_2d = sizeof(fillvar1_2d)/sizeof(fillvar1_2d[0]);
  int nvars2_2d = sizeof(fillvar2_2d)/sizeof(fillvar2_2d[0]);

  double bin_low1[] = {-15.5,-15.5,-15.5,-15.5,-15.5,-15.5,0};
  double bin_up1[] = {15.5,15.5,15.5,15.5,15.5,15.5,4000};
  int nbin1[] = {62,62,62/3,62/3,62,62,30};

  double bin_low2[] = {-15.5,-15.5,0,0,-3,-3,-3};
  double bin_up2[] = {15.5,15.5,0.5,0.5,3,3,3};
  int nbin2[] = {62,62,10,10,60,60,60};

    if(nvars1_2d != nvars2_2d || nvars1_2d != sizeof(bin_up1)/sizeof(bin_up1[0]) || nvars1_2d != sizeof(nbin1)/sizeof(nbin1[0]) ||  nvars1_2d != sizeof(nbin2)/sizeof(nbin2[0]) ||  nvars1_2d != sizeof(bin_up2)/sizeof(bin_up2[0]))
    {cout<<"No of variable and bin information dont match for 2d histograms"<<endl;exit(0);}

    /////// Plot 2D histogram to study with cuts
  for(int ie =0; ie < nenergies; ie++)
    {
      for(int iv =0; iv < nvars1_2d; iv++)
	{
	  TH2D* hist;
	  hist = new TH2D(("hist_nocuts_" + fillvar1_2d[iv] + "_" + fillvar1_2d[iv] + energies[ie]).c_str(),("hist_" + fillvar1_2d[iv] + "_"+ fillvar2_2d[iv]).c_str(),nbin1[iv],bin_low1[iv],bin_up1[iv],nbin2[iv],bin_low2[iv],bin_up2[iv]);
	  hist->GetXaxis()->SetTitle(fillvar1_2d[iv].c_str());  
	  hist->GetYaxis()->SetTitle(fillvar2_2d[iv].c_str());  

	  //get2dhist(energies[ie] + "_" + target_crystal + "_" + file_tag +".root",treename,fillvar2_2d[iv] + ":" + fillvar1_2d[iv],"n_h1X > 0 && n_h1Y > 0 && fabs(amp_max[C2]-amp_max["+central_crystal+"])/max(amp_max[C2],amp_max["+central_crystal+"]) < 0.1 && h1X > " + to_string(max_posx[ie]-2)  + " && h1X < " + to_string(max_posx[ie]+2) + " && h1Y > " + to_string(max_posy[ie]-2)  + " && h1Y < " + to_string(max_posy[ie]+2),1,hist);
	  get2dhist(energies[ie] + "_" + target_crystal + "_" + file_tag +".root",treename,fillvar2_2d[iv] + ":" + fillvar1_2d[iv],"n_h1X > 0 && n_h1Y > 0 && amp_max[C2]>0 && amp_max[C3] > 0 && fabs(amp_max["+reference_crystal+"]-amp_max["+central_crystal+"])/max(amp_max["+reference_crystal+"],amp_max["+central_crystal+"]) <" + to_string(amp_diff_cut[ie]),1,hist);
	  
	  if(iv == 2 || iv == 3) savehist2d(hist,"energy_" + energies[ie] + target_crystal + "_" + file_tag + "_no_cuts_" + fillvar1_2d[iv] + "_relamp",true);
	  else savehist2d(hist,"energy_" + energies[ie] + target_crystal + "_" + file_tag + "_no_cuts_" + fillvar1_2d[iv] + "_" +  fillvar2_2d[iv],false);
	}
    }



  //////////// Study 1-D histograms

  //double bin_low[] = {0,0,0,0,0,-0.5,-0.5,-0.5,0,0};
  //double bin_up[] = {1000,6000,1000,2000,2000,10.5,10.5,10.5,1,1};
  //int nbin[] = {100,100,100,100,100,11,11,11,30,30};
  //
  //if(nvars != sizeof(bin_low)/sizeof(bin_low[0]) || nvars != sizeof
  //  (bin_up)/sizeof(bin_up[0]) || nvars != sizeof(nbin)/sizeof(nbin[0]))
  // {cout<<"No of variable and bin information dont match"<<endl;exit(0);}
  //
  //for(int ie =0; ie < nenergies; ie++)
  // {
  //  for(int iv =0; iv < nvars; iv++)
  //	{
  //	  TH1D* hist;
  //	  hist = new TH1D(("hist_nocuts_" + fillvar[iv] + energies[ie]).c_str(),("hist" + fillvar[iv]).c_str(),nbin[iv],bin_low[iv],bin_up[iv]);
  //	  gethist(energies[ie] + "_" + target_crystal + "_"+file_tag+".root",treename,fillvar[iv],"",1,hist);
  //	  
  //	  string str = fillvar[iv];
  //	  plot_hist(hist,str,"energy_" + energies[ie] + target_crystal + "_" + "_no_cuts_" + fillvar[iv] + ".png");
  //	      
  //	}
  //}



  //////////// End of Study 1-D histograms

  
  for(int ig = 0; ig <1; ig++)
    {

      TGraphErrors *time_res[3];
      //time_res[0] = new TGraphErrors("time_Res_0","Time resolution",nbins_aeff);
      //time_res[1] = new TGraphErrors("time_Res_1","Time resolution",nbins_aeff);
      //time_res[2] = new TGraphErrors("time_Res_2","Time resolution",nbins_aeff);

      time_res[0] = new TGraphErrors(nbins_aeff);
      time_res[1] = new TGraphErrors(nbins_aeff);
      time_res[2] = new TGraphErrors(nbins_aeff);

      TGraphErrors *aeff_vs_beamenergy[5];
      aeff_vs_beamenergy[0] = new TGraphErrors(nenergies);
      aeff_vs_beamenergy[1] = new TGraphErrors(nenergies);
      aeff_vs_beamenergy[2] = new TGraphErrors(nenergies);
      aeff_vs_beamenergy[3] = new TGraphErrors(nenergies);
      aeff_vs_beamenergy[4] = new TGraphErrors(nenergies);
      
      for(int iv = 0; iv < 1;iv++)
	{

	  string str2 = "A_{eff} / #sigma_{n} of " +  central_crystal  + " & " + ref_crystal[iv];
	  string baseline_cut;
	  TH1D * histallengy = NULL;

	  time_res[iv]->GetXaxis()->SetTitle(("A_{eff} / #sigma_{n} of " +  central_crystal  + " & " + ref_crystal[iv]).c_str());
	  time_res[iv]->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");
  	      
	  for(int ib = 0; ib < nbins_aeff; ib++)
	    {
	      string evt_sel_cuts[nenergies];
	      string cut[nenergies];
	      TH1D *hist;

	      double bincenter=0;
	      double bincenter_error=0;
	      TH1D* hist_bin;

	      hist = new TH1D(("hist_" + to_string(iv) + basicfitvarname[iv]).c_str(),("dt_"+central_crystal+"_MCP1").c_str(),40,-3,3);
	      double median = 0;
	      double median_error = 0;
	      vector<double> values;
	      
	      hist_bin = new TH1D(("hist_" + to_string(iv) + binvar[iv]).c_str(),binvar[iv].c_str(),60,Aeff_bin_edges[max(0,ib-1)],Aeff_bin_edges[min(ib+2,nbins_aeff)]);
	      hist_bin->GetXaxis()->SetTitle(("#frac{A_{eff}}{#sigma_{n}} of "+reference_crystal+" and "+central_crystal+"").c_str());

	      for(int ie =0; ie < nenergies; ie++)
		{
		  baseline_cut = "amp_max["+reference_crystal+"]>100 && amp_max["+central_crystal+"]> 100 && n_h1X > 0 && n_h1Y > 0 && fabs(amp_max["+reference_crystal+"]-amp_max["+central_crystal+"])/max(amp_max["+reference_crystal+"],amp_max["+central_crystal+"]) <" + to_string(amp_diff_cut[ie]) +  "&& h1X > " + to_string(h1Xcut_min_cut[ie])  + " && h1X < " + to_string(h1Xcut_max_cut[ie]) + " && h1Y > " + to_string(h1Ycut_min_cut[ie])  + " && h1Y < " + to_string(h1Xcut_max_cut[ie]) + " && " + binvar[iv] + " > 80";

		  if(ib == 0)
		    {
		      //baseline_cut = "";
		      TH1D* histall = new TH1D(("hist_all_" + basicvarname[iv]).c_str(),("hist" + basicvarname[iv]).c_str(),50,0,1500);
		      gethist(energies[ie]+"_"+target_crystal+"_"+file_tag+".root",treename,binvar[iv],baseline_cut,1,histall);
		      if(histall->Integral() >  100)
			{
			  aeff_vs_beamenergy[iv]->SetPoint(ie,stof(energies[ie]),histall->GetMean());
			  aeff_vs_beamenergy[iv]->GetYaxis()->SetTitle(basicvarname[iv].c_str());
			  aeff_vs_beamenergy[iv]->GetXaxis()->SetTitle("Beam energy (in GeV)");
			  
			  if(histallengy == NULL)
			    histallengy = (TH1D*)histall->Clone(("hist_allernergy_" + basicvarname[iv]).c_str());
			  else
			    histallengy->Add(histall);
			  plot_hist(histall,str2,energies[ie]+"_energy_time_res_gainall_" + basicvarname[iv] + ".png");
			}
		      string amp_max_var[2] = {"amp_max[C3]","amp_max[C2]"};
		      //string amp_max_var[2] = {"amp_max[C3]/pow(b_rms[C3]","amp_max[C2]/pow(b_rms[C2]"};
		      for(int i = 0; i < 2 ; i++)
			{
			  string str3 = amp_max_var[i];
			  TH1D* histall2 = new TH1D(("hist_all2_" + amp_max_var[i]).c_str(),("hist" + amp_max_var[i]).c_str(),50,0,6000);
			  gethist(energies[ie]+"_"+target_crystal+"_"+file_tag+".root",treename,amp_max_var[i],baseline_cut,1,histall2);
			  aeff_vs_beamenergy[i+3]->SetPoint(ie,stof(energies[ie]),histall2->GetMean());
			  aeff_vs_beamenergy[i+3]->GetYaxis()->SetTitle(amp_max_var[i].c_str());
			  aeff_vs_beamenergy[i+3]->GetXaxis()->SetTitle("Beam energy (in GeV)");
			  
			  plot_hist(histall2,str3,energies[ie]+"_energy_time_res_gainall_" + amp_max_var[i] + ".png");
			    
			}
		    }

		  //string evt_sel_cuts = "n_h1X > 0 && n_h1Y > 0 && h1X > " + to_string(max_posx[ie]-2)  + " && h1X < " + to_string(max_posx[ie]+2) + " && h1Y > " + to_string(max_posy[ie]-2)  + " && h1Y < " + to_string(max_posy[ie]+e) ;
	      
		  //string evt_sel_cuts = "amp_max[" + ref_crystal[iv] + "] > 50 && amp_max["+central_crystal+"] > 50 &&  gain["+central_crystal+"] == " + gain_cut[ig];
		  //evt_sel_cuts[ie] = "fabs(amp_max[" + ref_crystal[iv] + "]-amp_max[" + central_crystal + "])/max(amp_max[" + ref_crystal[iv] + "],amp_max[" + central_crystal + "]) < " + to_string(amp_diff_cut[ie]) + " && amp_max[" + ref_crystal[iv] + "] > 0 && amp_max["+central_crystal+"] > 0 && gain["+central_crystal+"] == " + gain_cut[ig];
		  
		  
		  evt_sel_cuts[ie] = baseline_cut;
		  
		  cut[ie] = binvar[iv] + " > " + to_string(int(Aeff_bin_edges[ib])) + " && " + binvar[iv] + " < " + to_string(int(Aeff_bin_edges[ib+1])) + " && " + evt_sel_cuts[ie];
			      
		  gethist(energies[ie]+"_"+ target_crystal + "_" + file_tag +".root",treename,basicfitvar[iv],cut[ie],1,hist);
		  gethist(energies[ie]+"_"+ target_crystal + "_" + file_tag +".root",treename,binvar[iv],cut[ie],1,hist_bin);
		  double *median_witherror  = CalculateMedian(energies[ie]+"_"+ target_crystal + "_" + file_tag +".root",treename,binvar[iv],cut[ie],values);
		  median = median_witherror[0];
		  median_error = median_witherror[1];
		  cout<<endl<<"No of entries in vector = "<<values.size()<<" and now medain is "<<median<<" after adding file "<<energies[ie]<<" GeV"<<endl;
		}
	      
	      if(ib == 0 && histallengy != NULL && histallengy->Integral() >  100)
		plot_hist(histallengy,str2,"all_energy_time_res_gainall_" + basicvarname[iv] + ".png");
	      
	      TH1D * hists_give[1] = {hist};
	      string str_give[1] = {""};
	      plot_1dhists(1,hists_give,str_give,"all_energy_"+central_crystal+"_gain" + gain_cut[ig]+ "_" + basicfitvarname[iv] + "_inbinof_" + basicvarname[iv] + "_" + to_string(int(Aeff_bin_edges[ib])) + "_" + to_string(int(Aeff_bin_edges[ib+1])) + ".png");
	      
	      hists_give[0] = hist_bin;
	      plot_1dhists(1,hists_give,str_give,"all_energy_"+central_crystal+"_gain" + gain_cut[ig]+ "_" + basicvarname[iv] + "_inbinof_" + basicvarname[iv] + "_" + to_string(int(Aeff_bin_edges[ib])) + "_" + to_string(int(Aeff_bin_edges[ib+1])) + ".png");

	      //bincenter = hist_bin->GetMean();
	      bincenter = median;
	      bincenter_error = median_error;

	      string cut_descp = "All" + basicvarname[iv] + " > " + to_string(int(Aeff_bin_edges[ib])) + " && " + basicvarname[iv] + " < " + to_string(int(Aeff_bin_edges[ib+1]));
		  
	      int binmax = hist->GetMaximumBin();
	      double xmean = hist->GetXaxis()->GetBinCenter(binmax);
	      
	      //double std = hist->GetStdDev();	  
	      double std;
	      std = abs(xmean - hist->GetXaxis()->GetBinCenter(hist->FindFirstBinAbove(hist->GetBinContent(binmax)/2.5)));
	      std = hist->GetStdDev();
	      //double std = xmean/10;
	      
	      double xmin,xmax;
	      xmin = double(int((xmean - 3.5*max(0.0001,std))*100))/100; xmax = double(int((xmean + 3.5*max(0.01,std))*100))/100;
	      
	      /*    //if(iv == 0 &&  ib == 0  )	      {xmin = -1.2; xmax = -0.1;}
	      if(iv == 0 &&  ib == 2  )	      {xmin = -1.3; xmax = -0.05;}
	      //if(iv == 0 &&  ib > 2  )	      {xmin = -1.2; xmax = -0.2;}
	      if(iv == 0 &&  ib == 4  )	      {xmin = -1; xmax = 0.4;}
	      //if(iv == 0 &&  ib == 5  )	      {xmin = -1; xmax = 0.4;}
	      if(iv == 0 &&  ib == 5  )	      {xmin = -3.2; xmax = -0.4;}
	      */
	      /*if(iv == 0 &&  ib == 0  )	      {xmin = -1.15; xmax = 0.45;}
	      if(iv == 0 &&  ib == 1  )	      {xmin = -0.9; xmax = 0.2;}
	      if(iv == 0 &&  ib >1  )	      {xmin = -0.7; xmax = 0.0;}
	      if(iv == 0 &&  ib == 4  )	      {xmin = -0.75; xmax = -0.05;}
	      if(iv == 0 &&  ib == 5  )	      {xmin = -0.6; xmax = -0.1;}*/
	      //if(iv == 0 && ib < 2 && ib >0 )	      {	      xmin = double(int((xmean - 3*max(0.01,std))*100))/100; xmax = double(int((xmean + 3*max(0.01,std))*100))/100;}

	      if(iv == 0 &&  ib == 0  )	      {xmin = -1.5; xmax = 1;}
	      if(iv == 0 && ib ==nbins_aeff-3 )	      {	      xmin = double(int((xmean - 0.8*max(0.01,std))*100))/100; xmax = double(int((xmean + 0.8*max(0.01,std))*100))/100;}
	      if(iv == 0 && ib ==2 )	      {	      xmin = -0.8; xmax = 0.4;}
	      if(iv == 0 && ib ==nbins_aeff-3 )	      {	      xmin = double(int((xmean - 0.7*max(0.01,std))*100))/100; xmax = double(int((xmean + 0.7*max(0.01,std))*100))/100;}
	      if(iv == 0 && ib >nbins_aeff-3 )	      {	      xmin = double(int((xmean - 0.6*max(0.01,std))*100))/100; xmax = double(int((xmean + 0.6*max(0.01,std))*100))/100;}
	      // if(iv == 0 && ib == nbins_aeff-2 )	      {	      xmin = double(int((xmean - 0.8*max(0.01,std))*100))/100; xmax = double(int((xmean + 0.8*max(0.01,std))*100))/100;}
	      // if(iv == 0 && ib == nbins_aeff-1 )	      {	      xmin = double(int((xmean - 0.1*max(0.01,std))*100))/100; xmax = double(int((xmean + 0.1*max(0.01,std))*100))/100;}
	      if(iv == 0 && ib ==5 )	      {	      xmin = -0.75; xmax = 0.15;}
	      if(iv == 0 && ib ==6 )	      {	      xmin = double(int((xmean - 0.55*max(0.01,std))*100))/100; xmax = double(int((xmean + 0.55*max(0.01,std))*100))/100;}
if(iv == 0 && ib == 0 )	      xmin = -0.8;
	      //xmin = -3 ; xmax = 3;
	      if(xmin < -3) {xmin = -3;}
	      if(xmax > 3 ) { xmax = 3;}
	      
	      int nbins = 30;
	      //if(iv == 0 && ib >3 ) nbins = 20;
	      //if(iv == 0 && ib > 2 ) nbins = 20;
	      if(iv == 0 && ib == nbins_aeff-1 ) nbins = 20;
	      TH1D *hist_fit = new TH1D((basicfitvarname[iv] + to_string(iv)).c_str(),basicfitvarname[iv].c_str(),nbins,xmin,xmax);
	      
	      //hist_fit->Draw();
	      
	      hist_fit->Sumw2();
	      for(int ie =0; ie < nenergies; ie++)
		gethist(energies[ie]+"_"+target_crystal+"_"+file_tag+".root",treename,basicfitvar[iv],cut[ie],1,hist_fit);
	      
	      //cout<<"At energy "<<energies[ie]<<", fit of "<<fillvar[iv]<<" range choosen = "<<xmin<<" - "<<xmax<<" with max at "<<xmean<<" std dev = "<<std<<endl;
	      
	      //double *r_values = get_resolution_value_2gauss_diffmean(hist_fit,cut_descp,xmin,xmax);
	      //double resolution_values = (r_values[0]*r_values[6] + r_values[4]*(1-r_values[6]) ) /sqrt(2.0);
	      //double resolution_error = (r_values[1]*r_values[6] + r_values[5]*(1-r_values[6]) ) /sqrt(2.0);
	      
	      double *r_values = get_resolution_value_gauss_const(hist_fit,cut_descp,xmin,xmax);
	      double resolution_values = r_values[0] /sqrt(2.0);
	      double resolution_error = r_values[1]/sqrt(2.0);
	      
	      //double* r_values = get_resolution_value_doublecrystalball(hist_fit,cut_descp,xmin,xmax);
	      //double resolution_values = r_values[4]/sqrt(2.0);
	      //double resolution_error = r_values[5]/sqrt(2.0);
	      
	      if(ib ==100)
	      continue;
	      time_res[iv]->SetPoint(ib,bincenter,resolution_values*1000.0);
	      time_res[iv]->SetPointError(ib,bincenter_error,resolution_error*1000.0);
	    }
	  
	  //TH1D * hists[1] = {time_res[iv]};
	  //string str[1] = {""};
	  //plot_1dhists_withfit(1,hists,str,"all_energy_"+central_crystal+"_gain" + gain_cut[ig]+ "_timeres_vs" + basicfitvarname[iv] + ".png");
	  plot_tgraph(time_res[iv],"all_energy_"+central_crystal+"_gain" + gain_cut[ig]+ "_timeres_vs" + basicfitvarname[iv] + ".png");
	}
    
      plot_tgraph(aeff_vs_beamenergy[0],"all_energy_"+central_crystal+"_gain" + gain_cut[ig]+ "_beam_energy_vs_" + basicvarname[0] + ".png",true,"linear");
      plot_tgraph(aeff_vs_beamenergy[3],"all_energy_"+central_crystal+"_gain" + gain_cut[ig]+ "_beam_energy_vs_amp_maxc3_sigma.png",true,"linear");
      plot_tgraph(aeff_vs_beamenergy[4],"all_energy_"+central_crystal+"_gain" + gain_cut[ig]+ "_beam_energy_vs_amp_maxC2_sigma.png",true,"linear");
    }


  ////////////   Time resolution in bins of energy
  
  for(int ifunc = 15; ifunc < 6; ifunc++)
    {
      hist_timeres_vs_energy_C2_ref_crystal0[0]->Reset("ICESM");
      hist_timeres_vs_energy_C2_ref_crystal0[1]->Reset("ICESM");

      hist_timeres_vs_energy_C2_ref_crystal1[0]->Reset("ICESM");
      hist_timeres_vs_energy_C2_ref_crystal1[1]->Reset("ICESM");

      hist_timeres_vs_energy_C2_ref_crystal2[0]->Reset("ICESM");
      hist_timeres_vs_energy_C2_ref_crystal2[1]->Reset("ICESM");

      TGraphErrors *time_res[3];
      //time_res[0] = new TGraphErrors("time_Res_0","Time resolution",nbins_aeff);
      //time_res[1] = new TGraphErrors("time_Res_1","Time resolution",nbins_aeff);
      //time_res[2] = new TGraphErrors("time_Res_2","Time resolution",nbins_aeff);

      time_res[0] = new TGraphErrors(nenergies);
      time_res[1] = new TGraphErrors(nenergies);
      time_res[2] = new TGraphErrors(nenergies);


      vector<tgraphpair>time_reso_C2_C3_array;
      vector<tgraphpair>time_reso_C2_C1_array;
      
      for(int ig =0; ig < 1; ig++)
	{
	  for(int iv =0; iv < 1; iv++) //nvars
	    {
	      string evt_sel_cuts;
	      
	      TH1D * hist;
	      TH1D * hist_fit[nvars];
	      
	      double resolution_values[3] = {0,0,0};
	      double resolution_error[3] = {0,0,0};
	      
	      double sigma_amp[4] = {0,0,0};
	      double mean_amp[4] = {0,0,0};

	      double mean_aeff;
	      double mean_aeff_error;
	      double mean_aeff_error_L;
	      double mean_aeff_error_R;
	      TH1D* hist_bin;

	      for(int ie =0; ie < nenergies; ie++)
		{
		  
		  string beam_energy = energies[ie] + " GeV and gain of " + central_crystal + " = "  + gain_cut[ig];
		  evt_sel_cuts = amp_max_cut[ie] + " && gain[" + central_crystal + "] == " + gain_cut[ig];
	      
		  time_res[iv]->GetXaxis()->SetTitle(("A_{eff} / #sigma_{n} of " +  central_crystal  + " & " + ref_crystal[iv]).c_str());
		  time_res[iv]->GetYaxis()->SetTitle("Time resolution of Crystal (in ps)");

	  
		  evt_sel_cuts = "amp_max["+reference_crystal+"]>100 && amp_max["+central_crystal+"]> 100 && n_h1X > 0 && n_h1Y > 0 && fabs(amp_max["+reference_crystal+"]-amp_max["+central_crystal+"])/max(amp_max["+reference_crystal+"],amp_max["+central_crystal+"]) <" + to_string(amp_diff_cut[ie]) +  "&& h1X > " + to_string(h1Xcut_min_cut[ie])  + " && h1X < " + to_string(h1Xcut_max_cut[ie]) + " && h1Y > " + to_string(h1Ycut_min_cut[ie])  + " && h1Y < " + to_string(h1Xcut_max_cut[ie]) + " && " + binvar[iv] + " > 80";
		  
		  //cut[ie] = binvar[iv] + " > " + to_string(int(Aeff_bin_edges[ib])) + " && " + binvar[iv] + " < " + to_string(int(Aeff_bin_edges[ib+1])) + " && " + evt_sel_cuts[ie];

		  double xbinmin = 30*(ie+1),xbinmax = 500*(ie+3)/2;
		  int bins_bin = 35;
		  if(ie > 1) bins_bin =30; 
		  if(ie > 2) bins_bin =20; 
		  hist_bin = new TH1D(("hist_" + to_string(iv) + binvar[iv]).c_str(),basicvarname[iv].c_str(),bins_bin,xbinmin,xbinmax);
		    
		  gethist(energies[ie]+"_"+ target_crystal + "_" + file_tag +".root",treename,binvar[iv],evt_sel_cuts/*cut[ie]*/,1,hist_bin);
		  
		  double *mean_aeff_values  = get_resolution_value_doublecrystalball(hist_bin,beam_energy,xbinmin,xbinmax);
		  mean_aeff = mean_aeff_values[2];
		  mean_aeff_error_L = 1*mean_aeff_values[4];
		  mean_aeff_error_R = 1*mean_aeff_values[0];
		  mean_aeff_error = (mean_aeff_error_R + mean_aeff_error_L ) / 2.0;

		  
		  //evt_sel_cuts += " && " + binvar[iv] + " > " + to_string(int(Aeff_bin_edges[ib])) + " && " + binvar[iv] + " < " + to_string(int(Aeff_bin_edges[ib+1])); 
		  hist = new TH1D(("hist_" + fillvar[iv] + energies[ie] + gain_cut[ig] ).c_str(),"dt",100,-3,3);

		  double multiplier = 1;
		  evt_sel_cuts = evt_sel_cuts + " && " + binvar[iv] + " > " + to_string(mean_aeff - multiplier*mean_aeff_error_L) + " && " + binvar[iv] + " < " + to_string(mean_aeff + multiplier*mean_aeff_error_R);
		  gethist(energies[ie]+"_"+ target_crystal + "_" + file_tag +".root",treename,basicfitvar[iv],evt_sel_cuts,1,hist);
		  
		  //if(hist->Integral()<1000)
		  //continue;
	      	  
		  int binmax = hist->GetMaximumBin();
		  double xmean = hist->GetXaxis()->GetBinCenter(binmax);
		  
		  //double std = hist->GetStdDev();	  
		  double std;
		  std = abs(xmean - hist->GetXaxis()->GetBinCenter(hist->FindFirstBinAbove(hist->GetBinContent(binmax)/2.5)));
		  if(iv == 2)  std = hist->GetStdDev();
		  //double std = xmean/10;
		  
		  double xmin,xmax;
		  xmin = double(int((xmean - 5.*max(0.01,std))*100))/100; xmax = double(int((xmean + 5*max(0.01,std))*100))/100;

		  if(iv == 0 &&  ie == 1  )	      {xmin = -0.9; xmax = 0;}
		  //if(iv == 0 &&  ib > 2  )	      {xmin = -1.2; xmax = -0.2;}
		  //if(iv == 0 &&  ie == 3  )	      {xmin = -0.4; xmax = -0.;}
		  if(iv == 0 &&  ie == 3  )	      {xmin = -1.9; xmax = -0.55;}
		  //if(iv == 0 &&  ib == 5  )	      {xmin = -1; xmax = 0.4;}
		  if(iv == 0 &&  ie == 5  )	      {xmin = -3.2; xmax = -0.4;}

	      
		  if(xmin < -3 && iv < 2) {xmin = -3;}
		  if(xmax > 3 && iv < 2) { xmax = 3;}
		  
		  if(iv > 2) hist_fit[iv] = new TH1D((fillvar[iv] + energies[ie] + gain_cut[ig] + "_CBgauss").c_str(),(fillvar[iv] + "_CB_modgauss").c_str(),50,xmin,xmax);
		  else hist_fit[iv] = new TH1D((fillvar[iv] + energies[ie] + gain_cut[ig] + func_name[ifunc]).c_str(),(fillvar[iv] + "_"  + func_name[ifunc]).c_str(),30,xmin,xmax);
		  
		  hist_fit[iv]->Sumw2();
		  
		  gethist(energies[ie] + "_" + target_crystal + "_"+file_tag+".root",treename,fillvar[iv],evt_sel_cuts,1,hist_fit[iv]);

		  
		  cout<<"At energy "<<energies[ie]<<", fit of "<<fillvar[iv]<<" range choosen = "<<xmin<<" - "<<xmax<<" with max at "<<xmean<<" std dev = "<<std<<endl;
		  
		  double* r_values;
		  if(iv < 3)
		    { //"mod_gauss","2gauss_diffmean","2gauss_samemean","DCB","gauss","gauss_constbkg","mod_gauss_constbkg"
		      if(ifunc == 0)
			{
			  r_values = get_resolution_value_modgauss(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[0] / sqrt(2.0);
			  resolution_error[iv] = r_values[1]/sqrt(2.0);
			}
		      else if(ifunc == 1)
			{
			  r_values = get_resolution_value_2gauss_diffmean(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = (r_values[0]*r_values[6] + r_values[4]*(1-r_values[6]) ) /sqrt(2.0);
			  resolution_error[iv] = (r_values[1]*r_values[6] + r_values[5]*(1-r_values[6]) ) /sqrt(2.0);
			}
		      else if(ifunc == 2)
			{
			  r_values = get_resolution_value_2gauss_samemean(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = (r_values[0]*r_values[6] + r_values[4]*(1-r_values[6]) )/sqrt(2.0);
			  resolution_error[iv] = (r_values[1]*r_values[6] + r_values[5]*(1-r_values[6]) )/sqrt(2.0);
			}
		      else if(ifunc == 3)
			{
			  r_values = get_resolution_value_doublecrystalball(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[4]/sqrt(2.0);
			  resolution_error[iv] = r_values[5]/sqrt(2.0);
			}
		      else if(ifunc == 4)
			{
			  r_values = get_resolution_value_gauss(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[0]/sqrt(2.0);
			  resolution_error[iv] = r_values[1]/sqrt(2.0);
			}
		      else if(ifunc == 5)
			{
			  r_values = get_resolution_value_gauss_const(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[0]/sqrt(2.0);
			  resolution_error[iv] = r_values[1]/sqrt(2.0);
			}
		      else if(ifunc == 6)
			{
			  r_values = get_resolution_value_modgauss_const(hist_fit[iv],beam_energy,xmin,xmax);
			  resolution_values[iv] = r_values[0]/sqrt(2.0);
			  resolution_error[iv] = r_values[1]/sqrt(2.0);
			}
		    }
	      

		  time_res[iv]->SetPoint(ie,mean_aeff,resolution_values[0]*1000.0);
		  time_res[iv]->SetPointError(ie,mean_aeff_error,resolution_error[1]*1000.0);
	      
		  hist_timeresvs_ampeff_C2_C3->SetPoint(ie,mean_aeff,resolution_values[0]*1000.0);
		  hist_timeresvs_ampeff_C2_C3->SetPointError(ie,mean_aeff_error,resolution_error[1]*1000.0);
		}
	      
	      //hist_timeres_vs_energy_C2_ref_crystal0[ig]->SetBinContent(ie+1, resolution_values[0]*1000);
	      //hist_timeres_vs_energy_C2_ref_crystal0[ig]->SetBinError(ie+1, resolution_error[0]*1000);
	      //hist_timeres_vs_energy_C2_ref_crystal0[ig]->GetXaxis()->SetBinLabel(ie+1,energies[ie].c_str());

	      //hist_timeres_vs_energy_C2_ref_crystal1[ig]->SetBinContent(ie+1, resolution_values[1]*1000);
	      //hist_timeres_vs_energy_C2_ref_crystal1[ig]->SetBinError(ie+1,resolution_error[1]*1000);
	      //hist_timeres_vs_energy_C2_ref_crystal1[ig]->GetXaxis()->SetBinLabel(ie+1,energies[ie].c_str());

	      //hist_timeres_vs_energy_C2_ref_crystal2[ig]->SetBinContent(ie+1, resolution_values[2]*1000);
	      //hist_timeres_vs_energy_C2_ref_crystal2[ig]->SetBinError(ie+1,resolution_error[2]*1000);
	      //hist_timeres_vs_energy_C2_ref_crystal2[ig]->GetXaxis()->SetBinLabel(ie+1,energies[ie].c_str());
	   
	      cout<<"Time resolution of crystal is "<<resolution_values[0]*1000<<" or "<<resolution_values[1]*1000<<" with error "<<resolution_error[0]*1000<< " and "<<resolution_error[1]*1000<<endl;


	      plot_tgraph(time_res[iv],"all_energy_"+central_crystal+"_gain" + gain_cut[ig]+ "_timeres_vs" + basicfitvarname[iv] + "_energywisebin.png");

	    }

	  plot_tgraph(hist_timeresvs_ampeff_C2_C3,"defall_energy_"+central_crystal+"_gain" + gain_cut[ig]+ "_timeres_vsAeff_.png");

	}
    }	  
		  
}

void plot_energy_resolution()
{
  string central_crystal;
#ifdef target_C2
  central_crystal = "C2";
#elif defined(target_C3)
  central_crystal = "C3";
#endif

  string fillvar[] = {"E5x5","E3x3","amp_max["+central_crystal+"]","pedestal["+central_crystal+"]"};
  int nvars = sizeof(fillvar)/sizeof(fillvar[0]);
  TH1D * hist_energyfit[4];
  hist_energyfit[0] = new TH1D("E5x5","E5x5",200,0,10000);
  hist_energyfit[1] = new TH1D("E3x3","E3x3",200,0,10000);
  hist_energyfit[2] = new TH1D("amp","amp",200,200,10000);
  hist_energyfit[3] = new TH1D("ped","ped",120,0,20);
  //TH1D * hist_res = new TH1D("dt_C2_MCP1","dt_C2_MCP1",60,0,6.238);

  TH1D *hist_ped = new TH1D("ped","ped",20,0,30);

  TH1D * hist_ampvs_energy[2];
  hist_ampvs_energy[0] = new TH1D("e5x5_energy","E5x5 vs energy",nenergies,0,nenergies);
  hist_ampvs_energy[0]->GetXaxis()->SetTitle("Energy of beam (in GeV)");
  hist_ampvs_energy[0]->GetYaxis()->SetTitle("Mean E5x5 (in ADC)"); ///#sigma(E5x5)/amp_max

  hist_ampvs_energy[1] = new TH1D("e3x3_energy","E3x3 vs energy",nenergies,0,nenergies);
  hist_ampvs_energy[1]->GetXaxis()->SetTitle("Energy of beam (in GeV)");
  hist_ampvs_energy[1]->GetYaxis()->SetTitle("Mean  E3x3 (in ADC)");

  TH1D * hist_sigmavs_energy[2];
  hist_sigmavs_energy[0] = new TH1D("e5x5_energy","E5x5 vs energy",nenergies,0,nenergies);
  hist_sigmavs_energy[0]->GetXaxis()->SetTitle("Energy of beam (in GeV)");
  hist_sigmavs_energy[0]->GetYaxis()->SetTitle("#sigma(E5x5) / mean(E5x5)"); ///#sigma(E5x5)/amp_max

  hist_sigmavs_energy[1] = new TH1D("e3x3_energy","E3x3 vs energy",nenergies,0,nenergies);
  hist_sigmavs_energy[1]->GetXaxis()->SetTitle("Energy of beam (in GeV)");
  hist_sigmavs_energy[1]->GetYaxis()->SetTitle("#sigma(E3x3) / mean(E3x3)");

  TH1D * hist_sigmavs_energy_gain10[2];
  hist_sigmavs_energy_gain10[0] = new TH1D("e5x5_energy_gain10","E5x5 vs energy",nenergies,0,nenergies);
  hist_sigmavs_energy_gain10[0]->GetXaxis()->SetTitle("Energy of beam (in GeV)");
  hist_sigmavs_energy[0]->GetYaxis()->SetTitle("#sigma(E5x5) / mean(E5x5)"); ///#sigma(E5x5)/amp_max

  hist_sigmavs_energy_gain10[1] = new TH1D("e3x3_energy_gain10","E3x3 vs energy",nenergies,0,nenergies);
  hist_sigmavs_energy_gain10[1]->GetXaxis()->SetTitle("Energy of beam (in GeV)");
  hist_sigmavs_energy_gain10[1]->GetYaxis()->SetTitle("#sigma(E3x3) / mean(E3x3)");

  TGraphErrors * tgraph_e3x3_alpha = new TGraphErrors(nenergies+2);
  tgraph_e3x3_alpha->GetXaxis()->SetTitle("For E3x3");
  tgraph_e3x3_alpha->GetYaxis()->SetTitle("Alpha of CB");

  TGraphErrors * tgraph_e3x3_N = new TGraphErrors(nenergies+2);
  tgraph_e3x3_N->GetXaxis()->SetTitle("For E3x3");
  tgraph_e3x3_N->GetYaxis()->SetTitle("N of CB");

  TGraphErrors * tgraph_e5x5_alpha = new TGraphErrors(nenergies+2);
  tgraph_e5x5_alpha->GetXaxis()->SetTitle("For E5x5");
  tgraph_e5x5_alpha->GetYaxis()->SetTitle("Alpha of CB");

  TGraphErrors * tgraph_e5x5_N = new TGraphErrors(nenergies+2);
  tgraph_e5x5_N->GetXaxis()->SetTitle("For E5x5");
  tgraph_e5x5_N->GetYaxis()->SetTitle("N of CB");

  for(int ig =0; ig < 2; ig++)
    {
      for(int ie =0; ie < nenergies; ie++)
  	{
  	  string beam_energy = energies[ie] + " GeV and gain of " + central_crystal +" = "  + gain_cut[ig];
  	  string evt_sel_cuts;

	  evt_sel_cuts = amp_max_cut[ie] + " && gain["+central_crystal+"] == " + gain_cut[ig];
	  cout<<"Starting for energy"<<energies[ie]<<endl;

	  TH1D * hist_fit[int(nvars)];
	  //double resolution_values[3];
	      
	  TH1D * hists_list[2] = {NULL,NULL};
	  //string str_list[3] = {" "};
	  double amp_mean[3] = {0.0};
	  double amp_mean_error[3] = {0.0};

	  double sigma[3] = {0.0};
	  double sigma_error[3] = {0.0};
	      
	  for(int iv =0; iv < nvars; iv++)
	    {
	      //if(iv==2 )
	      //continue;
	      hist_energyfit[iv]->Reset("ICESM");
	      gethist(energies[ie] + "_" + central_crystal + "2x2_"+file_tag+".root",treename,fillvar[iv],evt_sel_cuts,1,hist_energyfit[iv]);
	      //hists_list[iv] = (TH1D*)hist_energyfit[int(min(2,iv))]->Clone();

	      if(hist_energyfit[iv]->Integral()<100)
		continue;

	      int binmax = hist_energyfit[iv]->GetMaximumBin();
	      double xmean = hist_energyfit[iv]->GetXaxis()->GetBinCenter(binmax);
		  
	      double std = hist_energyfit[iv]->GetStdDev();	  
	      //double std = abs(xmean - hist_energyfit[iv]->GetXaxis()->GetBinCenter(hist_energyfit[iv]->FindFirstBinAbove(hist_energyfit[iv]->GetBinContent(binmax)/2.5)));
	      //double std = xmean/10;
		  
	      //double xmin = max(1.0,double(int((xmean - 3*max(0.01,std))*100))/100), xmax = double(int((xmean + 2*max(0.01,std))*100))/100;
	      double xmin = 0.89*xmean, xmax = 1.1*xmean;

	      int nbins = 60;
	      if(iv == 0 && ie ==0)
		{
		  xmin = 0.85*xmean;
		  xmax = 1.1*xmean;
		}
	      if( iv == 2)
		{
		  xmin = max(200.0,xmean - 4*std);
		  xmax = min(xmean + 800.0,xmean + 2.75*std);
		}
	      if( iv == 3)
		{
		  xmin = 0;
		  xmax = 20;
		  nbins = 20;
		} 
	      //double xmin = 0, xmax = 8000;
	      //double xmin = max(0.,xmean - 2*std), xmax = min(6.238,xmean + 2*std);
	      hist_fit[iv] = new TH1D((fillvar[iv] + energies[ie]).c_str(),fillvar[iv].c_str(),nbins,xmin,xmax);
	      //hist_fit[iv]->Draw();


	      cout<<"At energy "<<energies[ie]<<", fit of "<<fillvar[iv]<<" range choosen = "<<xmin<<" - "<<xmax<<" with max at "<<xmean<<" std dev = "<<std<<endl;
	      gethist(energies[ie] + "_" + central_crystal + "2x2_"+file_tag+".root",treename,fillvar[iv],evt_sel_cuts,1,hist_fit[iv]);
	  
	  
	      //remove_zerobincontents(hist_fit[iv]);
	      double * rvalues;
	      if(iv == 0) /// E5x5
		{
		  rvalues = get_resolution_value_crystalball_gauss(hist_fit[iv],beam_energy,xmin,xmax);   ///  get_resolution_value_crystalball   get_resolution_value_doublecrystalball
		  sigma[iv] = rvalues[0]*rvalues[8] + rvalues[12]*(1-rvalues[8]);
		  sigma_error[iv] = rvalues[1]*rvalues[8] + rvalues[13]*(1-rvalues[8]);
		  amp_mean[iv] = rvalues[2]*rvalues[8] + rvalues[10]*(1-rvalues[8]);
		  amp_mean_error[iv] = rvalues[3]*rvalues[8] + rvalues[11]*(1-rvalues[8]);
		  tgraph_e5x5_alpha->SetPoint(ie,ie+1,rvalues[4]);
		  tgraph_e5x5_alpha->SetPointError(ie,0.01,rvalues[5]);
		  tgraph_e5x5_N->SetPoint(ie,ie+1,rvalues[6]);
		  tgraph_e5x5_N->SetPointError(ie,0.01,rvalues[7]);
		}
	  
	      if(iv == 1) /// E3x3
		{
		  rvalues = get_resolution_value_crystalball_gauss(hist_fit[iv],beam_energy,xmin,xmax);
		  sigma[iv] = rvalues[0]*rvalues[8] + rvalues[12]*(1-rvalues[8]);
		  sigma_error[iv] = rvalues[1]*rvalues[8] + rvalues[13]*(1-rvalues[8]);
		  amp_mean[iv] = rvalues[2]*rvalues[8] + rvalues[10]*(1-rvalues[8]);
		  amp_mean_error[iv] = rvalues[3]*rvalues[8] + rvalues[11]*(1-rvalues[8]);
		  tgraph_e3x3_alpha->SetPoint(ie,ie+1,rvalues[4]);
		  tgraph_e3x3_alpha->SetPointError(ie,0.01,rvalues[5]);
		  tgraph_e3x3_N->SetPoint(ie,ie+1,rvalues[6]);
		  tgraph_e3x3_N->SetPointError(ie,0.01,rvalues[7]);
		}
	      if(iv == 2) /// amp_max[C2]
		{
		  rvalues = get_resolution_value_crystalball(hist_fit[iv],beam_energy,xmin,xmax);
		  sigma[iv] = rvalues[0];
		  sigma_error[iv] = rvalues[1];
		  amp_mean[iv] = rvalues[2];
		  amp_mean_error[iv] = rvalues[3];
		}

	      if(iv == 3) 
		{
		  rvalues = get_resolution_value_gauss(hist_fit[iv],beam_energy,xmin,xmax);
		  hist_ped->Add(hist_fit[iv]);
		  //amp_mean[iv] = rvalues[2];
		  //amp_mean_error[iv] = rvalues[3];
		}
	  
	      // if(iv == 0 || iv == 1)
	      // 	{
	      // 	  TH2D *hist2d_energyvsX = new TH2D(("hist2d_" + fillvar[iv] + "vsX" + energies[ie]).c_str(),"hist2d_energyvsX",120,-15,15,120,xmin,xmax);
	      // 	  TH2D *hist2d_energyvsY = new TH2D(("hist2d_" + fillvar[iv] + "vsY" + energies[ie]).c_str(),"hist2d_energyvsY",60,-15,15,60,xmin,xmax);
	      
	      // 	  hist2d_energyvsX->GetXaxis()->SetTitle("X");
	      // 	  hist2d_energyvsX->GetYaxis()->SetTitle(fillvar[iv].c_str());

	      // 	  hist2d_energyvsY->GetXaxis()->SetTitle("Y");
	      // 	  hist2d_energyvsY->GetYaxis()->SetTitle(fillvar[iv].c_str());

	      // 	  get2dhist(energies[ie] + "_" + central_crystal + "2x2_"+file_tag+".root",treename, fillvar[iv] + ":h1X",evt_sel_cuts,1,hist2d_energyvsX);
	      // 	  get2dhist(energies[ie] + "_" + central_crystal + "2x2_"+file_tag+".root",treename,fillvar[iv] + ":h1Y",evt_sel_cuts,1,hist2d_energyvsY);
	      // 	  savehist2d(hist2d_energyvsX,"hist2d_" +fillvar[iv] + "h1X_atengergy" + energies[ie] + "_atgain" + gain_cut[ig],false);	  
	      // 	  savehist2d(hist2d_energyvsY,"hist2d_" +fillvar[iv] + "vsh1Y_atengergy" + energies[ie] + "_atgain" + gain_cut[ig],false);	  

	
	      // 	  TH1D * hists_profile[2] = {hist2d_energyvsX->ProfileX(),hist2d_energyvsY->ProfileX()};
	      
	      // 	  string str_profile[2] = {"on X","on Y"};
	      // 	  plot_1dhists(2,hists_profile,str_profile,"profile_vs_energy" + energies[ie] + "_atgain" + gain_cut[ig] + ".png");
	      // 	}
	    }

	  if(ig == 0)
	    {
	      hist_ampvs_energy[0]->SetBinContent(ie+1,amp_mean[0]);
	      hist_ampvs_energy[0]->SetBinError(ie+1,amp_mean_error[0]);
	      hist_ampvs_energy[0]->GetXaxis()->SetBinLabel(ie+1,energies[ie].c_str());
	  
	      hist_ampvs_energy[1]->SetBinContent(ie+1,amp_mean[1]);
	      hist_ampvs_energy[1]->SetBinError(ie+1,amp_mean_error[1]);
	      hist_ampvs_energy[1]->GetXaxis()->SetBinLabel(ie+1,energies[ie].c_str());

	      hist_sigmavs_energy[0]->SetBinContent(ie+1,sigma[0]/amp_mean[0]);
	      hist_sigmavs_energy[0]->SetBinError(ie+1,sigma[0]/amp_mean[0]*(sigma_error[0]/sigma[0] + amp_mean_error[0]/amp_mean[0]));
	      hist_sigmavs_energy[0]->GetXaxis()->SetBinLabel(ie+1,energies[ie].c_str());
	  
	      hist_sigmavs_energy[1]->SetBinContent(ie+1,sigma[1]/amp_mean[1]);
	      hist_sigmavs_energy[1]->SetBinError(ie+1,sigma[1]/amp_mean[1]*(sigma_error[1]/sigma[1] + amp_mean_error[1]/amp_mean[1]));
	      hist_sigmavs_energy[1]->GetXaxis()->SetBinLabel(ie+1,energies[ie].c_str());
	    }
	  else if(ig == 1)
	    {
	      hist_ampvs_energy[0]->SetBinContent(ie+1, gain_switch_calriabation_factor * amp_mean[0]);
	      hist_ampvs_energy[0]->SetBinError(ie+1, gain_switch_calriabation_factor * amp_mean_error[0]);
	      hist_ampvs_energy[0]->GetXaxis()->SetBinLabel(ie+1,energies[ie].c_str());
	  
	      hist_ampvs_energy[1]->SetBinContent(ie+1, gain_switch_calriabation_factor * amp_mean[1]);
	      hist_ampvs_energy[1]->SetBinError(ie+1, gain_switch_calriabation_factor * amp_mean_error[1]);
	      hist_ampvs_energy[1]->GetXaxis()->SetBinLabel(ie+1,energies[ie].c_str());

	      hist_sigmavs_energy_gain10[0]->SetBinContent(ie+1,sigma[0]/amp_mean[0]);
	      hist_sigmavs_energy_gain10[0]->SetBinError(ie+1,sigma[0]/amp_mean[0]*(sigma_error[0]/sigma[0] + amp_mean_error[0]/amp_mean[0]));
	      hist_sigmavs_energy_gain10[0]->GetXaxis()->SetBinLabel(ie+1,energies[ie].c_str());
	  
	      hist_sigmavs_energy_gain10[1]->SetBinContent(ie+1,sigma[1]/amp_mean[1]);
	      hist_sigmavs_energy_gain10[1]->SetBinError(ie+1,sigma[1]/amp_mean[1]*(sigma_error[1]/sigma[1] + amp_mean_error[1]/amp_mean[1]));
	      hist_sigmavs_energy_gain10[1]->GetXaxis()->SetBinLabel(ie+1,energies[ie].c_str());
	    }
	}
      double * rvalues = get_resolution_value_gauss(hist_ped,"all",0,20);
      
      plot_CB_params(tgraph_e3x3_alpha);
      plot_CB_params(tgraph_e3x3_N);
      plot_CB_params(tgraph_e5x5_alpha);
      plot_CB_params(tgraph_e5x5_N);

    }

  /*TH1D * hists_e3x3[1] = {hist_ampvs_energy[1]};
  string str_e3x3[1] = {" "};
  plot_1dhists_withfit(1,hists_e3x3,str_e3x3,"e3x3_vs_energy_gaincombined.png");
  
  TH1D * hists_e5x5[1] = {hist_ampvs_energy[0]};
  string str_e5x5[1] = {" "};
  plot_1dhists_withfit(1,hists_e5x5,str_e5x5,"e5x5_vs_energy_gaincombined.png");
*/
  
  TH1D * hists_e3x3_res[2] = {hist_sigmavs_energy[1],hist_sigmavs_energy_gain10[1]};
  string str_e3x3_res[2] = {"gain 1 of C2","gain 10 of C2"};
  plot_1dhists_withfit(2,hists_e3x3_res,str_e3x3_res,"e3x3_energy_resolution_vs_energy_gaincombined.png");
  
  TH1D * hists_e5x5_res[2] = {hist_sigmavs_energy[0],hist_sigmavs_energy_gain10[0]};
  string str_e5x5[2] = {"gain 1 of C2","gain 10 of C2"};
  plot_1dhists_withfit(2,hists_e5x5_res,str_e5x5,"e5x5_energy_resolution_vs_energy_gaincombined.png");

  
    /*  RooRealVar x("x", "x", -10, 10);
	RooRealVar mean("mean","mean of gaussian",0,-10,10) ;
	RooRealVar sigma("sigma","width of gaussian",3) ;
	RooGaussian gauss("gauss","gaussian PDF",x,mean,sigma);

	std::unique_ptr<RooDataSet> data{gauss.generate(x, 10000)};
	cout<<get_resolution_value_gauss((TH1D*)data->createHistogram("test",x,Binning(100,-10,10)),"test",-10,10);-*/
}

void basic_plots()
{
  string target_crystal = "C32x2";
  string central_crystal = "C3";
  
  string fillvar[] = {"amp_max[C2]","amp_max[C3]","amp_max[C1]","amp_max[MCP1]","amp_max[MCP2]","gain[C2]","gain[C3]","gain[C1]","R13","R15"};

  int nvars = sizeof(fillvar)/sizeof(fillvar[0]);

  double bin_low[] = {0,0,0,0,0,-0.5,-0.5,-0.5,0,0};
  double bin_up[] = {1000,6000,1000,2000,2000,10.5,10.5,10.5,1,1};
  int nbin[] = {100,100,100,100,100,11,11,11,30,30};

  if(nvars != sizeof(bin_low)/sizeof(bin_low[0]) || nvars != sizeof
     (bin_up)/sizeof(bin_up[0]) || nvars != sizeof(nbin)/sizeof(nbin[0]))
    {cout<<"No of variable and bin information dont match"<<endl;exit(0);}

  for(int ie =10; ie < nenergies; ie++)
    {
      for(int iv =0; iv < nvars; iv++)
	{
	  TH1D* hist;
	  hist = new TH1D(("hist_nocuts_" + fillvar[iv] + energies[ie]).c_str(),("hist" + fillvar[iv]).c_str(),nbin[iv],bin_low[iv],bin_up[iv]);
	  gethist(energies[ie] + "_" + target_crystal + "_"+file_tag+".root",treename,fillvar[iv],"",1,hist);
	  
	  string str = fillvar[iv];
	  plot_hist(hist,str,"energy_" + energies[ie] + target_crystal + "_" + "_no_cuts_" + fillvar[iv] + ".png");
	      
	}
    }

  string fillvar1_2d[] = {"h1X","h2X","h1X","h1Y"};
  string fillvar2_2d[] = {"h1Y","h2Y","amp_max[" + central_crystal +"]","amp_max[" + central_crystal +"]"};

  string cuts[] ={"","","n_h1X>0","n_h1Y>0"};
  
  int nvars1_2d = sizeof(fillvar1_2d)/sizeof(fillvar1_2d[0]);
  int nvars2_2d = sizeof(fillvar2_2d)/sizeof(fillvar2_2d[0]);

  double bin_low1[] = {-15.5,-15.5,-15.5,-15.5};
  double bin_up1[] = {15.5,15.5,15.5,15.5};
  int nbin1[] = {62,62,62,62};

  double bin_low2[] = {-15.5,-15.5,0,0};
  double bin_up2[] = {15.5,15.5,3000,3000};
  int nbin2[] = {62,62,30,30};

  if(nvars1_2d != nvars2_2d || nvars1_2d != sizeof(bin_up1)/sizeof(bin_up1[0]) || nvars1_2d != sizeof(nbin1)/sizeof(nbin1[0]) ||  nvars1_2d != sizeof(nbin2)/sizeof(nbin2[0]) ||  nvars1_2d != sizeof(bin_up2)/sizeof(bin_up2[0]))
    {cout<<"No of variable and bin information dont match for 2d histograms"<<endl;exit(0);}

  for(int ie =10; ie < nenergies; ie++)
    {
      for(int iv =0; iv < nvars1_2d; iv++)
	{
	  TH2D* hist;
	  hist = new TH2D(("hist_nocuts_" + fillvar1_2d[iv] + "_" + fillvar1_2d[iv] + energies[ie]).c_str(),("hist_" + fillvar1_2d[iv] + "_"+ fillvar2_2d[iv]).c_str(),nbin1[iv],bin_low1[iv],bin_up1[iv],nbin2[iv],bin_low2[iv],bin_up2[iv]);
	  hist->GetXaxis()->SetTitle(fillvar1_2d[iv].c_str());  
	  hist->GetYaxis()->SetTitle(fillvar2_2d[iv].c_str());  

	  get2dhist(energies[ie] + "_" + target_crystal + "_"+file_tag+".root",treename,fillvar2_2d[iv] + ":" + fillvar1_2d[iv],cuts[iv],1,hist);
	  
	  savehist2d(hist,"energy_" + energies[ie] + target_crystal + "_" + "_no_cuts_" + fillvar1_2d[iv] + "_" +  fillvar2_2d[iv],false);	  	      
	}
    }


  /* string cuts[] = {"gain[C2]","gain[C3]","gain[C1]"};
  for(int ig = 0; ig < 2; ig ++ )
    {
      for(int ie =0; ie < nenergies; ie++)
	{
	  for(int iv =0; iv < 2; iv++)
	    {
	      if(iv > 2 ) continue;
	      TH1D* hist;
	      hist = new TH1D(("hist_nocuts_" + fillvar[iv] + energies[ie]).c_str(),("hist" + fillvar[iv]).c_str(),nbin[iv],bin_low[iv],bin_up[iv]);
	      gethist(energies[ie] + "_"+file_tag+".root",treename,fillvar[iv],cuts + "=="+ to_string(gain_cut[ig]),1,hist);
	      
	      if(hist->Integra() > 100)
	      {
	      string str = fillvar[iv];
	      plot_hist(hist,str,"energy_" + energies[ie] + "_gaincut_" + gain_cut[ig] + "_" + fillvar[iv] + ".png");
	      }
	    }
	}
	}*/


  ////// for noise calculation
  string ped_nosie_var[4] = {"pedestal[C2]","pedestal[C3]","amp_max[C2]","amp_max[C3]"};
  for(int iv =0; iv < 4; iv++)
    {
      TH1D* hist;
      double xmin = 0, xmax = 8;
      if(iv < 2) {xmin = 8; xmax = 15;}
      if(iv > 1) {xmin = 0.5; xmax = 4.5;}
      hist = new TH1D(("hist_nocuts_" + ped_nosie_var[iv] + "_allenergy").c_str(),("hist" + ped_nosie_var[iv]).c_str(),30,xmin,xmax);
      
      gethist("pedestal_pedestal_digisync_VFE_fix_test_files.root",treename,ped_nosie_var[iv],"",1,hist);

      string str = ped_nosie_var[iv];

      double *r_values = get_resolution_value_gauss(hist,"All energy using pedestral spills",xmin,xmax); //crystalball
      
    }
}

void calculate_intercalritationfactor()
{
  
  string fillvar[] = {"R13","R15","amp_max[C2]"};
  string r13_var = "amp_max[C2] / (amp_max[C1] + amp_max[C2] + amp_max[C3] + amp_max[B1] + amp_max[B2] + amp_max[B3] + amp_max[D1] + amp_max[D2] + amp_max[D3])";
  string r15_var = "amp_max[C2] / (amp_max[C1] + amp_max[C2] + amp_max[C3] + amp_max[C4] + amp_max[C5] + amp_max[A1] + amp_max[A2] + amp_max[A3] + amp_max[A4] + amp_max[A5] + amp_max[B1] + amp_max[B2] + amp_max[B3] + amp_max[B4] + amp_max[B5] + amp_max[D1] + amp_max[D2] + amp_max[D3] + amp_max[D4] + amp_max[D5] + amp_max[E1] + amp_max[E2] + amp_max[E3] + amp_max[E4] + amp_max[E5])";
  string e8_var = "(amp_max[C1] + amp_max[C3] + amp_max[B1] + amp_max[B2] + amp_max[B3] + amp_max[D1] + amp_max[D2] + amp_max[D3])";
  string e24_var = "(amp_max[C1] + amp_max[C3] + amp_max[C4] + amp_max[C5] + amp_max[A1] + amp_max[A2] + amp_max[A3] + amp_max[A4] + amp_max[A5] + amp_max[B1] + amp_max[B2] + amp_max[B3] + amp_max[B4] + amp_max[B5] + amp_max[D1] + amp_max[D2] + amp_max[D3] + amp_max[D4] + amp_max[D5] + amp_max[E1] + amp_max[E2] + amp_max[E3] + amp_max[E4] + amp_max[E5])";
  string e9_var = "(amp_max[C1]+ amp_max[C2] + amp_max[C3] + amp_max[B1] + amp_max[B2] + amp_max[B3] + amp_max[D1] + amp_max[D2] + amp_max[D3])";
  string e25_var = "(amp_max[C1] + amp_max[C2] + amp_max[C3] + amp_max[C4] + amp_max[C5] + amp_max[A1] + amp_max[A2] + amp_max[A3] + amp_max[A4] + amp_max[A5] + amp_max[B1] + amp_max[B2] + amp_max[B3] + amp_max[B4] + amp_max[B5] + amp_max[D1] + amp_max[D2] + amp_max[D3] + amp_max[D4] + amp_max[D5] + amp_max[E1] + amp_max[E2] + amp_max[E3] + amp_max[E4] + amp_max[E5])";

  string e24_var_gain10 = "( (gain[C1] == 1 ? 0 : amp_max[C1]) + (gain[C3] == 1 ? 0 : amp_max[C3]) + (gain[C4] == 1 ? 0 : amp_max[C4]) + (gain[C5] == 1 ? 0 : amp_max[C5]) + (gain[A1] == 1 ? 0 : amp_max[A1]) + (gain[A2] == 1 ? 0 : amp_max[A2]) + (gain[A3] == 1 ? 0 : amp_max[A3]) + (gain[A4] == 1 ? 0 : amp_max[A4]) + (gain[A5] == 1 ? 0 : amp_max[A5]) + (gain[B1] == 1 ? 0 : amp_max[B1]) + (gain[B2] == 1 ? 0 : amp_max[B2]) + (gain[B3] == 1 ? 0 : amp_max[B3]) + (gain[B4] == 1 ? 0 : amp_max[B4]) + (gain[B5] == 1 ? 0 : amp_max[B5]) + (gain[D1] == 1 ? 0 : amp_max[D1]) + (gain[D2] == 1 ? 0 : amp_max[D2]) + (gain[D3] == 1 ? 0 : amp_max[D3]) + (gain[D4] == 1 ? 0 : amp_max[D4]) + (gain[D5] == 1 ? 0 : amp_max[D5]) + (gain[E1] == 1 ? 0 : amp_max[E1]) + (gain[E2] == 1 ? 0 : amp_max[E2]) + (gain[E3] == 1 ? 0 : amp_max[E3]) + (gain[E4] == 1 ? 0 : amp_max[E4]) + (gain[E5] == 1 ? 0 : amp_max[E5]) )";

  string gain_e8_cut = "(gain[C1] == 1 && gain[C3] == 1 && gain[B1] == 1 && gain[B2] == 1 && gain[B3] == 1 && gain[D1] == 1 && gain[D2] == 1 && gain[D3] == 1)";
  string gain_e24_cut = "(gain[C1] == 1 && gain[C3] == 1 && gain[C4] == 1 && gain[C5] == 1 && gain[A1] == 1 && gain[A2] == 1 && gain[A3] == 1 && gain[A4] == 1 && gain[A5] == 1 && gain[B1] == 1 && gain[B2] == 1 && gain[B3] == 1 && gain[B4] == 1 && gain[B5] == 1 && gain[D1] == 1 && gain[D2] == 1 && gain[D3] == 1 && gain[D4] == 1 && gain[D5] == 1 && gain[E1] == 1 && gain[E2] == 1 && gain[E3] == 1 && gain[E4] == 1 && gain[E5] == 1)";

  string filename = "100_"+file_tag+".root";
  string filename_gain1 = "gaitchsiwtch.root";
  string filename_gain10 = "forced_gain_100gev.root";

  double xmin = 500, xmax = 3000;
  int nbins = 30;
 
  vector<tgraphpair> calfact_r13_array;
  vector<tgraphpair> calfact_r19_array;
  vector<tgraphpair> calfact_e8_array;
  vector<tgraphpair> calfact_e24_array;

  vector<double>r13cuts;
  vector<double>r15cuts;
  vector<double>e8cuts;
  vector<double>e24cuts;
  
  int ncuts = 5;
  float r13cut_startval = 0.6;  float r13cut_diffval = 0.04;
  float r15cut_startval = 0.6;  float r15cut_diffval = 0.04;
  float e8_cut_startval = 1000;  float e8_cut_diffval = 200;
  float e24_cut_startval = 1000;  float e24_cut_diffval = 200;

  TGraphErrors * tgraph_cal_factor_r13 = new TGraphErrors(ncuts*1);
  tgraph_cal_factor_r13->GetXaxis()->SetTitle("Amp max when gain = 1");
  tgraph_cal_factor_r13->GetYaxis()->SetTitle("Amp max when gain = 10");

  TGraphErrors * tgraph_cal_factor_r19 = new TGraphErrors(ncuts*1);
  tgraph_cal_factor_r19->GetXaxis()->SetTitle("Amp max when gain = 1");
  tgraph_cal_factor_r19->GetYaxis()->SetTitle("Amp max when gain = 10");

  TGraphErrors * tgraph_cal_factor_e8 = new TGraphErrors(ncuts*1);
  tgraph_cal_factor_e8->GetXaxis()->SetTitle("Amp max when gain = 1");
  tgraph_cal_factor_e8->GetYaxis()->SetTitle("Amp max when gain = 10");

  TGraphErrors * tgraph_cal_factor_e24 = new TGraphErrors(ncuts*1);
  tgraph_cal_factor_e24->GetXaxis()->SetTitle("Amp max when gain = 1");
  tgraph_cal_factor_e24->GetYaxis()->SetTitle("Amp max when gain = 10");
 
   ///////////////// Plot basic histogram to study

  string basicvar[8] = {"amp_max[C2]",r13_var,r15_var,e8_var,e24_var,e9_var,e25_var,e24_var_gain10};
  string basicvarname[8] = {"amp_max[C2]","R13","R15","E8","E24","E9","E25","E24(gain-10)"};
    
  double bin_low[8] = {0,0,0,0,0,0,0,0};
  double bin_up[8] = {xmax,1,1,3000,3000,4500,4500,3000};
  int nbin[8] = {int(nbins*xmax/(xmax-xmin)), int(1.0/r15cut_diffval),  int(1.0/r15cut_diffval),  int(bin_up[3]/e8_cut_diffval),  int(bin_up[4]/e24_cut_diffval), 30, 30, 30};
  
  for(int iv =0; iv < 8; iv++)
    {
      string str = basicvarname[iv];
      
      TH1D* hist1 = new TH1D(("hist_nocuts_" + basicvarname[iv]).c_str(),("hist" + basicvarname[iv]).c_str(),nbin[iv],bin_low[iv],bin_up[iv]);
      gethist(filename_gain1,calc_treename,basicvar[iv],"",1,hist1);
      plot_hist(hist1,str,"calribation_no_cuts_gain1_" + basicvarname[iv] + ".png");

      hist1->Reset("ICES");
      if(iv ==3) gethist(filename_gain1,calc_treename,basicvar[iv],gain_e8_cut,1,hist1);
      else if(iv ==4) gethist(filename_gain1,calc_treename,basicvar[iv],gain_e24_cut,1,hist1);
      if( iv == 3 || iv ==4 ) plot_hist(hist1,str,"calribation_gaincut_gain1_" + basicvarname[iv] + ".png");
      
      TH1D* hist10 = new TH1D(("hist_nocuts_" + basicvarname[iv]).c_str(),("hist" + basicvarname[iv]).c_str(),nbin[iv],bin_low[iv],bin_up[iv]);
      gethist(filename_gain10,calc_treename,basicvar[iv],"",1,hist10);
      plot_hist(hist10,str,"calribation_no_cuts_gain10_" + basicvarname[iv] + ".png");

      hist10->Reset("ICES");
      if(iv ==3)  gethist(filename_gain10,calc_treename,basicvar[iv],gain_e8_cut,1,hist10);
      else if(iv ==4)  gethist(filename_gain10,calc_treename,basicvar[iv],gain_e24_cut,1,hist10);
      if( iv == 3 || iv ==4 ) plot_hist(hist10,str,"calribation_gaincut_gain10_" + basicvarname[iv] + ".png");

      if(iv == 5 )
	{
	  double *r_values;
	  
	  TH1D* hist1 = new TH1D(("hist_gaincut_" + basicvarname[iv]).c_str(),("hist" + basicvarname[iv]).c_str(),60,1000,4000);
	  gethist(filename_gain1,calc_treename,basicvar[iv],gain_e8_cut,1,hist1);
	  plot_hist(hist1,str,"calribation_distfore8_gain1_" + basicvarname[iv] + ".png");
	  
	  TH1D* hist10 = new TH1D(("hist_gaincut_" + basicvarname[iv]).c_str(),("hist" + basicvarname[iv]).c_str(),60,1000,4000);
	  gethist(filename_gain10,calc_treename,basicvar[iv],gain_e8_cut,1,hist10);
	  plot_hist(hist10,str,"calribation_distfore8_gain10_" + basicvarname[iv] + ".png");

	  str = basicvarname[0];
	    
	  TH1D* hist12 = new TH1D(("hist_gaincut_fore9_" + basicvarname[0]).c_str(),("hist_gaincut_fore9_" + basicvarname[0]).c_str(),60,1000,3000);
	  gethist(filename_gain1,calc_treename,basicvar[0],gain_e8_cut,1,hist12);
	  r_values = get_resolution_value_crystalball(hist12,"100 GeV (Gain calribation). Gain of C2 = 1 ",1000,3000);
			      
	  plot_hist(hist12,str,"calribation_distfore8_gain1_" + basicvarname[0] + ".png");
	  
	  TH1D* hist102 = new TH1D(("hist_gaincut_fore9_" + basicvarname[0]).c_str(),("hist_gaincut_fore9_" + basicvarname[0]).c_str(),60,1000,3000);
	  gethist(filename_gain10,calc_treename,basicvar[0],gain_e8_cut,1,hist102);
	  r_values = get_resolution_value_crystalball(hist102,"100 GeV (Gain calribation). Gain of C2 = 10 ",1000,3000);
	  plot_hist(hist102,str,"calribation_distfore8_gain10_" + basicvarname[0] + ".png");

	}

      if(iv == 6 )
	{
	  double *r_values;
	  int binmax;
	  double xmean;
	  
	  TH1D* hist1 = new TH1D(("hist_gaincut_" + basicvarname[iv]).c_str(),("hist" + basicvarname[iv]).c_str(),60,1000,4000);
	  gethist(filename_gain1,calc_treename,basicvar[iv],gain_e24_cut,1,hist1);
	  plot_hist(hist1,str,"calribation_distfore24_gain1_" + basicvarname[iv] + ".png");
	  
	  TH1D* hist10 = new TH1D(("hist_gaincut_" + basicvarname[iv]).c_str(),("hist" + basicvarname[iv]).c_str(),60,1000,4000);
	  gethist(filename_gain10,calc_treename,basicvar[iv],gain_e24_cut,1,hist10);
	  plot_hist(hist10,str,"calribation_distfore24_gain10_" + basicvarname[iv] + ".png");

	  str = basicvarname[0];
	  
	  TH1D* hist12 = new TH1D(("hist_gaincut_fore25_" + basicvarname[0]).c_str(),("hist_gaincut_fore25_" + basicvarname[0]).c_str(),60,1000,3000);
	  gethist(filename_gain1,calc_treename,basicvar[0],gain_e24_cut,1,hist12);
	  r_values = get_resolution_value_crystalball(hist12,"100 GeV (Gain calribation). Gain of C2 = 1 ",1000,3000);
	  plot_hist(hist12,str,"calribation_distfore24_gain1_" + basicvarname[0] + ".png");

	  TH1D* hist102 = new TH1D(("hist_gaincut_fore25_" + basicvarname[0]).c_str(),("hist_gaincut_fore25_" + basicvarname[0]).c_str(),60,1000,3000);
	  gethist(filename_gain10,calc_treename,basicvar[0],gain_e24_cut,1,hist102);
	  r_values = get_resolution_value_crystalball(hist102,"100 GeV (Gain calribation). Gain of C2 = 10 ",1000,3000);
	  plot_hist(hist102,str,"calribation_distfore24_gain10_" + basicvarname[0] + ".png");

	  binmax = hist1->GetMaximumBin();
	  xmean = hist1->GetXaxis()->GetBinCenter(binmax);
	  cout<<"For gain 1, mean = "<<xmean<<endl;

	  binmax = hist10->GetMaximumBin();
	  xmean = hist10->GetXaxis()->GetBinCenter(binmax);
	  cout<<"For gain 10, mean = "<<xmean<<endl;

	  string e25_var_corr = "(amp_max[C1] + amp_max[C2] + (3025.0/2575.0)*amp_max[C3] + amp_max[C4] + amp_max[C5] + amp_max[A1] + amp_max[A2] + amp_max[A3] + amp_max[A4] + amp_max[A5] + amp_max[B1] + amp_max[B2] + amp_max[B3] + amp_max[B4] + amp_max[B5] + amp_max[D1] + amp_max[D2] + amp_max[D3] + amp_max[D4] + amp_max[D5] + amp_max[E1] + amp_max[E2] + amp_max[E3] + amp_max[E4] + amp_max[E5])";
	  basicvar[iv] = e25_var_corr;
	  str = "E25 with calribation factor applied on E25";
	  TH1D* hist1_corr = new TH1D(("hist_gaincut_factor_corrected_" + basicvarname[iv]).c_str(),("hist" + basicvarname[iv]).c_str(),60,1000,4000);
	  gethist(filename_gain10,calc_treename,basicvar[iv],gain_e24_cut,1,hist1_corr);
	  plot_hist(hist1_corr,str,"calribation_corrected_distfore24_gain10_" + basicvarname[iv] + ".png");
	  
	}

    }

  ///////////////////  Plot tgraph for calribation factor
  
  // for(int ir = 0; ir < ncuts; ir++)
  //   {
  //     double sigma[8], sigma_error[8], amp_mean[8], amp_mean_error[8];

  //     //// Using R13 variables
  //     float lcut13 = r13cut_startval + ir*r13cut_diffval;
  //     float hcut13 = lcut13 + r13cut_diffval;

  //     string sel_conditions13,evt_sel_cut13;
  //     sel_conditions13 = to_string(lcut13) + " < R13 < " + to_string(hcut13) + " and gain of C2 = 1";
  //     evt_sel_cut13    = to_string(lcut13) + " < " + r13_var  + " && " + r13_var + " < " + to_string(hcut13) + " && gain[C2] == 1";
      
  //     TH1D* hist1;      
  //     hist1 = new TH1D(("hist_" + fillvar[2] + "_r13gain1").c_str(),("hist_" + fillvar[2] + "_r13_gain1").c_str(),nbins,xmin,xmax);
  //     gethist(filename_gain1,calc_treename,fillvar[2],evt_sel_cut13,1,hist1);
		       
  //     double* rvalues = get_resolution_value_crystalball(hist1,sel_conditions13,xmin,xmax);
  //     sigma[0] = rvalues[0];
  //     sigma_error[0] = rvalues[1];
  //     amp_mean[0] = rvalues[2];
  //     amp_mean_error[0] = rvalues[3];
      
  //     sel_conditions13 = to_string(lcut13) + " < R13 < " + to_string(hcut13) + " and gain of C2 = 10";
  //     evt_sel_cut13    = to_string(lcut13) + " < " + r13_var  + " && " + r13_var + " < " + to_string(hcut13) + " && gain[C2] == 10";
  //     TH1D* hist10;
  //     hist10 = new TH1D(("hist_" + fillvar[2] + "_r13gain10").c_str(),("hist_" + fillvar[2] + "_r13_gain10").c_str(),nbins,xmin,xmax);
  //     gethist(filename_gain10,calc_treename,fillvar[2],evt_sel_cut13,1,hist10);
      
  //     double* rvalues2 = get_resolution_value_crystalball(hist10,sel_conditions13,xmin,xmax);
  //     sigma[1] = rvalues2[0];
  //     sigma_error[1] = rvalues2[1];
  //     amp_mean[1] = rvalues2[2];
  //     amp_mean_error[1] = rvalues2[3];

  //     //// Using R15 variables     
  //     float lcut15 = r15cut_startval + ir*r15cut_diffval;
  //     float hcut15 = lcut15 + r15cut_diffval;

  //     string sel_conditions15,evt_sel_cut15;
  //     sel_conditions15 = to_string(lcut15) + " < R15 < " + to_string(hcut15) + " and gain of C2 = 1";
  //     evt_sel_cut15    = to_string(lcut15) + " < " + r15_var  + " && " + r15_var + " < " + to_string(hcut15) + " && gain[C2] == 1";
      
  //     TH1D* hist15;      
  //     hist15 = new TH1D(("hist_" + fillvar[2] + "_r15gain1").c_str(),("hist_" + fillvar[2] + "_r15_gain1").c_str(),nbins,xmin,xmax);
  //     gethist(filename_gain1,calc_treename,fillvar[2],evt_sel_cut15,1,hist15);
		       
  //     double* rvalues3 = get_resolution_value_crystalball(hist15,sel_conditions15,xmin,xmax);
  //     sigma[2] = rvalues3[0];
  //     sigma_error[2] = rvalues3[1];
  //     amp_mean[2] = rvalues3[2];
  //     amp_mean_error[2] = rvalues3[3];
      

  //     sel_conditions15 = to_string(lcut15) + " < R15 < " + to_string(hcut15) + " and gain of C2 = 10";
  //     evt_sel_cut15    = to_string(lcut15) + " < " + r15_var  + " && " + r15_var + " < " + to_string(hcut15) + " && gain[C2] == 10";
  //     TH1D* hist105;
  //     hist105 = new TH1D(("hist_" + fillvar[2] + "_r15gain10").c_str(),("hist_" + fillvar[2] + "_r15_gain10").c_str(),nbins,xmin,xmax);
  //     gethist(filename_gain10,calc_treename,fillvar[2],evt_sel_cut15,1,hist105);
      
  //     double* rvalues4 = get_resolution_value_crystalball(hist105,sel_conditions15,xmin,xmax);
  //     sigma[3] = rvalues4[0];
  //     sigma_error[3] = rvalues4[1];
  //     amp_mean[3] = rvalues4[2];
  //     amp_mean_error[3] = rvalues4[3];

  //     //// Using E8 variables
  //     float lcut_e8 = e8_cut_startval + ir*e8_cut_diffval;
  //     float hcut_e8 = lcut_e8 + e8_cut_diffval;

  //     string sel_conditions_e8,evt_sel_cut_e8;
  //     sel_conditions_e8 = to_string(lcut_e8) + " < E8 < " + to_string(hcut_e8) + " and gain of C2 = 1";
  //     evt_sel_cut_e8    = to_string(lcut_e8) + " < " + e8_var  + " && " + e8_var + " < " + to_string(hcut_e8) + " && gain[C2] == 1";
      
  //     TH1D* hist1_e8;      
  //     hist1_e8 = new TH1D(("hist_" + fillvar[2] + "_e8_gain1").c_str(),("hist_" + fillvar[2] + "_e8_gain1").c_str(),nbins,xmin,xmax);
  //     gethist(filename_gain1,calc_treename,fillvar[2],evt_sel_cut_e8,1,hist1_e8);
		       
  //     double* rvalues5 = get_resolution_value_crystalball(hist1_e8,sel_conditions_e8,xmin,xmax);
  //     sigma[4] = rvalues5[0];
  //     sigma_error[4] = rvalues5[1];
  //     amp_mean[4] = rvalues5[2];
  //     amp_mean_error[4] = rvalues5[3];
      
  //     sel_conditions_e8 = to_string(lcut_e8) + " < E8 < " + to_string(hcut_e8) + " and gain of C2 = 10";
  //     evt_sel_cut_e8    = to_string(lcut_e8) + " < " + e8_var  + " && " + e8_var + " < " + to_string(hcut_e8) + " && gain[C2] == 10";
  //     TH1D* hist10_e8;
  //     hist10_e8 = new TH1D(("hist_" + fillvar[2] + "_e8_gain10").c_str(),("hist_" + fillvar[2] + "_e8_gain10").c_str(),nbins,xmin,xmax);
  //     gethist(filename_gain10,calc_treename,fillvar[2],evt_sel_cut_e8,1,hist10_e8);
      
  //     double* rvalues6 = get_resolution_value_crystalball(hist10_e8,sel_conditions_e8,xmin,xmax);

  //     sigma[5] = rvalues6[0];
  //     sigma_error[5] = rvalues6[1];
  //     amp_mean[5] = rvalues6[2];
  //     amp_mean_error[5] = rvalues6[3];

  //     //// Using E24 variables
      
  //     float lcut_e24 = e24_cut_startval + ir*e24_cut_diffval;
  //     float hcut_e24 = lcut_e24 + e24_cut_diffval;

  //     string sel_conditions_e24,evt_sel_cut_e24;
  //     sel_conditions_e24 = to_string(lcut_e24) + " < E24 < " + to_string(hcut_e24) + " and gain of C2 = 1";
  //     evt_sel_cut_e24    = to_string(lcut_e24) + " < " + e24_var  + " && " + e24_var + " < " + to_string(hcut_e24) + " && gain[C2] == 1";

  //     TH1D* hist1_e24;      
  //     hist1_e24 = new TH1D(("hist_" + fillvar[2] + "_e24_gain1").c_str(),("hist_" + fillvar[2] + "_e24_gain1").c_str(),nbins,xmin,xmax);
  //     gethist(filename_gain1,calc_treename,fillvar[2],evt_sel_cut_e24,1,hist1_e24);

  //     double* rvalues7 = get_resolution_value_crystalball(hist1_e24,sel_conditions_e24,xmin,xmax);
  //     sigma[6] = rvalues7[0];
  //     sigma_error[6] = rvalues7[1];
  //     amp_mean[6] = rvalues7[2];
  //     amp_mean_error[6] = rvalues7[3];
      
  //     sel_conditions_e24 = to_string(lcut_e24) + " < E24 < " + to_string(hcut_e24) + " and gain of C2 = 10";
  //     evt_sel_cut_e24    = to_string(lcut_e24) + " < " + e24_var  + " && " + e24_var + " < " + to_string(hcut_e24) + " && gain[C2] == 10";
  //     TH1D* hist10_e24;
  //     hist10_e24 = new TH1D(("hist_" + fillvar[2] + "_e24_gain10").c_str(),("hist_" + fillvar[2] + "_e24_gain10").c_str(),nbins,xmin,xmax);
  //     gethist(filename_gain10,calc_treename,fillvar[2],evt_sel_cut_e24,1,hist10_e24);
      
  //     double* rvalues8 = get_resolution_value_crystalball(hist10_e24,sel_conditions_e24,xmin,xmax);
  //     sigma[7] = rvalues8[0];
  //     sigma_error[7] = rvalues8[1];
  //     amp_mean[7] = rvalues8[2];
  //     amp_mean_error[7] = rvalues8[3];
      
  //     //delete rvalues,rvalues2,rvalues3,rvalues4,rvalues5,rvalues6,rvalues7,rvalues8;

  //     ///// Make pair of points for tgraph
      
  //     tgraphpair t13;
  //     t13.x = amp_mean[0];
  //     t13.y = amp_mean[1];
  //     t13.xerror = amp_mean_error[0];
  //     t13.yerror = amp_mean_error[1];
  //     if(hist1->Integral() > 100 && hist10->Integral() > 100)
  // 	{
  // 	  r13cuts.push_back(lcut13);
  // 	  calfact_r13_array.push_back(t13);
  // 	}

  //     tgraphpair t19;
  //     t19.x = amp_mean[2];
  //     t19.y = amp_mean[3];
  //     t19.xerror = amp_mean_error[2];
  //     t19.yerror = amp_mean_error[3];
  //     if(hist15->Integral() > 100 && hist105->Integral() > 100)
  // 	{
  // 	  r15cuts.push_back(lcut15);
  // 	  calfact_r19_array.push_back(t19);
  // 	}
      
  //     tgraphpair te8;
  //     te8.x = amp_mean[4];
  //     te8.y = amp_mean[5];
  //     te8.xerror = amp_mean_error[4];
  //     te8.yerror = amp_mean_error[5];
  //     if(hist1_e8->Integral() > 100 && hist10_e8->Integral() > 100)
  // 	{
  // 	  e8cuts.push_back(lcut_e8);
  // 	  calfact_e8_array.push_back(te8);
  // 	}

  //     tgraphpair te24;
  //     te24.x = amp_mean[6];
  //     te24.y = amp_mean[7];
  //     te24.xerror = amp_mean_error[6];
  //     te24.yerror = amp_mean_error[7];
  //     if(hist1_e24->Integral() > 100 && hist10_e24->Integral() > 100)
  // 	{
  // 	  e24cuts.push_back(lcut_e24);
  // 	  calfact_e24_array.push_back(te24);
  // 	}

  //   }
  // cout<<"End plotting of fir histograms"<<endl;

  // sort(calfact_r13_array.begin(),calfact_r13_array.end(), comparestruct);	
  // sort(calfact_r19_array.begin(),calfact_r19_array.end(), comparestruct);	
  // sort(calfact_e8_array.begin(),calfact_e8_array.end(), comparestruct);	
  // sort(calfact_e24_array.begin(),calfact_e24_array.end(), comparestruct);	
  
  // int n =0 ;
  // for(tgraphpair values : calfact_r13_array)
  //   {
  //     tgraph_cal_factor_r13->SetPoint(n, values.x, values.y);
  //     tgraph_cal_factor_r13->SetPointError(n, values.xerror, values.yerror);
  //     n++;
  //   }
  // n = 0;
  // for(tgraphpair values : calfact_r19_array)
  //   {
  //     tgraph_cal_factor_r19->SetPoint(n, values.x, values.y);
  //     tgraph_cal_factor_r19->SetPointError(n, values.xerror, values.yerror);
  //     n++;
  //   }
  // n = 0;
  // for(tgraphpair values : calfact_e8_array)
  //   {
  //     tgraph_cal_factor_e8->SetPoint(n, values.x, values.y);
  //     tgraph_cal_factor_e8->SetPointError(n, values.xerror, values.yerror);
  //     n++;
  //   }
  // n = 0;
  // for(tgraphpair values : calfact_e24_array)
  //   {
  //     tgraph_cal_factor_e24->SetPoint(n, values.x, values.y);
  //     tgraph_cal_factor_e24->SetPointError(n, values.xerror, values.yerror);
  //     n++;
  //   }

  // plot_tgraph(tgraph_cal_factor_r13,"calcibration_factor_using_r13_energy100.png");
  // plot_tgraph(tgraph_cal_factor_r19,"calcibration_factor_using_r19_energy100.png");
  // plot_tgraph(tgraph_cal_factor_e8,"calcibration_factor_using_e8_energy100.png");
  // plot_tgraph(tgraph_cal_factor_e24,"calcibration_factor_using_e24_energy100.png");

  // cout<<"\nLower values for r13 cuts are :-   ";
  // for(int ir = 0; ir < (int)r13cuts.size(); ir++)
  //  cout<<r13cuts[ir]<<", ";
  // cout<<" with width = "<< r13cut_diffval;
  
  // cout<<"\nLower values for r15 cuts are :-   ";
  // for(int ir = 0; ir < (int)r15cuts.size(); ir++)
  //  cout<<r15cuts[ir]<<", ";
  // cout<<" with width = "<< r15cut_diffval;
 
  // cout<<"\nLower values for E8 cuts are :-   ";
  // for(int ir = 0; ir < (int)e8cuts.size(); ir++)
  //  cout<<e8cuts[ir]<<", ";
  // cout<<" with width = "<< e8_cut_diffval;
  
  // cout<<"\nLower values for E24 cuts are :-   ";
  // for(int ir = 0; ir < (int)e24cuts.size(); ir++)
  //  cout<<e24cuts[ir]<<", ";
  // cout<<" with width = "<< e24_cut_diffval;
  
}


void plot_position_resolution()
{

  string target_crystal = "C32x2"; 
  string energy_weighted_posX = "Sum$((max(8,amp_max_crystal[])-8)*posX_crystal[])/Sum$((max(8,amp_max_crystal[])-8))";
  string energy_weighted_posY = "Sum$((max(8,amp_max_crystal[])-8)*posY_crystal[])/Sum$((max(8,amp_max_crystal[])-8))";
  string logenergy_weighted_posX = "Sum$(log((max(8,amp_max_crystal[])-7))*posX_crystal[])/Sum$(log((max(8,amp_max_crystal[])-7)))";
  string logenergy_weighted_posY = "Sum$(log((max(8,amp_max_crystal[])-7))*posY_crystal[])/Sum$(log((max(8,amp_max_crystal[])-7)))";

  //string energy_weighted_posX = "Sum$(amp_max_crystal[]*posX_crystal[])/Sum$(amp_max_crystal)";
  //string energy_weighted_posY = "Sum$(amp_max_crystal[]*posY_crystal[])/Sum$(amp_max_crystal)";
  //string logenergy_weighted_posX = "Sum$(log(amp_max_crystal[])*posX_crystal[])/Sum$(log(amp_max_crystal))";
  //string logenergy_weighted_posY = "Sum$(log(amp_max_crystal[])*posY_crystal[])/Sum$(log(amp_max_crystal))";

  string fillvar[] = {"h1X","h1Y","trackX","trackY",energy_weighted_posX,energy_weighted_posY,logenergy_weighted_posX,logenergy_weighted_posY};
  string vartitle[] = {"h1X","h1Y","trackX","trackY","energy weighted X","energy weighted Y","log energy weighted X","log energy weighted Y"};
  string cuts[] = {"n_h1X>0","n_h1Y>0","n_tracks>0","","n_tracks>0","","","",""};
  
  int nvars = sizeof(fillvar)/sizeof(fillvar[0]);
  
  double bin_low[] = {-15,-15,-15,-15,-15,-15,-15,-15};
  double bin_up[] = {15,15,15,15,15,15,15,15};
  int nbin[] = {30,30,30,30,30,30,30,30};
  
  if(nvars != sizeof(bin_low)/sizeof(bin_low[0]) || nvars != sizeof
     (bin_up)/sizeof(bin_up[0]) || nvars != sizeof(nbin)/sizeof(nbin[0]))
    {cout<<"No of variable and bin information dont match"<<endl;exit(0);}

  for(int ie =0; ie < nenergies; ie++)
    {
      for(int iv =0; iv < nvars; iv++)
	{
	  TH1D* hist;
	  hist = new TH1D(("hist_nocuts_" + fillvar[iv] + energies[ie]).c_str(),("hist" + fillvar[iv]).c_str(),nbin[iv],bin_low[iv],bin_up[iv]);
	  
	  if(iv > 3) cuts[iv] = amp_max_cut[ie];
	  gethist(energies[ie] + "_" + target_crystal + "_"+file_tag+".root",treename,fillvar[iv],cuts[iv],1,hist);
	  hist->GetXaxis()->SetTitle(vartitle[iv].c_str());

	  
	  string str = vartitle[iv];
	  plot_hist(hist,str,"for_pos_res_energy_" + energies[ie] + "_C3_" + vartitle[iv] + ".png");
	      
	}
    }
  string fillvar1_2d[] = {"trackX","trackY",energy_weighted_posX,energy_weighted_posY,energy_weighted_posX,energy_weighted_posY,logenergy_weighted_posX,logenergy_weighted_posY,logenergy_weighted_posX,logenergy_weighted_posY};
  string fillvar2_2d[] = {"h1X","h1Y","h1X","h1Y","trackX","trackY","h1X","h1Y","trackX","trackY"};
  string cuts_2d[] = {"n_h1X>0 && n_tracks>0","n_h1Y>0 && n_tracks>0","n_h1X>0","n_h1Y>0","n_tracks>0","n_tracks>0","n_h1X>0","n_h1Y>0","n_tracks>0","n_tracks>0"};

  string titlevar1[] = {"trackX","trackY","Energy weighted X","Energy weighted Y","Energy weighted X","Energy weighted Y","log(energy) weighted X","log(energy) weighted Y","log(energy) weighted X","log(energy) weighted Y"};
  int nvars1_2d = sizeof(fillvar1_2d)/sizeof(fillvar1_2d[0]);
  int nvars2_2d = sizeof(fillvar2_2d)/sizeof(fillvar2_2d[0]);
  
  double bin_low1[] = {-15.5,-15.5,-15.5,-15.5,-15.5,-15.5,-15.5,-15.5,-15.5,-15.5};
  double bin_up1[] = {15.5,15.5,15.5,15.5,15.5,15.5,15.5,15.5,15.5,15.5};
  int nbin1[] = {63,63,63,63,63,63,63,63,63,63};

  double bin_low2[] = {-15.5,-15.5,-15.5,-15.5,-15.5,-15.5,-15.5,-15.5,-15.5,-15.5};
  double bin_up2[] = {15.5,15.5,15.5,15.5,15.5,15.5,15.5,15.5,15.5,15.5};
  int nbin2[] = {63,63,63,63,63,63,63,63,63,63};
  
  if(nvars1_2d != nvars2_2d || nvars1_2d != sizeof(bin_up1)/sizeof(bin_up1[0]) || nvars1_2d != sizeof(nbin1)/sizeof(nbin1[0]) ||  nvars1_2d != sizeof(nbin2)/sizeof(nbin2[0]) ||  nvars1_2d != sizeof(bin_up2)/sizeof(bin_up2[0]))
    {cout<<"No of variable and bin information dont match for 2d histograms"<<endl;exit(0);}

  for(int ie =0; ie < nenergies; ie++)
    {
      for(int iv =0; iv < nvars1_2d; iv++)
	{
	  TH2D* hist;
	  hist = new TH2D(("hist_nocuts_" + fillvar1_2d[iv] + "_" + fillvar2_2d[iv] + energies[ie]).c_str(),("hist_" + fillvar1_2d[iv] + "_"+ fillvar2_2d[iv]).c_str(),nbin1[iv],bin_low1[iv],bin_up1[iv],nbin2[iv],bin_low2[iv],bin_up2[iv]);
	  hist->GetYaxis()->SetTitle(fillvar2_2d[iv].c_str());  
	  hist->GetXaxis()->SetTitle(titlevar1[iv].c_str());
	  
	  
	  get2dhist(energies[ie] + "_" + target_crystal + "_"+file_tag+".root",treename,fillvar2_2d[iv] + ":" + fillvar1_2d[iv],"",1,hist);
	  
	  savehist2d(hist,"For_pos_res_energy_" + energies[ie] + target_crystal + "_" + "_no_cuts_" + fillvar2_2d[iv] + "_" +  titlevar1[iv],false);
	  
	  string xtitle;
	  xtitle = titlevar1[iv];
	  
	  plot_hist(hist->ProfileX(),xtitle,"For_pos_res_energy_" + energies[ie] + target_crystal + "_" + "_profile" + fillvar2_2d[iv] + "_inbinsof_" +  titlevar1[iv]+".png");

	  xtitle = fillvar1_2d[iv];
	  plot_hist(hist->ProfileY(),xtitle,"For_pos_res_energy_" + energies[ie] + target_crystal + "_" + "_profile" + titlevar1[iv] + "_inbinsof_" +  fillvar2_2d[iv] +".png");

	}
    }
}

void for_mor_radius()
{
  string fillvar[] = {"R13","R15","R35"};

  int nvars = sizeof(fillvar)/sizeof(fillvar[0]);

  double bin_low[] = {0.5,0.5,0.5};
  double bin_up[] = {1,1,1};
  int nbin[] = {100,100,100};

  if(nvars != sizeof(bin_low)/sizeof(bin_low[0]) || nvars != sizeof(bin_up)/sizeof(bin_up[0]) || nvars != sizeof(nbin)/sizeof(nbin[0]))
    {cout<<"No of variable and bin information dont match"<<endl;exit(0);}

  for(int ie =0; ie < nenergies; ie++)
    {
      for(int iv =0; iv < nvars; iv++)
	{
	  //if(iv != 0 || ie != 1)
	    //continue;
	  TH1D* hist;
	  hist = new TH1D(("hist_mrx_" + fillvar[iv] + energies[ie]).c_str(),("hist" + fillvar[iv]).c_str(),nbin[iv],bin_low[iv],bin_up[iv]);
	  gethist(energies[ie] + "_"+file_tag+".root",treename,fillvar[iv],"",1,hist);
	  
	  int binmax = hist->GetMaximumBin();
	  double xmean = hist->GetXaxis()->GetBinCenter(binmax);
	  
	  double std = hist->GetStdDev();	  
	  double xmin, xmax;

	  if(iv == 0)
	    {
	      xmin = 0.65;
	      xmax = 0.9;
	    }
	  else if(iv == 1)
	    {
	      xmin = 0.6;
	      xmax = 0.9;
	    }
	  else	  if(iv == 2)
	    {
	      xmin = 0.85;
	      xmax = 1.0;
	    }

	  xmin = max(0.00,xmean - 2.10*std);
	  xmax = min(1.00,xmean + 1.20*std);
	  
	  if(iv == 1 && ie == 0)
	    {
	      xmin = 0.45;
	      xmax = 0.78;
	    }
	  //xmin = 0.8;
	  //xmax = 1.1;
	  int nbins = 60;
	  
	  hist = new TH1D(("hist_mr_" + fillvar[iv] + energies[ie]).c_str(),("hist" + fillvar[iv]).c_str(),nbins,xmin,xmax);
	  gethist(energies[ie] + "_"+file_tag+".root",treename,fillvar[iv],"",1,hist);

	  cout<<endl<<endl<<xmean<<" "<<xmin<<" "<<xmax<<" "<<std<<endl<<endl;

	  double * r_values = get_resolution_value_crystalball_gauss(hist,energies[ie] + " GeV Moliere raidus",xmin,xmax);
	      
	}
    }

}


void perform_diff_event_check()
{
  bool make_plots = false;
  bool make_lookuptable = true;

  gSystem->Exec("mkdir ./plots/desync_plots/");
  ofstream outputFile1,outputFile2;

  if(make_lookuptable)
    {
      outputFile1.open("plots/desync_lookup_table_mcp1.txt",std::ofstream::trunc);
      outputFile1 << "Run Number" <<setw(20)<<"Spill Number"<<setw(20)<<"Event numbers"<<setw(20)<<"offset"<< endl;

      outputFile2.open("plots/desync_lookup_table_mcp2.txt",std::ofstream::trunc);
      outputFile2 << "Run Number" <<setw(20)<<"Spill Number"<<setw(20)<<"Event numbers"<<setw(20)<<"offset"<< endl;
    }
  
  for(int ie =0; ie < nenergies; ie++)
    {
      
      gSystem->Exec(("mkdir ./plots/desync_plots/" + energies[ie]).c_str());
      UInt_t          run;
      UInt_t          spill;
      UInt_t          event;
      UInt_t          trg;
      Int_t           PED;
      Int_t           laser;
      Int_t           phys;
      Int_t           nhits;
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
      Float_t         gain[28];
      Float_t         pedestal[28];
      Float_t         b_slope[28];
      Float_t         b_rms[28];
      Float_t         time[56];
      Float_t         time_error[56];
      Float_t         time_slope[56];
      Float_t         period[28];
      Float_t         maximum[28];
      Float_t         time_maximum[28];
      Float_t         amp_max[28];
      Float_t         time_max[28];
      Float_t         fit_ampl[28];
      Float_t         fit_time[28];
      Float_t         fit_terr[28];
      Float_t         fit_period[28];
      Float_t         amp_max_C2;
      Float_t         E1x1;
      Float_t         E3x3;
      Float_t         E5x5;
      Float_t         R13;
      Float_t         R15;
      Float_t         R35;
      Float_t         dt_C2_MCP1;
      Float_t         dt_C2_MCP2;
      Float_t         dt_MCP1_MCP2;
      Float_t         dt_C2_C1;
      Float_t         dt_C2_C3;
      Float_t         dt_C2_B2;
      Float_t         dt_C2_D2;
      Float_t         dt_C2_B3;
      Float_t         dt_C2_D3;
      Float_t         dt_C2_B1;
      Float_t         dt_C2_D1;
      Float_t         t_C2_C1;
      Float_t         t_C2_C3;
      Float_t         gain_C2;
      Float_t         gain_MCP1;
      Float_t         gain_MCP2;
      Float_t         dt_C2_C1_phase_correct;
      Float_t         dt_C2_C3_phase_correct;
      Float_t         dt_C2_B1_phase_correct;
      Float_t         dt_C2_B2_phase_correct;
      Float_t         dt_C2_B3_phase_correct;
      Float_t         dt_C2_D1_phase_correct;
      Float_t         dt_C2_D2_phase_correct;
      Float_t         dt_C2_D3_phase_correct;
      Int_t           n_h1X;
      Float_t         h1X[10];   //[n_h1X]
      Int_t           n_h2X;
      Float_t         h2X[10];   //[n_h2X]
      Int_t           n_h1Y;
      Float_t         h1Y[11];   //[n_h1Y]
      Int_t           n_h2Y;
      Float_t         h2Y[9];   //[n_h2Y]
      Int_t           n_tracks;
      Float_t         trackX[5];   //[n_tracks]
      Float_t         trackY[5];   //[n_tracks]
      Float_t         amp_max_crystal[25];
      Float_t         posX_crystal[25];
      Float_t         posY_crystal[25];
      Float_t         dt_C3_MCP1;
      Float_t         dt_C3_MCP2;
      Float_t         dt_C3_C2_phase_correct;
      Float_t         dt_C3_C4_phase_correct;
      Float_t         dt_C3_B4_phase_correct;
      Float_t         dt_C3_B2_phase_correct;
      Float_t         dt_C3_B3_phase_correct;
      Float_t         dt_C3_D4_phase_correct;
      Float_t         dt_C3_D2_phase_correct;
      Float_t         dt_C3_D3_phase_correct;
	
      TFile* fileIn;
      if(make_lookuptable) fileIn= new TFile((energies[ie] + "_C32x2_digisync_VFE_fix_test_files.root").c_str(), "read");
      //fileIn = new TFile((energies[ie] + "_C32x2_digisync_VFE_fix_files.root").c_str(), "read");
      //fileIn = new TFile((energies[ie] + "_C32x2_VFEfixe_syncfix_v2.root").c_str(), "read");
      if(make_plots) fileIn = new TFile((energies[ie] + "_C32x2_digisync_rearranged_events.root").c_str(), "read");
      
      TTree *Tin = (TTree*)fileIn->Get("new_tree");

      Tin->SetBranchAddress("run", &run);
      Tin->SetBranchAddress("spill", &spill);
      Tin->SetBranchAddress("event", &event);
      Tin->SetBranchAddress("trg", &trg);
      Tin->SetBranchAddress("PED", &PED);
      Tin->SetBranchAddress("laser", &laser);
      Tin->SetBranchAddress("phys", &phys);
      Tin->SetBranchAddress("nhits", &nhits);
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
      Tin->SetBranchAddress("b_slope", b_slope);
      Tin->SetBranchAddress("b_rms", b_rms);
      Tin->SetBranchAddress("time", time);
      Tin->SetBranchAddress("time_error", time_error);
      Tin->SetBranchAddress("time_slope", time_slope);
      Tin->SetBranchAddress("period", period);
      Tin->SetBranchAddress("maximum", maximum);
      Tin->SetBranchAddress("time_maximum", time_maximum);
      Tin->SetBranchAddress("amp_max", amp_max);
      Tin->SetBranchAddress("time_max", time_max);
      Tin->SetBranchAddress("fit_ampl", fit_ampl);
      Tin->SetBranchAddress("fit_time", fit_time);
      Tin->SetBranchAddress("fit_terr", fit_terr);
      Tin->SetBranchAddress("fit_period", fit_period);
      Tin->SetBranchAddress("amp_max_C2", &amp_max_C2);
      Tin->SetBranchAddress("E1x1", &E1x1);
      Tin->SetBranchAddress("E3x3", &E3x3);
      Tin->SetBranchAddress("E5x5", &E5x5);
      Tin->SetBranchAddress("R13", &R13);
      Tin->SetBranchAddress("R15", &R15);
      Tin->SetBranchAddress("R35", &R35);
      Tin->SetBranchAddress("dt_C2_MCP1", &dt_C2_MCP1);
      Tin->SetBranchAddress("dt_C2_MCP2", &dt_C2_MCP2);
      Tin->SetBranchAddress("dt_MCP1_MCP2", &dt_MCP1_MCP2);
      Tin->SetBranchAddress("dt_C2_C1", &dt_C2_C1);
      Tin->SetBranchAddress("dt_C2_C3", &dt_C2_C3);
      Tin->SetBranchAddress("dt_C2_B2", &dt_C2_B2);
      Tin->SetBranchAddress("dt_C2_D2", &dt_C2_D2);
      Tin->SetBranchAddress("dt_C2_B3", &dt_C2_B3);
      Tin->SetBranchAddress("dt_C2_D3", &dt_C2_D3);
      Tin->SetBranchAddress("dt_C2_B1", &dt_C2_B1);
      Tin->SetBranchAddress("dt_C2_D1", &dt_C2_D1);
      Tin->SetBranchAddress("t_C2_C1", &t_C2_C1);
      Tin->SetBranchAddress("t_C2_C3", &t_C2_C3);
      Tin->SetBranchAddress("gain_C2", &gain_C2);
      Tin->SetBranchAddress("gain_MCP1", &gain_MCP1);
      Tin->SetBranchAddress("gain_MCP2", &gain_MCP2);
      Tin->SetBranchAddress("dt_C2_C1_phase_correct", &dt_C2_C1_phase_correct);
      Tin->SetBranchAddress("dt_C2_C3_phase_correct", &dt_C2_C3_phase_correct);
      Tin->SetBranchAddress("dt_C2_B1_phase_correct", &dt_C2_B1_phase_correct);
      Tin->SetBranchAddress("dt_C2_B2_phase_correct", &dt_C2_B2_phase_correct);
      Tin->SetBranchAddress("dt_C2_B3_phase_correct", &dt_C2_B3_phase_correct);
      Tin->SetBranchAddress("dt_C2_D1_phase_correct", &dt_C2_D1_phase_correct);
      Tin->SetBranchAddress("dt_C2_D2_phase_correct", &dt_C2_D2_phase_correct);
      Tin->SetBranchAddress("dt_C2_D3_phase_correct", &dt_C2_D3_phase_correct);
      Tin->SetBranchAddress("n_h1X", &n_h1X);
      Tin->SetBranchAddress("h1X", h1X);
      Tin->SetBranchAddress("n_h2X", &n_h2X);
      Tin->SetBranchAddress("h2X", h2X);
      Tin->SetBranchAddress("n_h1Y", &n_h1Y);
      Tin->SetBranchAddress("h1Y", h1Y);
      Tin->SetBranchAddress("n_h2Y", &n_h2Y);
      Tin->SetBranchAddress("h2Y", h2Y);
      Tin->SetBranchAddress("n_tracks", &n_tracks);
      Tin->SetBranchAddress("trackX", trackX);
      Tin->SetBranchAddress("trackY", trackY);
      Tin->SetBranchAddress("amp_max_crystal", amp_max_crystal);
      Tin->SetBranchAddress("posX_crystal", posX_crystal);
      Tin->SetBranchAddress("posY_crystal", posY_crystal);
      Tin->SetBranchAddress("dt_C3_MCP1", &dt_C3_MCP1);
      Tin->SetBranchAddress("dt_C3_MCP2", &dt_C3_MCP2);
      Tin->SetBranchAddress("dt_C3_C2_phase_correct", &dt_C3_C2_phase_correct);
      Tin->SetBranchAddress("dt_C3_C4_phase_correct", &dt_C3_C4_phase_correct);
      Tin->SetBranchAddress("dt_C3_B4_phase_correct", &dt_C3_B4_phase_correct);
      Tin->SetBranchAddress("dt_C3_B2_phase_correct", &dt_C3_B2_phase_correct);
      Tin->SetBranchAddress("dt_C3_B3_phase_correct", &dt_C3_B3_phase_correct);
      Tin->SetBranchAddress("dt_C3_D4_phase_correct", &dt_C3_D4_phase_correct);
      Tin->SetBranchAddress("dt_C3_D2_phase_correct", &dt_C3_D2_phase_correct);
      Tin->SetBranchAddress("dt_C3_D3_phase_correct", &dt_C3_D3_phase_correct);      
      int nentries = Tin->GetEntries();

      if(make_plots)
	{

	  for(int ievt_crystal = 5; ievt_crystal < 6; ievt_crystal++)
	    {
	      for(int ievt_MCP = 0; ievt_MCP < 11; ievt_MCP++)
		{

		  int diff_evt = ievt_crystal - ievt_MCP;
		  string diff_evt_str;
		  if(diff_evt > 0)	    diff_evt_str = "#plus" + to_string(abs(diff_evt));
		  else if(diff_evt < 0)	    diff_evt_str = "#minus" + to_string(abs(diff_evt));
		  else 	 diff_evt_str = "";

		  TH1D *hist_C3_MCP1[1];
		  hist_C3_MCP1[0] = new TH1D("hist_dt_C3_MCP1","dt_C3_MCP1",50,0,6238);
		  hist_C3_MCP1[0]->GetXaxis()->SetTitle(("time difference between C3(N) and MCP1(N" + diff_evt_str+ ")").c_str());  
		  hist_C3_MCP1[0]->GetYaxis()->SetTitle();  

		  TH1D *hist_C3_MCP2[1];
		  hist_C3_MCP2[0] = new TH1D("hist_dt_C3_MCP2","dt_C3_MCP2",50,0,6238);
		  hist_C3_MCP2[0]->GetXaxis()->SetTitle(("time difference between C3(N) and MCP2(N" +diff_evt_str + ")").c_str());  
		  hist_C3_MCP2[0]->GetYaxis()->SetTitle();  


		  TH1D *hist_C3_MCP1_eventwise[1];
		  hist_C3_MCP1_eventwise[0] = new TH1D("hist_dt_C3_MCP1","dt_C3_MCP1",50,0,6238);
		  hist_C3_MCP1_eventwise[0]->GetXaxis()->SetTitle(("time difference between C3(N) and MCP1(N" + diff_evt_str+ ")").c_str());  
		  hist_C3_MCP1_eventwise[0]->GetYaxis()->SetTitle();  

		  TH1D *hist_C3_MCP2_eventwise[1];
		  hist_C3_MCP2_eventwise[0] = new TH1D("hist_dt_C3_MCP2","dt_C3_MCP2",50,0,6238);
		  hist_C3_MCP2_eventwise[0]->GetXaxis()->SetTitle(("time difference between C3(N) and MCP2(N" +diff_evt_str + ")").c_str());  
		  hist_C3_MCP2_eventwise[0]->GetYaxis()->SetTitle();  

		  TH2D *h3_spill_event_dt_C3_MCP1 = new TH2D("h3_spill_event_dt_C3_MCP1","h3_spill_event_dt_C3_MCP1",101,-0.5,100.05,4501,-0.5,4500.5);

		  h3_spill_event_dt_C3_MCP1->GetXaxis()->SetTitle("Spill");
		  h3_spill_event_dt_C3_MCP1->GetYaxis()->SetTitle("Event");
		  h3_spill_event_dt_C3_MCP1->GetZaxis()->SetTitle(("time difference between C3(N) and MCP1(N" + diff_evt_str+ ")").c_str());

		  TH2D * h3_spill_event_dt_C3_MCP2 = new TH2D("h3_spill_event_dt_C3_MCP2","h3_spill_event_dt_C3_MCP2",101,-0.5,100.5,4501,-0.5,4500.5);

		  h3_spill_event_dt_C3_MCP2->GetXaxis()->SetTitle("Spill");
		  h3_spill_event_dt_C3_MCP2->GetYaxis()->SetTitle("Event");
		  h3_spill_event_dt_C3_MCP2->GetZaxis()->SetTitle(("time difference between C3(N) and MCP2(N" + diff_evt_str + ")").c_str());

		  const int nspills_tocheck = 100;
		  TH2D *h3_spillwise_event_vs_dt_C3_MCP1[nspills_tocheck];
		  TH2D *h3_spillwise_event_vs_dt_C3_MCP2[nspills_tocheck];
	  
		  for(int ispill = 0; ispill < nspills_tocheck; ispill ++)
		    {

		      h3_spillwise_event_vs_dt_C3_MCP1[ispill] = new TH2D(("h3_spill_event_dt_C3_MCP1_"+ to_string(ispill)).c_str(),"h3_spill_event_dt_C3_MCP1",40,-0.5,4500.5,40,0,6500);

		      h3_spillwise_event_vs_dt_C3_MCP1[ispill]->GetXaxis()->SetTitle("Event number");
		      //h3_spillwise_event_vs_dt_C3_MCP1[ispill]->GetYaxis()->SetTitle("Event");
		      h3_spillwise_event_vs_dt_C3_MCP1[ispill]->GetYaxis()->SetTitle(("time difference between C3(N) and MCP1(N" + diff_evt_str+ ")").c_str());

		      h3_spillwise_event_vs_dt_C3_MCP2[ispill] = new TH2D(("h3_spill_event_dt_C3_MCP2_"+ to_string(ispill)).c_str(),"h3_spill_event_dt_C3_MCP2",40,-0.5,4500.5,40,0,6500);
	      
		      h3_spillwise_event_vs_dt_C3_MCP2[ispill]->GetXaxis()->SetTitle("Event number");
		      //h3_spillwise_event_vs_dt_C3_MCP2[ispill]->GetYaxis()->SetTitle("Event");
		      h3_spillwise_event_vs_dt_C3_MCP2[ispill]->GetYaxis()->SetTitle(("time difference between C3(N) and MCP2(N" + diff_evt_str + ")").c_str());
		
		    }

		  for (int iev=ievt_crystal; iev<nentries-max(0,diff_evt); iev++)
		    {
		      Tin->GetEntry(iev);
	      
		      float time_C3 = fit_time[C3];

		      if(iev+diff_evt < 0 || iev+diff_evt > nentries - 1) continue;
	      
		      Tin->GetEntry(iev+diff_evt);

		      //if(spill == 62 || spill == 29) continue;
		      float time_MCP1 = fit_time[MCP1];
		      float time_MCP2 = fit_time[MCP2];
		      float time_CLK = fit_time[CLK];
      
		      float dt_C3_MCP1_test = (time_C3 - time_MCP1 + time_CLK - (int((time_C3 - time_MCP1 + time_CLK)/6.238)*6.238))*1000;
		      float dt_C3_MCP2_test = (time_C3 - time_MCP2 + time_CLK - (int((time_C3 - time_MCP2 + time_CLK)/6.238)*6.238))*1000;

		      if(diff_evt == 0)
			{
			  dt_C3_MCP1_test = dt_C3_MCP1*1000;
			  dt_C3_MCP2_test = dt_C3_MCP2*1000;
			}
		      
		      hist_C3_MCP1[0]->Fill(dt_C3_MCP1_test);
		      hist_C3_MCP2[0]->Fill(dt_C3_MCP2_test);
		      h3_spill_event_dt_C3_MCP1->Fill(spill,event,dt_C3_MCP1_test);
		      h3_spill_event_dt_C3_MCP2->Fill(spill,event,dt_C3_MCP2_test);
		  
		      if(spill == 1 && diff_evt < 2 && diff_evt > -1)cout<< "1 - Finalized = "<<energies[ie] <<" GeV"<<setw(20)<<run <<setw(20)<< spill <<setw(20)<<event<<" "<<diff_evt<< "    "<<dt_C3_MCP1_test<<"  "<<time_C3<<"  "<<time_MCP1<<"  "<<time_CLK<<endl;

		      if(amp_max[MCP1] > 200) h3_spillwise_event_vs_dt_C3_MCP1[0]->Fill(event,dt_C3_MCP1_test);
		      if(amp_max[MCP2] > 200) h3_spillwise_event_vs_dt_C3_MCP2[0]->Fill(event,dt_C3_MCP2_test);

		      for(int ispill = 1; ispill < nspills_tocheck; ispill ++)
			{
			  if(spill == ispill && amp_max[MCP1] > 200) h3_spillwise_event_vs_dt_C3_MCP1[ispill]->Fill(event,dt_C3_MCP1_test);
			  if(spill == ispill && amp_max[MCP2] > 200) h3_spillwise_event_vs_dt_C3_MCP2[ispill]->Fill(event,dt_C3_MCP2_test);
			}
		    }
		  string str[1] = {""};

		  plot_1dhists(1,hist_C3_MCP1,str,"desync_plots/" + energies[ie] + "/" + energies[ie] + "_energy_dt_C3_MCP1_C3_N" + to_string(ievt_crystal) + "_MCP1_CLK_N" + to_string(ievt_MCP) + ".png");
		  plot_1dhists(1,hist_C3_MCP2,str,"desync_plots/" + energies[ie] + "/" + energies[ie] + "_energy_dt_C3_MCP1_C3_N" + to_string(ievt_crystal) + "_MCP2_CLK_N" + to_string(ievt_MCP) + ".png");

		  plotTH2D(h3_spill_event_dt_C3_MCP1,"desync_plots/" + energies[ie] + "/" + energies[ie] + "_energy_dt_C3_MCP1_C3_N" + to_string(ievt_crystal) + "_MCP1_CLK_N" + to_string(ievt_MCP) + "_spill_event");
		  plotTH2D(h3_spill_event_dt_C3_MCP2,"desync_plots/" + energies[ie] + "/" + energies[ie] + "_energy_dt_C3_MCP2_C3_N" + to_string(ievt_crystal) + "_MCP2_CLK_N" + to_string(ievt_MCP) + "_spill_event");

		  for(int ispill = 0; ispill < nspills_tocheck; ispill ++)
		    {
		      if(ispill == 0)
			{
			  plotTH2D(h3_spillwise_event_vs_dt_C3_MCP1[ispill],"desync_plots/" + energies[ie] + "/" + energies[ie] + "_energy_dt_C3_MCP1_C3_N" + to_string(ievt_crystal) + "_MCP1_CLK_N" + to_string(ievt_MCP) + "_vs_event_forallspill");
			  plotTH2D(h3_spillwise_event_vs_dt_C3_MCP2[ispill],"desync_plots/" + energies[ie] + "/" + energies[ie] + "_energy_dt_C3_MCP2_C3_N" + to_string(ievt_crystal) + "_MCP2_CLK_N" + to_string(ievt_MCP) + "_vs_event_forallspill");
			}
		      else
			{
			  if(h3_spillwise_event_vs_dt_C3_MCP1[ispill]->Integral() < 1) continue;
		      
			  gSystem->Exec(("mkdir ./plots/desync_plots/" + energies[ie] + "/" + to_string(ispill) ).c_str());
			  plotTH2D(h3_spillwise_event_vs_dt_C3_MCP1[ispill],"desync_plots/" + energies[ie] + "/" + to_string(ispill) + "/" + energies[ie] + "_energy_dt_C3_MCP1_C3_N" + to_string(ievt_crystal) + "_MCP1_CLK_N" + to_string(ievt_MCP) + "_vs_event_forspill" + to_string(ispill));
			  plotTH2D(h3_spillwise_event_vs_dt_C3_MCP2[ispill],"desync_plots/" + energies[ie] + "/" + to_string(ispill) + "/" + energies[ie] + "_energy_dt_C3_MCP2_C3_N" + to_string(ievt_crystal) + "_MCP2_CLK_N" + to_string(ievt_MCP) + "_vs_event_forspill" + to_string(ispill));
			}
		    }
		}


	    }
	}


      //////////////// LOOKUP table calculation

    
      if(make_lookuptable)
	{
	  cout<<"Start lookup table calculation"<<endl;
	  TFile* fileOut = new TFile((energies[ie] + "_C32x2_digisync_rearranged_events_v2.root").c_str(), "recreate");  
	  TTree* Tout = new TTree("new_tree", "new_tree");

	  Tout->Branch("run", &run, "run/i");
	  Tout->Branch("spill", &spill, "spill/i");
	  Tout->Branch("event", &event, "event/i");

	  Tout->Branch("trg", &trg, "trg/i");
	  Tout->Branch("PED", &PED, "PED/I");
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
	  Tout->Branch("b_slope", b_slope, "b_slope[28]/F");
	  Tout->Branch("b_rms", b_rms, "b_rms[28]/F");
	  Tout->Branch("time", time, "time[56]/F");
	  Tout->Branch("time_error", time_error, "time_error[56]/F");
	  Tout->Branch("time_slope", time_slope, "time_slope[56]/F");
	  Tout->Branch("period", period, "period[28]/F");
	  Tout->Branch("maximum", maximum, "maximum[28]/F");
	  Tout->Branch("time_maximum", time_maximum, "time_maximum[28]/F");
	  Tout->Branch("amp_max", amp_max, "amp_max[28]/F");
	  Tout->Branch("time_max", time_max, "time_max[28]/F");
	  Tout->Branch("fit_ampl", fit_ampl, "fit_ampl[28]/F");
	  Tout->Branch("fit_time", fit_time, "fit_time[28]/F");
	  Tout->Branch("fit_terr", fit_terr, "fit_terr[28]/F");
	  Tout->Branch("fit_period", fit_period, "fit_period[28]/F");
  
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

	  Tout->Branch("n_h1X", &n_h1X,"n_h1X/I");   
	  Tout->Branch("h1X", h1X,"h1X[n_h1X]/F");   

	  Tout->Branch("n_h2X", &n_h2X,"n_h2X/I");   
	  Tout->Branch("h2X", h2X,"h2X[n_h2X]/F");   

	  Tout->Branch("n_h1Y", &n_h1Y,"n_h1Y/I");   
	  Tout->Branch("h1Y", h1Y,"h1Y[n_h1Y]/F");   

	  Tout->Branch("n_h2Y", &n_h2Y,"n_h2Y/I");   
	  Tout->Branch("h2Y", h2Y,"h2Y[n_h2Y]/F");   

	  Tout->Branch("n_tracks", &n_tracks,"n_tracks/I");   
	  Tout->Branch("trackX", trackX,"trackX[n_tracks]/F");   
	  Tout->Branch("trackY", trackY,"trackY[n_tracks]/F");   
    
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


	  int max_ndiff_events_mcp = 20;       ////////// maximum number of offset to check
	  int number_events_tosync = 50;      ////////// maximum number of event to check in the one iteration of loop 
	  
	  int max_event_tocheck = number_events_tosync;
	  
	  for (int iev=0; iev<nentries/*- ndiff_mcp1*/; iev+=max(1,max_event_tocheck))
	    {
	      int ndiff_mcp1 = 0;
	      int ndiff_mcp2 = 0;

	      fileIn->cd();

	      float std = 100000;    /////// for calculating paramter to check the syncness

	      max_event_tocheck =  number_events_tosync;    /////////// number of events to check in this iteration of loop
	      
	      int start_eventno_tocheck;     /////////// store the event number of the first event checked in this iteration of loop
	      int end_eventno_tocheck;   /////////// store the event number of the last event checked in this iteration of loop

	      float expected_std = 1000;
	  
	      bool spill_changed =false;        /////////   to check are we at the end of spill (means number of events left in this spill is less than  variable "number_events_tosync")
	      bool end_of_events = false;       /////////   to check are we at the end of all runs (means number of events left in all runs is less than  variable "number_events_tosync")
	      UInt_t current_spill;             ////////////// to track the what is the current spill number
	      
	      while(std > expected_std)
		{

		  double dt[number_events_tosync];    ////////// to store the time values and calculate the standard devartion from them

		  for(int i=0; i < max_event_tocheck; i++)   //////////////// looping over the events to check in this iteration and caluclating the time differences
		    {
		  
		      Tin->GetEntry(iev+i);         ///////////////  get the events for crystal without shifting
		      
		      float time_C3 = fit_time[C3];
		      
		      if(i==0)current_spill = spill;   //////// store what is the curent spill number
		      
		      if(iev + i + ndiff_mcp1 > nentries-1)   //////// check if the shifted event is going out of range or not
			{
			  max_event_tocheck = i;    //////////////  set the maxium number to check in this iteration to be upto before than this event
			  end_of_events = true;
			}
			  
		      Tin->GetEntry(iev + i + ndiff_mcp1);   /////// get the shifted event
		      
		      if(current_spill != spill)     ////// check the spill is same as that of in the begining
			{
			  max_event_tocheck =  i;   //////////////  set the maxium number to check in this iteration to be upto before than this event
			  spill_changed = true;
			  //exit(0);
			  break;
			}
		      float time_MCP1 = fit_time[MCP1];
		      float time_CLK = fit_time[CLK];
      
		      float dt_C3_MCP1_test = (time_C3 - time_MCP1 + time_CLK - (int((time_C3 - time_MCP1 + time_CLK)/6.238)*6.238))*1000;
		      dt[i] = dt_C3_MCP1_test;
		    }

		  if(max_event_tocheck == 0)  break;
		  
		  Tin->GetEntry(iev+max_event_tocheck-1);  ////// get the event number of last checked event
		  end_eventno_tocheck = event;
		  
		  Tin->GetEntry(iev);   ////// get the event number of first checked event
		  start_eventno_tocheck =  event;
      
		  //std = calculate_std_from_array(dt,max_event_tocheck);
		  std = calculate_ithquantile_from_array(dt,max_event_tocheck,0.5);   ///////// calculate the inter quantile range
	 
		  //if(max_event_tocheck < number_events_tosync) break;
	      
		  //cout<< energies[ie] <<" GeV"<<setw(20)<<run <<setw(20)<< spill <<setw(20)<<iev<<" -"<<start_eventno_tocheck << " - " << end_eventno_tocheck<<setw(20)<< ndiff_mcp1 << "    "<<std<<endl;

		  if(ndiff_mcp1 > max_ndiff_events_mcp)   
		    std = 10;
		  if(std > expected_std)
		    ndiff_mcp1++;
		}
	  
	
	      outputFile1 << run <<setw(20)<< spill <<setw(20)<<iev<<" -"<<start_eventno_tocheck << " - " << end_eventno_tocheck<<setw(20)<< ndiff_mcp1 << endl;

	      if(ndiff_mcp1 < max_ndiff_events_mcp)
		{
		  TH1D *hist_C3_MCP2_diffeventwise[1];   /////// to store the distribution of checked time difference
		  hist_C3_MCP2_diffeventwise[0] = new TH1D("hist_dt_C3_MCP2","dt_C3_MCP2",50,0,6238);
		  hist_C3_MCP2_diffeventwise[0]->GetXaxis()->SetTitle("time difference between C3 and MCP1");  
		  hist_C3_MCP2_diffeventwise[0]->GetYaxis()->SetTitle();
	      
		  for(int i=0; i < max_event_tocheck; i++)  //////// start storing the events
		    {
		      fileIn->cd();	  

		      Tin->GetEntry(iev+i);
			
		      float time_C3 = fit_time[C3];
		      float time_C2 = fit_time[C2];

		      Tin->GetEntry(iev + i + ndiff_mcp1);

		      float time_MCP1 = fit_time[MCP1];
		      float time_MCP2 = fit_time[MCP2];
		      float time_CLK = fit_time[CLK];

		      float amp_max1 = amp_max[MCP1];
		      float amp_max2 = amp_max[MCP2];
           
		      float pedestal_MCP1 = pedestal[MCP1];
		      float b_slope_MCP1 = b_slope[MCP1];
		      float b_rms_MCP1 = b_rms[MCP1];
		      float time_error_MCP1 = time_error[MCP1];
		      float time_slope_MCP1 = time_slope[MCP1];
		      float period_MCP1 = period[MCP1];
		      float maximum_MCP1 = maximum[MCP1];
		      float time_maximum_MCP1 = time_maximum[MCP1];
		      float time_max_MCP1 = time_max[MCP1];
		      float fit_ampl_MCP1 = fit_ampl[MCP1];
		      float fit_time_MCP1 = fit_time[MCP1];
		      float fit_terr_MCP1 = fit_terr[MCP1];
		      float fit_period_MCP1 = fit_period[MCP1];
		      float gain_MCP1_out = gain_MCP1;
		      float gain_MCP2_out = gain_MCP2;
	 
		      float pedestal_MCP2 = pedestal[MCP2];
		      float b_slope_MCP2 = b_slope[MCP2];
		      float b_rms_MCP2 = b_rms[MCP2];
		      float time_error_MCP2 = time_error[MCP2];
		      float time_slope_MCP2 = time_slope[MCP2];
		      float period_MCP2 = period[MCP2];
		      float maximum_MCP2 = maximum[MCP2];
		      float time_maximum_MCP2 = time_maximum[MCP2];
		      float time_max_MCP2 = time_max[MCP2];
		      float fit_ampl_MCP2 = fit_ampl[MCP2];
		      float fit_time_MCP2 = fit_time[MCP2];
		      float fit_terr_MCP2 = fit_terr[MCP2];
		      float fit_period_MCP2 = fit_period[MCP2];
		      
		      float pedestal_CLK = pedestal[CLK];
		      float b_slope_CLK = b_slope[CLK];
		      float b_rms_CLK = b_rms[CLK];
		      float time_error_CLK = time_error[CLK];
		      float time_slope_CLK = time_slope[CLK];
		      float period_CLK = period[CLK];
		      float maximum_CLK = maximum[CLK];
		      float time_maximum_CLK = time_maximum[CLK];
		      float amp_max_CLK = amp_max[CLK];
		      float time_max_CLK = time_max[CLK];
		      float fit_ampl_CLK = fit_ampl[CLK];
		      float fit_time_CLK = fit_time[CLK];
		      float fit_terr_CLK = fit_terr[CLK];
		      float fit_period_CLK = fit_period[CLK];
		      float gain_CLK_out = gain[CLK];

		      float dt_C3_MCP1_out = (time_C3 - time_MCP1 + time_CLK - (int((time_C3 - time_MCP1 + time_CLK)/6.238)*6.238));
		      float dt_C3_MCP2_out = (time_C3 - time_MCP2 + time_CLK - (int((time_C3 - time_MCP2 + time_CLK)/6.238)*6.238));
		      if(spill == 1)cout<< "1 - Finalized = "<<energies[ie] <<" GeV"<<setw(20)<<run <<setw(20)<< spill <<setw(20)<<event<<" "<<ndiff_mcp1<< "    "<<dt_C3_MCP1_out<<"  "<<time_C3<<"  "<<time_MCP1<<"  "<<time_CLK<<" "<<amp_max[MCP1]<<" "<<amp_max[MCP2]<<endl;


		      Tin->GetEntry(iev+i);

		      //cout<< "Finalized = "<<energies[ie] <<" GeV"<<setw(20)<<run <<setw(20)<< spill <<setw(20)<<start_eventno_tocheck << " - " << start_eventno_tocheck + max_event_tocheck<<setw(20)<<n diff_mcp1 << "    "<<dt_C3_MCP1_out<<endl;
		      dt_C3_MCP1 = (time_C3 - time_MCP1 + time_CLK - (int((time_C3 - time_MCP1 + time_CLK)/6.238)*6.238));
		      dt_C3_MCP2 = (time_C3 - time_MCP2 + time_CLK - (int((time_C3 - time_MCP2 + time_CLK)/6.238)*6.238));
		      dt_C2_MCP1 = (time_C2 - time_MCP1 + time_CLK - (int((time_C2 - time_MCP1 + time_CLK)/6.238)*6.238));
		      dt_C2_MCP2 = (time_C2 - time_MCP2 + time_CLK - (int((time_C2 - time_MCP2 + time_CLK)/6.238)*6.238));
		      amp_max[MCP1] = amp_max1;
		      amp_max[MCP2] = amp_max2;

		      pedestal[MCP1] = pedestal_MCP1;
		      b_slope[MCP1] = b_slope_MCP1;
		      b_rms[MCP1] = b_rms_MCP1;
		      time[MCP1] =  time_MCP1;
		      time_error[MCP1] =  time_error_MCP1;
		      time_slope[MCP1] = time_slope_MCP1;
		      period[MCP1] =  period_MCP1;
		      maximum[MCP1] =  maximum_MCP1;
		      time_maximum[MCP1] =  time_maximum_MCP1;
		      time_max[MCP1] = time_max_MCP1;
		      fit_ampl[MCP1] =  fit_ampl_MCP1;
		      fit_time[MCP1] =  fit_time_MCP1;
		      fit_terr[MCP1] = fit_terr_MCP1;
		      fit_period[MCP1] =  fit_period_MCP1;
		      gain_MCP1 =  gain_MCP1_out;
		      gain_MCP2 = gain_MCP2_out;
	 
		      pedestal[MCP2] =  pedestal_MCP2;
		      b_slope[MCP2] = b_slope_MCP2;
		      b_rms[MCP2] =  b_rms_MCP2;
		      time[MCP2]=  time_MCP2 ;
		      time_error[MCP2] =  time_error_MCP2;
		      time_slope[MCP2] =  time_slope_MCP2;
		      period[MCP2] =  period_MCP2;
		      maximum[MCP2] =  maximum_MCP2;
		      time_maximum[MCP2] =  time_maximum_MCP2;
		      time_max[MCP2] =  time_max_MCP2;
		      fit_ampl[MCP2] =  fit_ampl_MCP2;
		      fit_time[MCP2] =  fit_time_MCP2;
		      fit_terr[MCP2] =  fit_terr_MCP2;
		      fit_period[MCP2] =  fit_period_MCP2;

		      pedestal[CLK] =  pedestal_CLK;
		      b_slope[CLK] =  b_slope_CLK ;
		      b_rms[CLK]=  b_rms_CLK;
		      time[CLK] =  time_CLK;
		      time_error[CLK] = time_error_CLK;
		      time_slope[CLK] = time_slope_CLK ;
		      period[CLK] =  period_CLK;
		      maximum[CLK] =  maximum_CLK ;
		      time_maximum[CLK] =  time_maximum_CLK ;
		      amp_max[CLK]  =  amp_max_CLK;
		      time_max[CLK] = time_max_CLK ;
		      fit_ampl[CLK] =  fit_ampl_CLK ;
		      fit_time[CLK] =  fit_time_CLK ;
		      fit_terr[CLK] =  fit_terr_CLK;
		      fit_period[CLK]  =  fit_period_CLK ;
		      gain[CLK] =  gain_CLK_out;
		      
		      cout<< "2 - Finalized = "<<energies[ie] <<" GeV"<<setw(20)<<run <<setw(20)<< spill <<setw(20)<<event<<" "<<ndiff_mcp1<< "    "<<dt_C3_MCP1_out<<"  "<<time_C3<<"  "<<time_MCP1<<"  "<<time_CLK<<" "<<amp_max[MCP1]<<" "<<amp_max[MCP2]<<endl;

		      hist_C3_MCP2_diffeventwise[0]->Fill(dt_C3_MCP1_out*1000);
		      fileOut->cd();	  
		      Tout->Fill();
		    }
		  string str[1] = {""};
		  //plot_1dhists(1,hist_C3_MCP2_diffeventwise,str,("desync_plots/" + energies[ie] + "/" + to_string(spill) +"/energy_dt_C3_MCP1_C3" + to_string(start_eventno_tocheck) + ".png").c_str());		
		}
	      
	      if(spill_changed)
		{
		  ndiff_mcp1 = 0;
		  spill_changed = false;
		}
	      if(end_of_events)
		break;

	    }

      
	  fileOut->Write();
	  fileOut->Close();
	}
    
    }
  outputFile1.close();
  outputFile2.close();
  
}

void calculate_clock_period()
{

  //string filename = "_C32x2_digisync_VFE_fix_test_files.root";
  string filename = "_C32x2_digisync_rearranged_events_v2.root";
  string fillvar1_2d[] = {"time_maximum[CLK]"};
  string fillvar2_2d[] = {"maximum[CLK]"};

  string cuts[] ={""};
  
  int nvars1_2d = sizeof(fillvar1_2d)/sizeof(fillvar1_2d[0]);
  int nvars2_2d = sizeof(fillvar2_2d)/sizeof(fillvar2_2d[0]);

  double bin_low2[] = {1250};
  double bin_up2[] = {1450};
  int nbin2[] = {10000};

  double bin_low1[] = {30.5};
  double bin_up1[] = {175.5};
  int nbin1[] = {145};

  if(nvars1_2d != nvars2_2d || nvars1_2d != sizeof(bin_up1)/sizeof(bin_up1[0]) || nvars1_2d != sizeof(nbin1)/sizeof(nbin1[0]) ||  nvars1_2d != sizeof(nbin2)/sizeof(nbin2[0]) ||  nvars1_2d != sizeof(bin_up2)/sizeof(bin_up2[0]))
    {cout<<"No of variable and bin information dont match for 2d histograms"<<endl;exit(0);}

  TH1D * hist_clock_phase = new TH1D("hist_clock_phase","clock_phase",60,4.1,8.1);
  vector<double> highs;
  vector<double> lows;
  for(int ie =0; ie < nenergies; ie++)
    {
      for(int iv =0; iv < nvars1_2d; iv++)
	{
	  TH2D* hist;
	  hist = new TH2D(("hist_nocuts_" + fillvar1_2d[iv] + "_" + fillvar1_2d[iv] + energies[ie]).c_str(),("hist_" + fillvar1_2d[iv] + "_"+ fillvar2_2d[iv]).c_str(),nbin1[iv],bin_low1[iv],bin_up1[iv],nbin2[iv],bin_low2[iv],bin_up2[iv]);
	  hist->GetXaxis()->SetTitle(fillvar1_2d[iv].c_str());  
	  hist->GetYaxis()->SetTitle(fillvar2_2d[iv].c_str());  

	  get2dhist(energies[ie] + filename,treename,fillvar2_2d[iv] + ":" + fillvar1_2d[iv],cuts[iv],1,hist);
	  
	  savehist2d(hist,"energy_" + energies[ie] + "C3_" + "_no_cuts_" + fillvar1_2d[iv] + "_" +  fillvar2_2d[iv],false);
	  TH1D * h = (TH1D*)hist->ProfileX();
	  plot_hist(h,"clock time","energy_" + energies[ie] + "_clock_time_vs_clock_amp.png");
	  for(int ib=1; ib < h->GetNbinsX()+1; ib++)
	    {
	      float bin_prev =h->GetBinContent(max(0,ib-1));
	      float bin_current = h->GetBinContent(ib);
	      float bin_next =h->GetBinContent(min(h->GetNbinsX()+1,ib+1));
	      if(bin_current > bin_prev && bin_current > bin_next)
		highs.push_back(h->GetXaxis()->GetBinCenter(ib));
	      if(bin_current < bin_prev && bin_current < bin_next)
		lows.push_back(h->GetXaxis()->GetBinCenter(ib));
	    }
	}
    }

  for(int i=0;i<highs.size()-1;i++)
    {
      cout<<highs[i+1] - highs[i]<<endl;
      hist_clock_phase->Fill(highs[i+1] - highs[i]);
    }
  for(int i=0;i<lows.size()-1;i++)
    {
      cout<<lows[i+1] - lows[i]<<endl;
      hist_clock_phase->Fill(lows[i+1] - lows[i]);
    }

  plot_hist(hist_clock_phase,"preiod of clock","all_energy_without_fitted_clockphase.png");

  string fillvar[] = {"fit_period[CLK]"};

  int nvars = sizeof(fillvar)/sizeof(fillvar[0]);
  
  double bin_low[] = {6.16};
  double bin_up[] = {6.4};
  int nbin[] = {30};

  if(nvars != sizeof(bin_low)/sizeof(bin_low[0]) || nvars != sizeof
     (bin_up)/sizeof(bin_up[0]) || nvars != sizeof(nbin)/sizeof(nbin[0]))
    {cout<<"No of variable and bin information dont match"<<endl;exit(0);}

  for(int ie =0; ie < nenergies; ie++)
    {
      for(int iv =0; iv < nvars; iv++)
	{
	  TH1D* hist;
	  hist = new TH1D(("hist_nocuts_" + fillvar[iv] + energies[ie]).c_str(),("hist" + fillvar[iv]).c_str(),nbin[iv],bin_low[iv],bin_up[iv]);
	  //gethist(energies[ie] + filename,treename,fillvar[iv],"",1,hist);
	  gethist(energies[ie] + filename,treename,fillvar[iv],"",1,hist);
	  
	  string str = fillvar[iv];
	  plot_hist(hist,str,"energy_" + energies[ie] + "_clock_period.png");

	  double *r_values = get_resolution_value_gauss(hist,"Clock fit period fit",bin_low[iv],bin_up[iv]); //crystalball
	}
    }


}

void functions()
{
  
  //plot_time_resolution();
  
  intercystal_timeresolution();
  
  //plot_energy_resolution();
  
  //basic_plots();
  
  //calculate_intercalritationfactor();
  
  //plot_position_resolution();

  //for_mor_radius();
  
  //perform_diff_event_check();

  //calculate_clock_period();
  
  
  // vector<double> values;
  
  // double *r = CalculateMedian("100_C3D2x2_" + file_tag +".root",treename,"1.0 / sqrt( pow(b_rms[C3]/amp_max[C3],2) + pow(b_rms[C2]/amp_max[C2],2))","1.0 / sqrt( pow(b_rms[C3]/amp_max[C3],2) + pow(b_rms[C2]/amp_max[C2],2)) > 550 && 1.0 / sqrt( pow(b_rms[C3]/amp_max[C3],2) + pow(b_rms[C2]/amp_max[C2],2)) < 1000 && amp_max[C2]>0 && amp_max[C3]> 0 && n_h1X > 0 && n_h1Y > 0 && fabs(amp_max[C2]-amp_max[C3])/max(amp_max[C2],amp_max[C3]) <0.300000&& h1X > -3.000000 && h1X < 4.000000 && h1Y > -2.000000 && h1Y < 4.000000 ",values);
}
