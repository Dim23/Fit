#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TROOT.h>
#include <TF1.h>
#include <TLine.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <Math/PdfFuncMathCore.h>
#include <Math/IntegratorOptions.h>

const int NO = 1;
int Nch0 ;
//Задаем постоянные;
const double sigma = 685, pi = TMath::Pi(), bmax = 18;
int color[10] = {1, 2, 3, 4, 42, 6, 46, 8, 9, 12};
double a1 = -4.125, a2 = 1.51, a3 = -3.29, a4, teta = 1.41, n_knee = 275.5;
int bin_cent[11] = {400, 191, 134, 93, 63, 41, 25, 14, 7, 3, 2};

int Nn;
int current_id;
TFile *d_outfile;

//Уссловная вероятность;
double PnbGaus(double *x, double *par)
{
    double Cb = x[0], n = par[0], teta = par[1], n_knee = par[2], a1 = par[3], a2 = par[4], a3 = par[5];
    double fn = n_knee * exp(a1 * Cb + a2 * pow(Cb, 2) + a3 * pow(Cb, 3)) / teta;
    return ROOT::Math::gamma_pdf(n, fn, teta);
};

double PnbGaus2(double *x, double *par)
{
    double n = x[0], Cb = par[0], teta = par[1], n_knee = par[2], a1 = par[3], a2 = par[4], a3 = par[5];
    double fn = n_knee * exp(a1 * Cb + a2 * pow(Cb, 2) + a3 * pow(Cb, 3)) / teta;
    return ROOT::Math::gamma_pdf(n, fn, teta);
};

double ftPn(double *x, double *par)
{
    // Parameters
    double teta = par[0];
    double n_knee = par[1];
    double a1 = par[2];
    double a2 = par[3];
    double a3 = par[4];

    // Variables
    double n = x[0];
    // Function
    TF1 *f = new TF1("f", PnbGaus, 0, 10, 6);
    f->SetParameters(n, teta, n_knee, a1, a2, a3);
    double func = f->Integral(0, 1);
    return func;
}

void Start(const char *fileadres = "/home/dim/FIT/data/NEWurqmd.root", const char *current_mult = "hNpart", const char *outadres = "/home/dim/FIT/outfile.root",int minNch=25)
{
    Nch0=minNch;
    TFile *file = new TFile(fileadres);
    TH1D *Gev = (TH1D *)file->Get(current_mult);
    //file->Close();
    bin_cent[0] = 1.1 * Gev->FindLastBinAbove();
    Nn = bin_cent[0];
    Gev->Scale(1 / Gev->Integral(NO, Gev->GetNbinsX(), "width"));
    Gev->SetTitle("");
    Gev->GetYaxis()->SetTitle("1/N dN_{ch}/dN");
    Gev->GetXaxis()->SetTitle("N_{ch}");
    Gev->GetXaxis()->SetRangeUser(0, Nn);
    Gev->GetYaxis()->SetRangeUser(1e-6, 0.5);
   
    // Y axis ratio plot settings
    Gev->GetYaxis()->SetNdivisions(505);
    Gev->GetYaxis()->SetTitleSize(18);
    Gev->GetYaxis()->SetTitleFont(43);
    Gev->GetYaxis()->SetTitleOffset(1.2);
    Gev->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    Gev->GetYaxis()->SetLabelSize(18);
    // X axis ratio plot settings
    Gev->GetXaxis()->SetLabelSize(0);
    TH1F *GevC = (TH1F *)Gev->Clone();
    GevC->GetXaxis()->SetTitleSize(18);
    GevC->GetXaxis()->SetTitleFont(43);
    GevC->GetXaxis()->SetTitleOffset(1.5);
    GevC->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    GevC->GetXaxis()->SetLabelSize(18);
    GevC->GetXaxis()->SetTitleOffset(3);
    GevC->GetYaxis()->SetTitle("ratio data/fit");
    GevC->GetXaxis()->SetTitle("N_{ch}");

    //Задаем функцию для фитирования Gaus и строим верхнию половину графика
    TF1 *fc2 = new TF1("fit_func", ftPn, Nch0, bin_cent[0], 5);
    fc2->SetParameters(teta, bin_cent[0] * 0.6, a1, a2, a3);
    TCanvas *c2 = new TCanvas("data fit result", "ratio data/fit", 650, 500);
    TPad *pad12 = new TPad("pad12", "pad12", 0, 0.31, 1, 1.0);
    pad12->SetBottomMargin(0.05);
    pad12->Draw();
    pad12->cd()->SetLogy();
    Gev->Fit(fc2, "RM");
    Gev->Draw();

    teta = (fc2->GetParameter(0));
    n_knee = (fc2->GetParameter(1));
    a1 = (fc2->GetParameter(2));
    a2 = (fc2->GetParameter(3));
    a3 = (fc2->GetParameter(4));
cout <<endl<<"Fit results"<<endl;
cout<<"teta="<<teta<<", n_knee="<<n_knee<<", a1="<<a1<<" ,a2="<<a2<<", a3="<<a3<<endl<<endl;
    TLine *line2 = new TLine(n_knee, 1e-6, n_knee, 0.5);
    line2->Draw("SAME");

    c2->cd(); // Go back to the main canvas before defining pad2
    TPad *pad22 = new TPad("pad22", "pad22", 0, 0, 1, 0.3);
    pad22->SetBottomMargin(0.3);
    pad22->SetTopMargin(0.01);
    pad22->Draw();
    pad22->cd();
    TF1 *NEWfc2 = new TF1;
    NEWfc2 = Gev->GetFunction("fit_func");
    GevC->Divide(NEWfc2);
    GevC->SetMinimum(0.4); // Define Y ..
    GevC->SetMaximum(1.6); // .. range
    GevC->SetStats(0);     // No statistics on lower plot
    GevC->Draw();
    TLine *line = new TLine(0, 1, Nn, 1);
    line->Draw("SAME");
    file->Close();
    d_outfile = new TFile(outadres,"recreate");
    c2->Write();
    fc2->Write();
}
double Pb(double *b, double *Nch)
{
    // Parameters
    double n0 = Nch[0];
    double nn = Nch[1];
    // Variables
    double cb = pi * b[0] * b[0] / sigma;
    // Function
    TF1 *Pnb = new TF1("Pnb", PnbGaus2, 0, Nn, 7);
    Pnb->SetParameters(cb, teta, n_knee, a1, a2, a3, a4);
    double IPnb = Pnb->Integral(n0, nn);
    TF1 *Pn = new TF1("Pn", ftPn, 0, Nn, 6);
    Pn->SetParameters(teta, n_knee, a1, a2, a3, a4);
    double IPn = Pn->Integral(n0, nn);

    return 2 * pi * b[0] * IPnb / (sigma * IPn);
}
double bmean(double n0, double nn)
{
    TF1 *fPb = new TF1("fPb", Pb, 0, bmax, 2);
    fPb->SetParameters(n0, nn);
    return fPb->Mean(0, bmax);
}

double bsigma(double n0, double nn)
{
    double mean_b = 0, mean_b2 = 0;
    TF1 *fPb = new TF1("fPb", Pb, 0, bmax, 2);
    fPb->SetParameters(n0, nn);
    mean_b = fPb->Mean(0, bmax);
    mean_b2 = fPb->Moment(2, 0, bmax);
    return sqrt(mean_b2 - (mean_b * mean_b));
}

double integ(double n0, double nn)
{
    TF1 *Pn = new TF1("Pn", ftPn, 0, Nn, 5);
    Pn->SetParameters(teta, n_knee, a1, a2, a3);
    double IPn = Pn->Integral(n0, nn);
    return IPn;
}

void Rebin(double N0)
{
    double cen[10] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
    int n0 = bin_cent[0] - 1;
    double integr = 0, norm;
    norm = integ(N0, bin_cent[0]);
    for (int i = 0; i < 10; i++)
    {
        while (integr < cen[i])
        {

            integr = integ(n0, bin_cent[0]) / norm;
            n0 = n0 - 1;
        }
        bin_cent[i + 1] = n0 + 1;
    }
cout<<"Centrality classes"<<endl;
    for (int i = 0; i < 11; i++)
    {
        cout << "bin" << i << " " << bin_cent[i] << " cent " << integ(bin_cent[i], bin_cent[0]) / norm << endl;
    }
}

void PlotMeanb()
{
    Rebin(NO);
    double Bmean[10], Bsigma[10], cent[10], centEr[10];
    cout<<endl<<"Mean B vs Centrality"<<endl;
    for (int i = 0; i < 10; i++)
    {
        Bmean[i] = bmean(bin_cent[i + 1], bin_cent[i]);
        Bsigma[i] = bsigma(bin_cent[i + 1], bin_cent[i]);
        cent[i] = i * 10 + 5;
        centEr[i] = 0;
        cout << Bmean[i] << " +- " << Bsigma[i] << endl;
    }

    TGraphErrors *Gr = new TGraphErrors(10, cent, Bmean, centEr, Bsigma);
    Gr->SetMarkerStyle(25);
    Gr->SetMarkerSize(1.5);
    Gr->SetMarkerColor(kBlack);
    Gr->SetLineColor(kBlack);
    Gr->SetLineWidth(2.);
    Gr->SetTitle("");
    Gr->GetXaxis()->SetTitle("Centrality, %");
    Gr->GetYaxis()->SetTitle("<b>, fm");
    Gr->SetTitle("B vs Centraliry");
    Gr->SetName("Fit_B_Mean");
    TCanvas *c = new TCanvas("canvas fit b", "B vs Centraliry", 650, 500);
    d_outfile->cd();
    Gr->Draw();
    Gr->Write();
    //c->Write();
    d_outfile->Close();
}

void GammaFit(const char *fileadres = "/home/dim/FIT/data/NEWurqmd.root", const char *current_mult = "hNpart", const char *outadres = "/home/dim/FIT/outfile.root",int minNch=25)
{
    Start(fileadres ,current_mult , outadres ,minNch);
    PlotMeanb();
}
