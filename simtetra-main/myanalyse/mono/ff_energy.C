#ifdef __CLING__
#pragma cling optimize(0)
#endif
void ff_energy()
{
//=========Macro generated from canvas: cEff/Efficiency vs neutron energy
//=========  (Wed Dec  3 15:50:38 2025) by ROOT version 6.28/06
   TCanvas *cEff = new TCanvas("cEff", "Efficiency vs neutron energy",561,609,900,600);
   gStyle->SetOptStat(0);
   cEff->Range(-1.325419,4.009142,1.471196,23.4606);
   cEff->SetFillColor(0);
   cEff->SetBorderMode(0);
   cEff->SetBorderSize(2);
   cEff->SetLogx();
   cEff->SetGridx();
   cEff->SetGridy();
   cEff->SetFrameBorderMode(0);
   cEff->SetFrameBorderMode(0);
   
   Double_t Graph0_fx1002[16] = { 0.09964721, 0.4966631, 0.9957405, 1.996321, 2.982838, 3.97375, 4.998701, 5.980176, 6.952062, 7.966816, 8.999687, 9.950094, 10.92227, 11.98942, 12.97344, 13.93792 };
   Double_t Graph0_fy1002[16] = { 20.08361, 20.08361, 19.10429, 17.41579, 15.59222, 14.61289, 12.9244, 12.18146, 10.99951, 10.39165, 9.142167, 8.365459, 8.432999, 7.89268, 7.487442, 7.352362 };
   Double_t Graph0_fex1002[16] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   Double_t Graph0_fey1002[16] = { 0.1416933, 0.1416192, 0.1381593, 0.1319621, 0.1248359, 0.1208429, 0.1135473, 0.1103494, 0.1048666, 0.1019166, 0.09548298, 0.09144944, 0.09179325, 0.08871866, 0.08658522, 0.08594766 };
   TGraphErrors *gre = new TGraphErrors(16,Graph0_fx1002,Graph0_fy1002,Graph0_fex1002,Graph0_fey1002);
   gre->SetName("Graph0");
   gre->SetTitle("Detection efficiency");
   gre->SetFillStyle(1000);
   gre->SetLineColor(51);
   gre->SetLineWidth(2);
   gre->SetMarkerColor(6);
   gre->SetMarkerStyle(6);
   gre->SetMarkerSize(1.1);
   
   TH1F *Graph_Graph01002 = new TH1F("Graph_Graph01002","Detection efficiency",100,0.09,15.39);
   Graph_Graph01002->SetMinimum(5.954288);
   Graph_Graph01002->SetMaximum(21.51546);
   Graph_Graph01002->SetDirectory(nullptr);
   Graph_Graph01002->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph01002->SetLineColor(ci);
   Graph_Graph01002->GetXaxis()->SetTitle("Neutron energy (MeV)");
   Graph_Graph01002->GetXaxis()->SetRange(0,101);
   Graph_Graph01002->GetXaxis()->SetMoreLogLabels();
   Graph_Graph01002->GetXaxis()->SetNoExponent();
   Graph_Graph01002->GetXaxis()->SetLabelFont(42);
   Graph_Graph01002->GetXaxis()->SetTitleOffset(1);
   Graph_Graph01002->GetXaxis()->SetTitleFont(42);
   Graph_Graph01002->GetYaxis()->SetTitle("Efficiency (%)");
   Graph_Graph01002->GetYaxis()->SetLabelFont(42);
   Graph_Graph01002->GetYaxis()->SetTitleFont(42);
   Graph_Graph01002->GetZaxis()->SetLabelFont(42);
   Graph_Graph01002->GetZaxis()->SetTitleOffset(1);
   Graph_Graph01002->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph01002);
   
   gre->Draw("APLE");
   
   TLegend *leg = new TLegend(0.6002227,0.7725694,0.9008909,0.9027778,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph0","Mean measured multiplicity (x100%)","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(6);
   entry->SetMarkerStyle(6);
   entry->SetMarkerSize(1.1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("NULL","N points = 16","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   
   TPaveText *pt = new TPaveText(0.3324276,0.9345833,0.6675724,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_LaTex = pt->AddText("Detection efficiency");
   pt->Draw();
   cEff->Modified();
   cEff->cd();
   cEff->SetSelected(cEff);
}
