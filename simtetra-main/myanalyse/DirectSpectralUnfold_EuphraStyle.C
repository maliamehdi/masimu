#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TF1.h>
#include <iostream>

// -----------------------------------------------------------------------------
// DirectSpectralUnfold_EuphraStyle
// Reproduction FIDÈLE de la méthode Euphra (Oberstedt/Billnert style)
// avec indexation C 0-based, décalages +1, utilisation de GetBinContent(j,i),
// fit de la queue exponentielle optionnel.
// -----------------------------------------------------------------------------
TH1D* DirectSpectralUnfold_EuphraStyle(const TH1*  hMeasIn,
                                       const TH2*  hRespIn,
                                       bool doTailFit         = true,
                                       double fitMin          = 2000,
                                       double fitMax          = 8000,
                                       int    firstRealBin    = 1,
                                       int    lastPhysicalBin = -1)
{
    if (!hMeasIn || !hRespIn) {
        std::cerr << "[DirectSpectralUnfold_EuphraStyle] ERROR null histos.\n";
        return nullptr;
    }

    // -------------------------------------------------------
    // Déterminer nombre de bins (attention: +2 pour under/overflow)
    // -------------------------------------------------------
    int NX = hRespIn->GetNbinsX();   // bins physiques X (measured)
    int NY = hRespIn->GetNbinsY();   // bins physiques Y (true)
    int N  = NY + 1;                 // Comme Euphra : tableau C de 0..NY

    if (lastPhysicalBin < 0) lastPhysicalBin = NY; // souvent = dernier bin de Y

    // -------------------------------------------------------
    // Construire tableau source[i] (spectre mesuré)
    // -------------------------------------------------------
    std::vector<double> source(N,0.0);

    for (int i = 0; i < N; i++) {
        // IMPORTANT : GetBinContent(i) avec i=0 récupère l’underflow ROOT
        source[i] = hMeasIn->GetBinContent(i);
    }

    // -------------------------------------------------------
    // Construire tableau réponse response[i][j]
    // Rappel : dans Euphra :
    //   response[i][j] = GetBinContent(j,i)
    //   => i = axe Y, j = axe X
    //   => j=0 → underflow X
    // -------------------------------------------------------
    std::vector<std::vector<double>> response(N,
                std::vector<double>(N+1, 0.0));

    for (int i = 0; i < N; i++) {        // ligne Y
        for (int j = 0; j < N+1; j++) {  // colonne X
            response[i][j] = hRespIn->GetBinContent(j,i);
        }
    }

    // -------------------------------------------------------
    // FIT DE LA QUEUE EXPONENTIELLE (optionnel)
    // -------------------------------------------------------
    if (doTailFit) {

        TF1* expo = new TF1("expoTail","expo", fitMin, fitMax);
        ((TH1*)hMeasIn)->Fit(expo, "RQ");  // quiet

        double b = expo->GetParameter(0);
        double a = expo->GetParameter(1);

        // Remplir les "trous"
        for (int i = firstRealBin; i <= lastPhysicalBin; ++i) {
            double x = hMeasIn->GetBinCenter(i);
            double val = TMath::Exp(b + a*x);
            source[i] = val;
        }
    }

    // -------------------------------------------------------
    // HISTOGRAMME DE SORTIE (spectre déplié)
    // -------------------------------------------------------
    TH1D* hOut = new TH1D("hUnfold_EuphraStyle",
                          "Direct unfolding (Euphra style);E_{true};Arb.",
                          NY,
                          hRespIn->GetYaxis()->GetXbins()->GetArray());
    hOut->SetDirectory(nullptr);

    // -------------------------------------------------------
    // ALGORITHME PRINCIPAL : du bin le plus haut vers 1
    // -------------------------------------------------------
    for (int i = lastPhysicalBin; i >= firstRealBin; --i)
    {
        if (source[i] == 0) continue;

        // Normalisation : alpha = source[i] / response[i][i+1]
        double denom = response[i][i+1];
        if (denom <= 0) continue;

        double alpha = source[i] / denom;

        // Stocker dans le spectre vrai
        hOut->SetBinContent(i, alpha);

        // Construire spectre simulé (ligne i)
        std::vector<double> normSpec(i+1, 0.0);
        for (int j = 0; j < i+1; j++) {
            normSpec[j] = response[i][j+1] * alpha;
        }

        // Soustraire du résiduel
        for (int j = 0; j < i; j++) {
            source[j] -= normSpec[j];
            if (source[j] < 1e-12) source[j] = 0.0;
        }
    }

    return hOut;
}