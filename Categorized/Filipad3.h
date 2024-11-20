#ifndef FILIPAD3_H
#define FILIPAD3_H

#include <iostream>
#include "TFile.h"
#include "TStyle.h"
#include "TH1D.h"

#include "TCanvas.h"
#include "TPad.h"
#include "TString.h"
#include "TLegend.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TLatex.h"

class Filipad3 {
public:
    Filipad3(int rows, int cols, int inID=1, float inRelSize=1.1, float inR = 0.4, int inXOffset = 100, int inYOffset=100, float inAspect=0.7, int ichop=5, int ichopPt=3);
    ~Filipad3();

    void DrawCanvas();
    TPad* GetPad(int row, int col);

private:
    int numRows;
    int numCols;
    int ID;
    float aspektCanvas;
    int sizeCanvas;
    int sdxCanvas;
    int sdyCanvas;
    float MarginLeft;
    float MarginBottom;
    float MarginRight;
    float MarginTop;
    TString mcpad;
    TCanvas *canvas;
    TPad ***pads;

    // Filipad2 settings
    float space;
    float ratio;
};

Filipad3::Filipad3(int rows, int cols, int inID, float inRelSize, float inR, int inXOffset, int inYOffset, float inAspect, int ichop, int ichopPt)
    : numRows(rows), numCols(cols), ID(inID), aspektCanvas(inAspect), sizeCanvas(300 * inRelSize),
      sdxCanvas(0), sdyCanvas(0), MarginLeft(0.15), MarginBottom(0.08), MarginRight(0.03), MarginTop(0.02),
      mcpad(Form("c%d", inID)), canvas(nullptr), pads(nullptr),
      space(0), ratio(inR)
{
    canvas = new TCanvas(mcpad, mcpad, sdxCanvas, sdyCanvas, sizeCanvas * aspektCanvas, sizeCanvas);
    canvas->SetName(Form("%s", mcpad.Data()));
    canvas->SetFillStyle(4000);
    canvas->SetFillColor(10);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    canvas->SetTopMargin(0.);
    canvas->SetBottomMargin(0.);

    pads = new TPad**[numRows];
    float padWidth = 1.0 / numCols;
    float padHeight = 1.0 / numRows;

    for (int row = 0; row < numRows; row++) {
        pads[row] = new TPad*[numCols];
        for (int col = 0; col < numCols; col++) {
            TString padName = Form("pad_%d_%d", row, col);
            float x1 = col * padWidth;
            float y1 = 1.0 - (row + 1) * padHeight;
            float x2 = (col + 1) * padWidth;
            float y2 = 1.0 - row * padHeight;

            pads[row][col] = new TPad(padName, padName, x1, y1, x2, y2);
            pads[row][col]->SetTopMargin(MarginTop / (1 - ratio));
            pads[row][col]->SetBottomMargin(0.0015);
            pads[row][col]->SetLeftMargin(MarginLeft);
            pads[row][col]->SetRightMargin(MarginRight);
        }
    }
}

Filipad3::~Filipad3() {
    for (int row = 0; row < numRows; row++) {
        for (int col = 0; col < numCols; col++) {
            delete pads[row][col];
        }
        delete[] pads[row];
    }
    delete[] pads;
    delete canvas;
}

void Filipad3::DrawCanvas() {
    canvas->cd();

    for (int row = 0; row < numRows; row++) {
        for (int col = 0; col < numCols; col++) {
            pads[row][col]->Draw();
        }
    }
}

TPad* Filipad3::GetPad(int row, int col) {
    if (row >= 0 && row < numRows && col >= 0 && col < numCols) {
        return pads[row][col];
    } else {
        return nullptr;
    }
}

#endif // FILIPAD3_H
