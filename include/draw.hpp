#ifndef NOPEG_ANALYSER_DRAW_HPP
#define NOPEG_ANALYSER_DRAW_HPP

#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPaveStats.h>
#include <TText.h>
#include <memory>


namespace Draw {
    namespace Axis {
        enum axis { X,
                    Y };
    }

    struct Color {
        [[maybe_unused]] inline static short lightBlue  = static_cast<short>(TColor::GetColor("#0C5DA5"));
        [[maybe_unused]] inline static short grey       = static_cast<short>(TColor::GetColor("#787878"));
        [[maybe_unused]] inline static short red        = static_cast<short>(TColor::GetColor("#FF2C00"));
        [[maybe_unused]] inline static short lightGreen = static_cast<short>(TColor::GetColor("#00B945"));
        [[maybe_unused]] inline static short yellow     = static_cast<short>(TColor::GetColor("#FF9500"));
        [[maybe_unused]] inline static short violet     = static_cast<short>(TColor::GetColor("#845B97"));
        [[maybe_unused]] inline static short black      = static_cast<short>(TColor::GetColor("#121415"));
        [[maybe_unused]] inline static short green      = static_cast<short>(TColor::GetColor("#006F29"));
    };


    [[nodiscard]] auto makeCanvas(const char *name = "", const char *title = "", const int ww = 800, const int wh = 600, const float bottomMargin = 0.12,
                                  const float topMargin = 0.12) -> std::shared_ptr<TCanvas> {
        auto canvas = std::make_shared<TCanvas>(name, title, ww, wh);
        canvas->SetLeftMargin(bottomMargin);
        canvas->SetBottomMargin(topMargin);
        return canvas;
    }

    template<typename T>
    auto setHistColor(const T &h, const short color, const float alpha = 0) -> void {
        h->SetLineColor(color);
        if (alpha > 0) h->SetFillColorAlpha(color, alpha);
    }

    template<typename T>
    auto setAxisLabel(const T &h, const Axis::axis axis, const std::string &labelAxis, const float offset = 0.8, const float size = 0.05) -> void {
        if (axis == Axis::X) {
            h->GetXaxis()->SetTitle(labelAxis.c_str());
            h->GetXaxis()->SetTitleOffset(offset);
            h->GetXaxis()->SetTitleSize(size);
        } else if (axis == Axis::Y) {
            h->GetYaxis()->SetTitle(labelAxis.c_str());
            h->GetYaxis()->SetTitleOffset(offset);
            h->GetYaxis()->SetTitleSize(size);
        }
    }
    auto tileStats(const std::unique_ptr<TObjArray> &boxes, const float X1NDC, const float X2NDC, const float Y1NDC, const float Y2NDC) -> void {
        const int N        = boxes->GetEntries();
        const int n_across = TMath::CeilNint(std::sqrt(N));
        const int n_down   = TMath::CeilNint(static_cast<float>(N) / static_cast<float>(n_across));
        const float width  = (X2NDC - X1NDC) / static_cast<float>(n_across);
        const float height = (Y2NDC - Y1NDC) / static_cast<float>(n_down);

        int across = 0;
        int down   = 0;
        for (Int_t i = 0; i < N; i++) {
            auto tps = dynamic_cast<TPaveStats *>(boxes->At(i));
            if (!tps) return;

            tps->SetX1NDC(X2NDC - static_cast<float>(across + 1) * width);
            tps->SetX2NDC(X2NDC - static_cast<float>(across) * width);
            tps->SetY1NDC(Y2NDC - static_cast<float>(down + 1) * height);
            tps->SetY2NDC(Y2NDC - static_cast<float>(down) * height);

            across++;
            if (across >= n_across) {
                across = 0;
                down++;
            }
        }
    }

    auto tileStats(const std::unique_ptr<TObjArray> &boxes) -> void {
        const int N           = boxes->GetEntries();
        const int n_across    = TMath::CeilNint(TMath::Sqrt(N));
        const int n_down      = TMath::CeilNint((Double_t) N / (Double_t) n_across);
        const float defX1     = 8.20e-01;
        const float defX2     = 9.80e-01;
        const float defY1     = 8.35e-01;
        const float defY2     = 9.95e-01;
        const float defWidth  = defX2 - defX1;
        const float defHeight = defY2 - defY1;

        tileStats(boxes, defX2 - static_cast<float>(n_across) * defWidth, defX2, defY2 - static_cast<float>(n_down) * defHeight, defY2);
    }

    template<typename T>
    auto setStatBoxes(const std::shared_ptr<TCanvas> &canvas, const std::initializer_list<T> &v) -> void {
        canvas->Update();
        const auto boxes = std::make_unique<TObjArray>();
        for (const auto &h: v) {
            const auto st = dynamic_cast<TPaveStats *>(h->GetListOfFunctions()->FindObject("stats"));
            st->SetTextColor(h->GetLineColor());
            st->Draw();
            boxes->Add(st);
        }
        tileStats(boxes);
        canvas->Modified();
    }

    struct ArgsH {
        const std::string &fileName    = "";
        const std::shared_ptr<TF1> &f1 = nullptr;
        const float x_min              = 0;
        const float x_max              = 0;
        const bool setStats            = true;
        const bool logx                = false;
        const bool logy                = false;
    };

    template<class H>
    auto draw(const H &h, const std::string &path, const std::string &label_xAxis, const ArgsH &args1 = {}) -> void {
        const auto canvas      = makeCanvas();
        std::string myFileName = args1.fileName;
        if (args1.fileName.empty()) myFileName = h->GetName();

        if (args1.f1) h->Fit(args1.f1.get(), "", "", args1.x_min, args1.x_max);
        h->Draw();
        h->SetMinimum(0);
        h->SetStats(args1.setStats);

        if (args1.logx) canvas->SetLogx();
        if (args1.logy) canvas->SetLogy();

        setHistColor(h, Color::lightBlue, 0.15);
        setAxisLabel(h, Axis::X, label_xAxis, 1);
        if (args1.setStats) setStatBoxes(canvas, {h});
        canvas->SaveAs((path + myFileName + ".pdf").c_str());
    }

    template<class H>
    auto draw(const H &h1, const H &h2, const std::string &path, const std::string &label_xAxis, const ArgsH &args1 = {}) -> void {
        const auto canvas      = makeCanvas();
        const auto hs          = std::make_shared<THStack>("hs", "");
        const auto h2Copy      = std::shared_ptr<TH1F>(dynamic_cast<TH1F *>(h2->Clone()));
        std::string myFileName = args1.fileName;
        if (args1.fileName.empty()) myFileName = h1->GetName();

        setHistColor(h1, Color::red);
        setHistColor(h2Copy, Color::lightBlue);

        h2Copy->SetLineStyle(2);

        h1->Scale(1. / h1->Integral(), "width");
        h2Copy->Scale(1. / h2Copy->Integral(), "width");

        hs->Add(h1.get(), "sames");
        hs->Add(h2Copy.get(), "sames");
        hs->Draw("hist nostack");

        setAxisLabel(hs, Axis::X, label_xAxis);
        setStatBoxes(canvas, {h1, h2Copy});

        if (args1.logx) canvas->SetLogx();
        if (args1.logy) canvas->SetLogy();

        auto legend = std::make_unique<TLegend>();
        legend->AddEntry(h1.get(), "old", "l");
        legend->AddEntry(h2Copy.get(), "new", "l");
        legend->Draw();

        canvas->SaveAs((path + myFileName + ".pdf").c_str());
    }
}// namespace Draw

#endif//NOPEG_ANALYSER_DRAW_HPP
