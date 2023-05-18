#ifndef ROOT_TFoamMT
#define ROOT_TFoamMT

#include <BS_thread_pool_light.hpp>
#include <Riostream.h>
#include <TFoam.h>
#include <TFoamCell.h>
#include <TFoamVect.h>
#include <TH1.h>
#include <TMutex.h>
#include <TRandom.h>
#include <TRefArray.h>
#include <thread>

class TFoamMT : public TFoam {
    using TFoam::TFoam;

public:
    TFoamMT(const std::string &name, Int_t nbThreads) : TFoam(name.c_str()), NbThreads(nbThreads), pool(nbThreads) {}
    void Explore(TFoamCell *) override;
    void MakeAlpha() override;
    void MakeLocalAlpha(std::vector<Double_t> &alpha);
    void LocalMCLoop(TFoamCell *, Int_t);

private:
    Int_t NbThreads;
    Int_t fEvMT = 0;

    std::atomic<Double_t> nevEff;
    std::array<Double_t, 5> fSum = {};

    TMutex MutexSum;
    TMutex MutexHistEdg;

    BS::thread_pool_light pool;
};

#endif// ROOT_TFoamMT
