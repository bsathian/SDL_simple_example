#include "StudyHitOccupancy.h"

StudyHitOccupancy::StudyHitOccupancy(const char* studyName, StudyHitOccupancy::StudyHitOccupancyMode mode_)
{

    studyname = studyName;
    mode = mode_;
    switch (mode)
    {
        case kStudyAll: modename = "all"; break;
        default: modename = "UNDEFINED"; break;
    }

}

void StudyHitOccupancy::bookStudy()
{
    // Book Histograms
    const int nlayers = NLAYERS;
    for (int ii = 0; ii < nlayers; ++ii)
    {
        ana.histograms.addHistogram(TString::Format("nhits_lower_by_layer%d", ii), 56, 0, 14000, [&, ii]() { return nhits_in_lower_modules_by_layer[ii]; } );
        ana.histograms.addHistogram(TString::Format("nhits_upper_by_layer%d", ii), 56, 0, 14000, [&, ii]() { return nhits_in_upper_modules_by_layer[ii]; } );
        ana.histograms.addHistogram(TString::Format("nhits_both_by_layer%d", ii),112, 0, 28000, [&, ii]() { return nhits_in_both_modules_by_layer[ii]; } );
    }

}

void StudyHitOccupancy::doStudy(SDL::Event& event, std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents)
{

    // Each do study is performed per event

    // First clear all the output variables that will be used to fill the histograms for this event
    for (int ii = 0; ii < NLAYERS; ++ii)
    {
        nhits_in_lower_modules_by_layer[ii] = 0;
        nhits_in_upper_modules_by_layer[ii] = 0;
        nhits_in_both_modules_by_layer[ii] = 0;
    }

    //***********************
    // Studying Hit Occupancy
    //***********************

    std::vector<SDL::Module*> moduleList = event.getModulePtrs();

    // Loop over tracklets in event
    for (auto& modulePtr : moduleList)
    {
        int subdet = modulePtr->subdet();
        if (subdet != 5)
            continue;

        // Get the nhit
        int nhit_in_module = modulePtr->getHitPtrs().size();

        // Get the layer index
        int ilayer = modulePtr->layer();

        // Is Lower
        int isLower = modulePtr->isLower();

        // Add to the counter
        if (isLower)
            nhits_in_lower_modules_by_layer[ilayer-1] += nhit_in_module;
        else
            nhits_in_upper_modules_by_layer[ilayer-1] += nhit_in_module;
        nhits_in_both_modules_by_layer[ilayer-1] += nhit_in_module;

    }

    // First clear all the output variables that will be used to fill the histograms for this event
    for (int ii = 0; ii < NLAYERS; ++ii)
    {
        // std::cout <<  " nhits_in_lower_modules_by_layer[ii]: " << nhits_in_lower_modules_by_layer[ii] <<  std::endl;
        // std::cout <<  " nhits_in_upper_modules_by_layer[ii]: " << nhits_in_upper_modules_by_layer[ii] <<  std::endl;
        // std::cout <<  " nhits_in_both_modules_by_layer[ii]: " << nhits_in_both_modules_by_layer[ii] <<  std::endl;
    }

}
