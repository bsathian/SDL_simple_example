#ifndef StudyHitOccupancy_h
#define StudyHitOccupancy_h

#include "SDL/Event.cuh"

#include "Study.h"

#include <vector>
#include <tuple>

#include "TString.h"
#include "trktree.h"
#include "constants.h"

#include "AnalysisConfig.h"

class StudyHitOccupancy : public Study
{

public:
    enum StudyHitOccupancyMode
    {
        kStudyAll = 1,
    };

    const char* studyname;
    StudyHitOccupancyMode mode;
    const char* modename;
    std::array<int, NLAYERS> nhits_in_lower_modules_by_layer;
    std::array<int, NLAYERS> nhits_in_upper_modules_by_layer;
    std::array<int, NLAYERS> nhits_in_both_modules_by_layer;

    StudyHitOccupancy(const char* studyName, StudyHitOccupancy::StudyHitOccupancyMode mode_);
    virtual void bookStudy();
    virtual void doStudy(SDL::Event& recoevent, std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents);

};

#endif
