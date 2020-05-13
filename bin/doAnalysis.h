#ifndef process_h
#define process_h

#include <vector>
#include <map>
#include <tuple>

#include "trktree.h"
#include "rooutil.h"
#include "cxxopts.h"

#include "SDL/Event.cuh" // SDL::Event
#include "SDL/Module.cuh" // SDL::Module
#include "SDL/PrintUtil.h" // SDL::out
#include "SDL/EndcapGeometry.h" // SDL::EndcapGeometry
#include "SDL/ModuleConnectionMap.h" // SDL::ModuleConnectionMap

// Efficiency study modules
#include "Study.h"
#include "StudyHitOccupancy.h"

#include "constants.h"

#include "AnalysisConfig.h"

//#include "trkCore.h"

void printModuleConnectionInfo(std::ofstream&);
// bool hasAll12HitsInBarrel(unsigned int);
void processModuleBoundaryInfo();
int layer(int lay, int det);

#endif
