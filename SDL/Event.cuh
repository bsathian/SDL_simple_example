#ifndef Event_h
#define Event_h

#include <vector>
#include <list>
#include <map>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <cuda_runtime.h>

#include "Module.cuh"
#include "Hit.cuh"
#include "MiniDoublet.cuh"
#include "PrintUtil.h"
#include "Algo.h"
#include "ModuleConnectionMap.h"
#include "GeometryUtil.cuh"

namespace SDL
{
    class Event
    {
        private:

            // map of modules (this holds the actual instances)
            std::map<unsigned int, Module*> modulesMapByDetId_;

            // list of hits (this holds the actual instances)
            std::list<Hit> hits_;

            // list of hit_boundaries (this holds the actual instances) this is only used for the 2S in the endcap
            std::list<Hit> hits_2s_edges_;

            // list of MiniDoublets (this holds the actual instances)
            std::list<MiniDoublet> miniDoublets_;

            // list of module pointers (hold only the pointers to the actual instances)
            std::vector<Module*> modulePtrs_;

            // list of lower module pointers (hold only the pointers to the actual instances)
            // (lower means, the module that is closer to the luminous region)
            std::vector<Module*> lowerModulePtrs_;

            // boolean to turn on debug mode
            SDL::LogLevel logLevel_;

            // diagnostic variables
            // # of hits in barrel
            std::array<unsigned int, 6> n_hits_by_layer_barrel_;

            // # of hits in endcap
            std::array<unsigned int, 5> n_hits_by_layer_endcap_;

            // # of hits in barrel in upper module
            std::array<unsigned int, 6> n_hits_by_layer_barrel_upper_;

            // # of hits in endcap in upper module
            std::array<unsigned int, 5> n_hits_by_layer_endcap_upper_;

            // # of pairs of hits considered for mini-doublet
            std::array<unsigned int, 6> n_miniDoublet_candidates_by_layer_barrel_;

            // # of pairs of hits considered for mini-doublet
            std::array<unsigned int, 6> n_miniDoublet_by_layer_barrel_;

            // # of pairs of hits considered for mini-doublet
            std::array<unsigned int, 5> n_miniDoublet_candidates_by_layer_endcap_;

            // # of pairs of hits considered for mini-doublet
            std::array<unsigned int, 5> n_miniDoublet_by_layer_endcap_;

            // Multiplicity of mini-doublet candidates considered in this event
            void incrementNumberOfHits(SDL::Module& module);

            // Multiplicity of mini-doublet candidates considered in this event
            void incrementNumberOfMiniDoubletCandidates(SDL::Module& module);

            // Multiplicity of mini-doublet formed in this event
            void incrementNumberOfMiniDoublets(SDL::Module& module);
            
            //CUDA stuff
            Module* modulesInGPU;
            Hit* hitsInGPU;
            Hit* hits2sEdgeInGPU;
            MiniDoublet* mdCandsGPU;
            MiniDoublet* mdsInGPU;
            int mdGPUCounter;
            void  initModulesInGPU();
            void initHitsInGPU();
            void initMDsInGPU();
            void miniDoubletGPUWrapper(SDL::MDAlgo algo);

            //counter variables
            int moduleMemoryCounter;
            int hitMemoryCounter;
            int hit2SEdgeMemoryCounter;
            int mdMemoryCounter;
 

        public:

            // cnstr/destr
            Event();
            ~Event();

            // Module related functions
            bool hasModule(unsigned int detId);
            Module* getModule(unsigned int detId);
            const std::vector<Module*> getModulePtrs() const;
            const std::vector<Module*> getLowerModulePtrs() const;


            // Set debug
            void setLogLevel(SDL::LogLevel logLevel=SDL::Log_Nothing);

            // Hit related functions
            void addHitToModule(Hit hit, unsigned int detId);

            // MiniDoublet related functions
            void addMiniDoubletToLowerModule(MiniDoublet md, unsigned int detId);
            void addMiniDoubletToEvent(MiniDoublet md, unsigned int detId);

            // Create mini doublets
            void createMiniDoublets(MDAlgo algo=Default_MDAlgo);

            // Create mini doublet for a module
            void createMiniDoubletsFromLowerModule(unsigned int detId, int maxMDCands, MDAlgo algo=Default_MDAlgo);


            // Multiplicity of Hits
            unsigned int getNumberOfHits();

            // Multiplicity of hits in this event
            unsigned int getNumberOfHitsByLayerBarrel(unsigned int);

            // Multiplicity of hits in this event
            unsigned int getNumberOfHitsByLayerEndcap(unsigned int);

            // Multiplicity of hits in this event for upper module
            unsigned int getNumberOfHitsByLayerBarrelUpperModule(unsigned int);

            // Multiplicity of hits in this event for upper module
            unsigned int getNumberOfHitsByLayerEndcapUpperModule(unsigned int);

            // Multiplicity of mini-doublets
            unsigned int getNumberOfMiniDoublets();

            // Multiplicity of mini-doublet candidates considered in this event
            unsigned int getNumberOfMiniDoubletCandidates();

            // Multiplicity of mini-doublet candidates considered in this event
            unsigned int getNumberOfMiniDoubletCandidatesByLayerBarrel(unsigned int);

            // Multiplicity of mini-doublet candidates considered in this event
            unsigned int getNumberOfMiniDoubletCandidatesByLayerEndcap(unsigned int);

            // Multiplicity of mini-doublet formed in this event
            unsigned int getNumberOfMiniDoubletsByLayerBarrel(unsigned int);

            // Multiplicity of mini-doublet formed in this event
            unsigned int getNumberOfMiniDoubletsByLayerEndcap(unsigned int);

            // cout printing
            friend std::ostream& operator<<(std::ostream& out, const Event& event);
            friend std::ostream& operator<<(std::ostream& out, const Event* event);

    };
}
//CUDA Kernels
__global__ void createMiniDoubletsInGPU(SDL::MiniDoublet* mdCands, int n, SDL::MDAlgo algo);

#endif
