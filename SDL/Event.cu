#include "Event.cuh"

//CUDA Kernel for Minidoublet creation
__global__ void createMiniDoubletsInGPU(SDL::MiniDoublet* mdCands, int n, SDL::MDAlgo algo)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for(int i = tid; i<n; i+= stride)
    {
        mdCands[i].runMiniDoubletAlgo(algo);
    }

}


SDL::Event::Event() : logLevel_(SDL::Log_Nothing)
{
    n_hits_by_layer_barrel_.fill(0);
    n_hits_by_layer_endcap_.fill(0);
    n_hits_by_layer_barrel_upper_.fill(0);
    n_hits_by_layer_endcap_upper_.fill(0);
    n_miniDoublet_candidates_by_layer_barrel_.fill(0);
    n_miniDoublet_by_layer_barrel_.fill(0);
    n_miniDoublet_candidates_by_layer_endcap_.fill(0);
    n_miniDoublet_by_layer_endcap_.fill(0);
    
}

SDL::Event::~Event()
{
}

bool SDL::Event::hasModule(unsigned int detId)
{
    if (modulesMapByDetId_.find(detId) == modulesMapByDetId_.end())
    {
        return false;
    }
    else
    {
        return true;
    }
}

void SDL::Event::setLogLevel(SDL::LogLevel logLevel)
{
    logLevel_ = logLevel;
}

void SDL::Event::initModulesInGPU()
{
    const int MODULE_MAX=50000;
    cudaMallocManaged(&modulesInGPU,MODULE_MAX * sizeof(SDL::Module));
}

SDL::Module* SDL::Event::getModule(unsigned int detId)
{
    // using std::map::emplace
    static int counter = 0;
    if(counter == 0)
    {
        initModulesInGPU();
    }
    std::pair<std::map<unsigned int, Module*>::iterator, bool> emplace_result = modulesMapByDetId_.emplace(detId,nullptr);
    // Retreive the module
    auto& inserted_or_existing = (*(emplace_result.first)).second;

    // If new was inserted, then insert to modulePtrs_ pointer list
    if (emplace_result.second) // if true, new was inserted
    {
        //cudaMallocManaged(&((*(emplace_result.first)).second),sizeof(SDL::Module));
         (*(emplace_result.first)).second = &modulesInGPU[counter];

        //*inserted_or_existing =SDL:: Module(detId);
        modulesInGPU[counter] = SDL::Module(detId);
        Module* module_ptr = inserted_or_existing;
        
        // Add the module pointer to the list of modules
        modulePtrs_.push_back(module_ptr);
        // If the module is lower module then add to list of lower modules
        if (module_ptr->isLower())
            lowerModulePtrs_.push_back(module_ptr);
       
       counter++;

    }

    return inserted_or_existing;
}

const std::vector<SDL::Module*> SDL::Event::getModulePtrs() const
{
    return modulePtrs_;
}

const std::vector<SDL::Module*> SDL::Event::getLowerModulePtrs() const
{
    return lowerModulePtrs_;
}

void SDL::Event::initHitsInGPU()
{
    const int HIT_MAX = 1000000;
    cudaMallocManaged(&hitsInGPU,HIT_MAX * sizeof(SDL::Hit));
    const int HIT_2S_MAX = 100000;
    cudaMallocManaged(&hits2sEdgeInGPU,HIT_2S_MAX * sizeof(SDL::Hit));
}


void SDL::Event::addHitToModule(SDL::Hit hit, unsigned int detId)
{
    // Add to global list of hits, where it will hold the object's instance
    static int counter = 0;
    static int counter2SEdge = 0;
    // And get the module (if not exists, then create), and add the address to Module.hits_
    //construct a cudaMallocManaged object and send that in, so that we won't have issues in the GPU
    if(counter == 0)
    {
        initHitsInGPU();
    }
    hitsInGPU[counter] = hit;
    hitsInGPU[counter].setModule(getModule(detId));
    getModule(detId)->addHit(&hitsInGPU[counter]);
    hits_.push_back(hitsInGPU[counter]);


    // Count number of hits in the event
    incrementNumberOfHits(*getModule(detId));

    // If the hit is 2S in the endcap then the hit boundary needs to be set
    if (getModule(detId)->subdet() == SDL::Module::Endcap and getModule(detId)->moduleType() == SDL::Module::TwoS)
    {
         
        hits2sEdgeInGPU[counter2SEdge] = GeometryUtil::stripHighEdgeHit(hitsInGPU[counter]);
        hits2sEdgeInGPU[counter2SEdge+1] = GeometryUtil::stripLowEdgeHit(hitsInGPU[counter]);
//        hits_2s_edges_.push_back(GeometryUtil::stripHighEdgeHit(&hits_.back()));
//        hits_.back().setHitHighEdgePtr(&(hits_2s_edges_.back()));
//        hits_2s_edges_.push_back(GeometryUtil::stripLowEdgeHit(*hitForGPU));
//        hits_.back().setHitLowEdgePtr(&(hits_2s_edges_.back()));
        hits_2s_edges_.push_back(hits2sEdgeInGPU[counter2SEdge]);
        hitsInGPU[counter].setHitHighEdgePtr(&hits2sEdgeInGPU[counter2SEdge]);

        hits_2s_edges_.push_back(hits2sEdgeInGPU[counter2SEdge+1]);
        hitsInGPU[counter].setHitLowEdgePtr(&hits2sEdgeInGPU[counter2SEdge+1]);

        counter2SEdge+= 2;
    }

    counter++;
}

void SDL::Event::createMiniDoublets(MDAlgo algo)
{
    // Loop over lower modules
    const int MAX_MD_CAND = 5000000;
    cudaMallocManaged(&mdCandsGPU,MAX_MD_CAND*sizeof(SDL::MiniDoublet));
    mdGPUCounter = 0;
    for (auto& lowerModulePtr : getLowerModulePtrs())
    {
        // Create mini doublets
        createMiniDoubletsFromLowerModule(lowerModulePtr->detId(), MAX_MD_CAND,algo);
    }
    if(mdGPUCounter < MAX_MD_CAND and mdGPUCounter > 0) //incomplete dudes from the final iteration
    {
        miniDoubletGPUWrapper(algo);
    }
}


void SDL::Event::miniDoubletGPUWrapper(SDL::MDAlgo algo)
{
    int nThreads = 256;
    int nBlocks = (mdGPUCounter % nThreads == 0) ? mdGPUCounter/nThreads : mdGPUCounter/nThreads + 1;
    createMiniDoubletsInGPU <<<nBlocks, nThreads>>> (mdCandsGPU,mdGPUCounter,algo);
    cudaError_t cudaerr = cudaDeviceSynchronize();
    if (cudaerr != cudaSuccess)
    {          
        std::cout<<"kernel launch failed with error : "<<cudaGetErrorString(cudaerr)<<std::endl;    
    }

    for(int i = 0; i < mdGPUCounter; i++)
    {
        auto mdCand = mdCandsGPU[i];
        if(mdCand.passesMiniDoubletAlgo(algo))
        {
            // Count the number of md formed
            SDL::Module& lowerModule = (Module&)((mdCand.lowerHitPtr())->getModule()); 
            incrementNumberOfMiniDoublets(lowerModule);

            if (lowerModule.subdet() == SDL::Module::Barrel)
            {
                addMiniDoubletToEvent(mdCand, lowerModule.detId());
            }
            else
            {
                addMiniDoubletToEvent(mdCand, lowerModule.detId());
            }
        }
    }
    mdGPUCounter = 0; 
   
}

void SDL::Event::createMiniDoubletsFromLowerModule(unsigned int detId, int maxMDCands,SDL::MDAlgo algo)
{
    // Get reference to the lower Module
    Module& lowerModule = *getModule(detId);

    // Get reference to the upper Module
    Module& upperModule = *getModule(lowerModule.partnerDetId());
    // Double nested loops
    // Loop over lower module hits
    //Number hardcoded from occupancy plots
    for (auto& lowerHitPtr : lowerModule.getHitPtrs())
    {
        // Get reference to lower Hit
        SDL::Hit& lowerHit = *lowerHitPtr;

        // Loop over upper module hits
        for (auto& upperHitPtr : upperModule.getHitPtrs())
        {

            // Get reference to upper Hit
            SDL::Hit& upperHit = *upperHitPtr;

            // Create a mini-doublet candidate
            SDL::MiniDoublet mdCand(lowerHitPtr, upperHitPtr);
            if(lowerModule.moduleType() == SDL::Module::PS and upperModule.moduleLayerType() == SDL::Module::Strip)
            {
                mdCand.setDrDz(tiltedGeometry.getDrDz(upperModule.detId())); 
            }
            else
            {
                mdCand.setDrDz(tiltedGeometry.getDrDz(lowerModule.detId()));

            }
            if(lowerModule.subdet() == SDL::Module::Endcap)
            {
                if(lowerModule.moduleType() == SDL::Module::PS and upperModule.moduleLayerType() == SDL::Module::Strip)
                {
                    mdCand.setLowerModuleSlope(SDL::endcapGeometry.getSlopeLower(upperModule.detId()));
                }
                else
                {
                    mdCand.setLowerModuleSlope(SDL::endcapGeometry.getSlopeLower(lowerModule.detId()));
                }
            }
            else
            {
                //FIXME: Might need some jugaad for nonexistent det Ids
                if(lowerModule.moduleType() == SDL::Module::PS and upperModule.moduleLayerType() == SDL::Module::Strip)
                {
                    mdCand.setLowerModuleSlope(SDL::tiltedGeometry.getSlope(upperModule.detId()));
                }
                else
                {
                    mdCand.setLowerModuleSlope(SDL::tiltedGeometry.getSlope(lowerModule.detId()));
                }
            }
//	        memcpy(&mdCandsGPU[mdGPUCounter],&mdCand,sizeof(SDL::MiniDoublet));
                mdCandsGPU[mdGPUCounter] = mdCand;
            mdGPUCounter++;

            if(mdGPUCounter == maxMDCands)
            {
                miniDoubletGPUWrapper(algo);
            }
            // Count the number of mdCand considered
            incrementNumberOfMiniDoubletCandidates(lowerModule);
        }
    }       
    // Run mini-doublet algorithm on mdCand (mini-doublet candidate)
           //after running MD algo'
}

// Multiplicity of mini-doublets
unsigned int SDL::Event::getNumberOfHits() { return hits_.size(); }

// Multiplicity of mini-doublets
unsigned int SDL::Event::getNumberOfMiniDoublets() { return miniDoublets_.size(); }


// Multiplicity of mini-doublet candidates considered in this event
unsigned int SDL::Event::getNumberOfMiniDoubletCandidates() { unsigned int n = 0; for (unsigned int i = 0; i < 6; ++i) {n += n_miniDoublet_candidates_by_layer_barrel_[i];} for (unsigned int i = 0; i < 5; ++i) {n += n_miniDoublet_candidates_by_layer_endcap_[i];} return n; }

// Multiplicity of mini-doublet formed in this event
unsigned int SDL::Event::getNumberOfMiniDoubletsByLayerBarrel(unsigned int ilayer) { return n_miniDoublet_by_layer_barrel_[ilayer]; }

// Multiplicity of mini-doublet formed in this event
unsigned int SDL::Event::getNumberOfMiniDoubletsByLayerEndcap(unsigned int ilayer) { return n_miniDoublet_by_layer_endcap_[ilayer]; }

// Multiplicity of hits in this event
void SDL::Event::incrementNumberOfHits(SDL::Module& module)
{
    int layer = module.layer();
    int isbarrel = (module.subdet() == SDL::Module::Barrel);

    // Only count hits in lower module
    if (not module.isLower())
    {
        if (isbarrel)
            n_hits_by_layer_barrel_upper_[layer-1]++;
        else
            n_hits_by_layer_endcap_upper_[layer-1]++;
    }
    else
    {
        if (isbarrel)
            n_hits_by_layer_barrel_[layer-1]++;
        else
            n_hits_by_layer_endcap_[layer-1]++;
    }
}

// Multiplicity of mini-doublet candidates considered in this event
void SDL::Event::incrementNumberOfMiniDoubletCandidates(SDL::Module& module)
{
    int layer = module.layer();
    int isbarrel = (module.subdet() == SDL::Module::Barrel);
    if (isbarrel)
        n_miniDoublet_candidates_by_layer_barrel_[layer-1]++;
    else
        n_miniDoublet_candidates_by_layer_endcap_[layer-1]++;
}

// Multiplicity of mini-doublet formed in this event
void SDL::Event::incrementNumberOfMiniDoublets(SDL::Module& module)
{
    int layer = module.layer();
    int isbarrel = (module.subdet() == SDL::Module::Barrel);
    if (isbarrel)
        n_miniDoublet_by_layer_barrel_[layer-1]++;
    else
        n_miniDoublet_by_layer_endcap_[layer-1]++;
}

void SDL::Event::initMDsInGPU()
{
    const int MD_MAX = 60000;
    cudaMallocManaged(&mdsInGPU,MD_MAX * sizeof(SDL::MiniDoublet));
}

void SDL::Event::addMiniDoubletToEvent(SDL::MiniDoublet md, unsigned int detId)
{
    static int counter = 0;
    if(counter == 0)
    {
        initMDsInGPU();
    }
    // Add to global list of mini doublets, where it will hold the object's instance

    // And get the module (if not exists, then create), and add the address to Module.hits_
    //construct a cudaMallocManaged object and send that in, so that we won't have issues in the GPU
    mdsInGPU[counter] = md;
    getModule(detId)->addMiniDoublet(&mdsInGPU[counter]);
    miniDoublets_.push_back(mdsInGPU[counter]);
    // And get the layer
    counter++;
}



namespace SDL
{
    std::ostream& operator<<(std::ostream& out, const Event& event)
    {

        out << "" << std::endl;
        out << "==============" << std::endl;
        out << "Printing Event" << std::endl;
        out << "==============" << std::endl;
        out << "" << std::endl;

        for (auto& modulePtr : event.modulePtrs_)
        {
            out << modulePtr;
        }


        return out;
    }

    std::ostream& operator<<(std::ostream& out, const Event* event)
    {
        out << *event;
        return out;
    }

}

