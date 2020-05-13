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

SDL::Module* SDL::Event::getModule(unsigned int detId)
{
    // using std::map::emplace
    std::pair<std::map<unsigned int, Module*>::iterator, bool> emplace_result = modulesMapByDetId_.emplace(detId,nullptr);

    // Retreive the module
    auto& inserted_or_existing = (*(emplace_result.first)).second;

    // If new was inserted, then insert to modulePtrs_ pointer list
    if (emplace_result.second) // if true, new was inserted
    {

        // The pointer to be added
        cudaMallocManaged(&inserted_or_existing,sizeof(Module));
        cudaDeviceSynchronize();

        *inserted_or_existing =SDL:: Module(detId);
         Module* module_ptr = inserted_or_existing;
        
        // Add the module pointer to the list of modules
        modulePtrs_.push_back(module_ptr);
        // If the module is lower module then add to list of lower modules
        if (module_ptr->isLower())
            lowerModulePtrs_.push_back(module_ptr);
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

void SDL::Event::addHitToModule(SDL::Hit hit, unsigned int detId)
{
    // Add to global list of hits, where it will hold the object's instance
    SDL::Hit *hitForGPU;
    cudaMallocManaged(&hitForGPU,sizeof(SDL::Hit));
    cudaDeviceSynchronize();

//    hits_.push_back(*hitForGPU);
       // And get the module (if not exists, then create), and add the address to Module.hits_
    //construct a cudaMallocManaged object and send that in, so that we won't have issues in the GPU
    *hitForGPU = hit;
    hitForGPU->setModule(getModule(detId));
    getModule(detId)->addHit(hitForGPU);
    hits_.push_back(*hitForGPU);

    // Count number of hits in the event
    incrementNumberOfHits(*getModule(detId));

    // If the hit is 2S in the endcap then the hit boundary needs to be set
    if (getModule(detId)->subdet() == SDL::Module::Endcap and getModule(detId)->moduleType() == SDL::Module::TwoS)
    {
        SDL::Hit *hit_2s_high_edge;
        SDL::Hit *hit_2s_low_edge;
        cudaMallocManaged(&hit_2s_high_edge,sizeof(SDL::Hit));
        cudaMallocManaged(&hit_2s_low_edge,sizeof(SDL::Hit));
        cudaDeviceSynchronize();
        
        *hit_2s_high_edge = GeometryUtil::stripHighEdgeHit(*hitForGPU);
        *hit_2s_low_edge = SDL::GeometryUtil::stripLowEdgeHit(*hitForGPU);
        hits_2s_edges_.push_back(*hit_2s_high_edge);
        hitForGPU->setHitHighEdgePtr(hit_2s_high_edge);
        hits_2s_edges_.push_back(*hit_2s_low_edge);
        hitForGPU->setHitLowEdgePtr(hit_2s_low_edge);
    }
}

void SDL::Event::createMiniDoublets(MDAlgo algo)
{
    for (auto& lowerModulePtr : getLowerModulePtrs())
    {
        createMiniDoubletsFromLowerModule(lowerModulePtr->detId(), algo);
    }
}

void SDL::Event::createMiniDoubletsFromLowerModule(unsigned int detId, SDL::MDAlgo algo)
{
    // Get reference to the lower Module
    Module& lowerModule = *getModule(detId);

    // Get reference to the upper Module
    Module& upperModule = *getModule(lowerModule.partnerDetId());
    // Double nested loops
    // Loop over lower module hits
    //Number hardcoded from occupancy plots
    SDL::MiniDoublet* mdCandsGPU;
    cudaMallocManaged(&mdCandsGPU,100 * 100 * sizeof(SDL::MiniDoublet));
    cudaDeviceSynchronize();
    int counter = 0;
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
                mdCand.setLowerModuleSlope(SDL::endcapGeometry.getSlopeLower(lowerModule.detId()));
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
	        memcpy(&mdCandsGPU[counter],&mdCand,sizeof(SDL::MiniDoublet));
            counter++;
            incrementNumberOfMiniDoubletCandidates(lowerModule);
        }
    }       
    int nThreads = 32;
    int nBlocks = (counter % nThreads == 0) ? counter/nThreads : counter/nThreads + 1;
    if(counter > 0)
    {
        createMiniDoubletsInGPU <<<nBlocks, nThreads>>> (mdCandsGPU,counter,algo);
        cudaError_t cudaerr = cudaDeviceSynchronize();
        if (cudaerr != cudaSuccess)
        {       
            std::cout<<"kernel launch failed with error : "<<cudaGetErrorString(cudaerr)<<std::endl;    

        }
    }
    for(int i = 0; i < counter; i++)
    {
        auto mdCand = mdCandsGPU[i];
        if(mdCand.passesMiniDoubletAlgo(algo))
        {
            // Count the number of md formed
            incrementNumberOfMiniDoublets(lowerModule);
            addMiniDoubletToEvent(mdCand,detId);
        }
          
    }
    cudaFree(mdCandsGPU);
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

void SDL::Event::addMiniDoubletToEvent(SDL::MiniDoublet md, unsigned int detId)
{
    // Add to global list of mini doublets, where it will hold the object's instance

    // And get the module (if not exists, then create), and add the address to Module.hits_
    //construct a cudaMallocManaged object and send that in, so that we won't have issues in the GPU
    SDL::MiniDoublet *mdForGPU;
    cudaMallocManaged(&mdForGPU,sizeof(SDL::MiniDoublet));
    cudaDeviceSynchronize();
    *mdForGPU = md;
    getModule(detId)->addMiniDoublet(mdForGPU);
    miniDoublets_.push_back(*mdForGPU);
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

