#include "MiniDoublet.cuh"

#define SDL_INF 123456789


SDL::MiniDoublet::MiniDoublet()
{
}

SDL::MiniDoublet::~MiniDoublet()
{
}

SDL::MiniDoublet::MiniDoublet(const MiniDoublet& md): lowerHitPtr_(md.lowerHitPtr()), upperHitPtr_(md.upperHitPtr())
                                                      ,passAlgo_(md.getPassAlgo())
                                                      ,lowerShiftedHit_(md.getLowerShiftedHit())
                                                      ,upperShiftedHit_(md.getUpperShiftedHit())
                                                      ,dz_(md.getDz())
                                                      ,shiftedDz_(md.getShiftedDz())
                                                      ,dphi_(md.getDeltaPhi())
                                                      ,dphi_noshift_(md.getDeltaPhiNoShift())
                                                      ,dphichange_(md.getDeltaPhiChange())
                                                      ,dphichange_noshift_(md.getDeltaPhiChangeNoShift())
                                                      ,miniCut_(md.getMiniCut())
{
    setAnchorHit();
}

SDL::MiniDoublet::MiniDoublet(SDL::Hit* lowerHitPtr, SDL::Hit* upperHitPtr) : lowerHitPtr_(lowerHitPtr), upperHitPtr_(upperHitPtr)
                                                      ,passAlgo_(0)
                                                      ,dz_(0)
                                                      ,shiftedDz_(0)
                                                      ,dphi_(0)
                                                      ,dphi_noshift_(0)
                                                      ,dphichange_(0)
                                                      ,dphichange_noshift_(0)
{
    setAnchorHit();
}


void SDL::MiniDoublet::setAnchorHit()
{
    const SDL::Module& lowerModule = lowerHitPtr()->getModule();

    // Assign anchor hit pointers based on their hit type
    if (lowerModule.moduleType() == SDL::Module::PS)
    {
        if (lowerModule.moduleLayerType() == SDL::Module::Pixel)
        {
            anchorHitPtr_ = lowerHitPtr();
        }
        else
        {
            anchorHitPtr_ = upperHitPtr();
        }
    }
    else
    {
        anchorHitPtr_ = lowerHitPtr();
    }

}

void SDL::MiniDoublet::setLowerModuleSlope(float val)
{
    slopeForHitShifting_ = val;
}

SDL::Hit* SDL::MiniDoublet::lowerHitPtr() const
{
    return lowerHitPtr_;
}

SDL::Hit* SDL::MiniDoublet::upperHitPtr() const
{
    return upperHitPtr_;
}

SDL::Hit* SDL::MiniDoublet::anchorHitPtr() const
{
    return anchorHitPtr_;
}

const int& SDL::MiniDoublet::getPassAlgo() const
{
    return passAlgo_;
}

const SDL::Hit& SDL::MiniDoublet::getLowerShiftedHit() const
{
    return lowerShiftedHit_;
}

const SDL::Hit& SDL::MiniDoublet::getUpperShiftedHit() const
{
    return upperShiftedHit_;
}

const float& SDL::MiniDoublet::getDz() const
{
    return dz_;
}

const float& SDL::MiniDoublet::getShiftedDz() const
{
    return shiftedDz_;
}

const float& SDL::MiniDoublet::getDeltaPhi() const
{
    return dphi_;
}

const float& SDL::MiniDoublet::getDeltaPhiChange() const
{
    return dphichange_;
}

const float& SDL::MiniDoublet::getDeltaPhiNoShift() const
{
    return dphi_noshift_;
}

const float& SDL::MiniDoublet::getDeltaPhiChangeNoShift() const
{
    return dphichange_noshift_;
}


void SDL::MiniDoublet::setLowerShiftedHit(float x, float y, float z, int idx)
{
    lowerShiftedHit_.setXYZ(x, y, z);
    lowerShiftedHit_.setIdx(idx);
}

void SDL::MiniDoublet::setUpperShiftedHit(float x, float y, float z, int idx)
{
    upperShiftedHit_.setXYZ(x, y, z);
    upperShiftedHit_.setIdx(idx);
}

void SDL::MiniDoublet::setDz(float dz)
{
    dz_ = dz;
}

void SDL::MiniDoublet::setDrDz(float drdz)
{
    drdz_ = drdz;
}

void SDL::MiniDoublet::setShiftedDz(float shiftedDz)
{
    shiftedDz_ = shiftedDz;
}

void SDL::MiniDoublet::setDeltaPhi(float dphi)
{
    dphi_ = dphi;
}

void SDL::MiniDoublet::setDeltaPhiChange(float dphichange)
{
    dphichange_ = dphichange;
}

void SDL::MiniDoublet::setDeltaPhiNoShift(float dphi)
{
    dphi_noshift_ = dphi;
}

void SDL::MiniDoublet::setDeltaPhiChangeNoShift(float dphichange)
{
    dphichange_noshift_ = dphichange;
}

void SDL::MiniDoublet::setMiniCut(float var)
{
    miniCut_ = var;
}

const float& SDL::MiniDoublet::getMiniCut() const
{
    return miniCut_;
}
bool SDL::MiniDoublet::passesMiniDoubletAlgo(SDL::MDAlgo algo) const
{
    // Each algorithm is an enum shift it by its value and check against the flag
    return passAlgo_ & (1 << algo);
}

void SDL::MiniDoublet::runMiniDoubletAlgo(SDL::MDAlgo algo, SDL::LogLevel logLevel)
{
    if (algo == SDL::AllComb_MDAlgo)
    {
        runMiniDoubletAllCombAlgo();
    }
    else if (algo == SDL::Default_MDAlgo)
    {
        runMiniDoubletDefaultAlgo(logLevel);
    }
    else
    {
        printf("Warning: Unrecognized mini-doublet algorithm! %d\n",algo);
        return;
    }
}

void SDL::MiniDoublet::runMiniDoubletAllCombAlgo()
{
    passAlgo_ |= (1 << SDL::AllComb_MDAlgo);
}

void SDL::MiniDoublet::runMiniDoubletDefaultAlgo(SDL::LogLevel logLevel)
{
    // Retreived the lower module object
      const SDL::Module& lowerModule = lowerHitPtr_->getModule();

    if (lowerModule.subdet() == SDL::Module::Barrel)
    {
        runMiniDoubletDefaultAlgoBarrel(logLevel);
    }
    else
    {
        runMiniDoubletDefaultAlgoEndcap(logLevel);
    }
}

void SDL::MiniDoublet::runMiniDoubletDefaultAlgoBarrel(SDL::LogLevel logLevel)
{
    // First get the object that the pointer points to
    const SDL::Hit& lowerHit = (*lowerHitPtr_);
    const SDL::Hit& upperHit = (*upperHitPtr_);

    // Retreived the lower module object
    const SDL::Module& lowerModule = lowerHitPtr_->getModule();

//TODO:Change these into regular arrays
    setMiniCut(-999);

    // There are series of cuts that applies to mini-doublet in a "barrel" region

    // Cut #1: The dz difference
    // Ref to original code: https://github.com/slava77/cms-tkph2-ntuple/blob/184d2325147e6930030d3d1f780136bc2dd29ce6/doubletAnalysis.C#L3067

    setDz(lowerHit.z() - upperHit.z());
    const float& dz = getDz();

    // const float dzCut = lowerModule.moduleType() == SDL::Module::PS ? 10.f : 1.5f; // Could be tighter for PS modules

    //*
    // const float dzCut = 10.f; // Could be tighter for PS modules
    // if (not (std::abs(dz) < dzCut)) // If cut fails continue
    //*

    //*
    const float dzCut = lowerModule.moduleType() == SDL::Module::PS ? 2.f : 10.f;
    // const bool isNotInvertedCrosser = lowerModule.moduleType() == SDL::Module::PS ? true : (lowerHit.z() * dz > 0); // Not used as this saves very little on combinatorics. but could be something we can add back later
    const float sign = ((dz > 0) - (dz < 0)) * ((lowerHit.z() > 0) - (lowerHit.z() < 0));
    const float invertedcrossercut = (abs(dz) > 2) * sign;
    if (not (abs(dz) < dzCut and invertedcrossercut <= 0)) // Adding inverted crosser rejection
    //*

    {

/*        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << lowerModule << std::endl;
            SDL::cout << "Debug: " << __FUNCTION__ << "()" << std::endl;
            SDL::cout << "upperHit: " << upperHit << std::endl;
            SDL::cout << "lowerHit: " << lowerHit << std::endl;
            SDL::cout << "dz : " << dz << std::endl;
            SDL::cout << "dzCut : " << dzCut << std::endl;
        }*/

        // did not pass default algo
        passAlgo_ &= (0 << SDL::Default_MDAlgo);
        return;
    }
    /*else
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << lowerModule << std::endl;
            SDL::cout << "Debug: " << __FUNCTION__ << "()" << std::endl;
            SDL::cout << "upperHit: " << upperHit << std::endl;
            SDL::cout << "lowerHit: " << lowerHit << std::endl;
            SDL::cout << "dz : " << dz << std::endl;
            SDL::cout << "dzCut : " << dzCut << std::endl;
        }
    }*/

    // Calculate the cut thresholds for the selection
    float miniCut = 0;
    if (lowerModule.moduleLayerType() == SDL::Module::Pixel)
        miniCut = dPhiThreshold(lowerHit, lowerModule);
    else
        miniCut = dPhiThreshold(upperHit, lowerModule);

    // Cut #2: dphi difference
    // Ref to original code: https://github.com/slava77/cms-tkph2-ntuple/blob/184d2325147e6930030d3d1f780136bc2dd29ce6/doubletAnalysis.C#L3085
    float xn = 0, yn = 0 , zn = 0;
    if (lowerModule.side() != SDL::Module::Center) // If barrel and not center it is tilted
    {
        // Shift the hits and calculate new xn, yn position
        float shiftedCoords[3];
        shiftStripHits(lowerHit, upperHit, lowerModule, shiftedCoords, logLevel);
        xn = shiftedCoords[0];
        yn = shiftedCoords[1];
        zn = shiftedCoords[2];

        // Lower or the upper hit needs to be modified depending on which one was actually shifted
        if (lowerModule.moduleLayerType() == SDL::Module::Pixel)
        {
            // SDL::Hit upperHitMod(upperHit);
            // upperHitMod.setXYZ(xn, yn, upperHit.z());
            // setDeltaPhi(lowerHit.deltaPhi(upperHitMod));
            setUpperShiftedHit(xn, yn, upperHit.z());
            setDeltaPhi(lowerHit.deltaPhi(getUpperShiftedHit()));
            setDeltaPhiNoShift(lowerHit.deltaPhi(upperHit));
        }
        else
        {
            // SDL::Hit lowerHitMod(lowerHit);
            // lowerHitMod.setXYZ(xn, yn, lowerHit.z());
            // setDeltaPhi(lowerHitMod.deltaPhi(upperHit));
            setLowerShiftedHit(xn, yn, lowerHit.z());
            setDeltaPhi(getLowerShiftedHit().deltaPhi(upperHit));
            setDeltaPhiNoShift(lowerHit.deltaPhi(upperHit));
        }
    }
    else
    {
        setDeltaPhi(lowerHit.deltaPhi(upperHit));
        setDeltaPhiNoShift(lowerHit.deltaPhi(upperHit));
    }

    setMiniCut(miniCut);

    if (not (abs(getDeltaPhi()) < miniCut)) // If cut fails continue
    {
        /*if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << lowerModule << std::endl;
            SDL::cout << "Debug: " << __FUNCTION__ << "()" << std::endl;
            SDL::cout << "upperHit: " << upperHit << std::endl;
            SDL::cout << "lowerHit: " << lowerHit << std::endl;
            SDL::cout << "fabsdPhi : " << getDeltaPhi() << std::endl;
            SDL::cout << "miniCut : " << miniCut << std::endl;
        }*/

        // did not pass default algo
        passAlgo_ &= (0 << SDL::Default_MDAlgo);
        return;
    }
    /*else
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << lowerModule << std::endl;
            SDL::cout << "Debug: " << __FUNCTION__ << "()" << std::endl;
            SDL::cout << "upperHit: " << upperHit << std::endl;
            SDL::cout << "lowerHit: " << lowerHit << std::endl;
            SDL::cout << "fabsdPhi : " << getDeltaPhi() << std::endl;
            SDL::cout << "miniCut : " << miniCut << std::endl;
        }
    }*/

    // Cut #3: The dphi change going from lower Hit to upper Hit
    // Ref to original code: https://github.com/slava77/cms-tkph2-ntuple/blob/184d2325147e6930030d3d1f780136bc2dd29ce6/doubletAnalysis.C#L3076
    if (lowerModule.side() != SDL::Module::Center)
    {
        // When it is tilted, use the new shifted positions
        if (lowerModule.moduleLayerType() == SDL::Module::Pixel)
        {
            // SDL::Hit upperHitMod(upperHit);
            // upperHitMod.setXYZ(xn, yn, upperHit.z());
            // dPhi Change should be calculated so that the upper hit has higher rt.
            // In principle, this kind of check rt_lower < rt_upper should not be necessary because the hit shifting should have taken care of this.
            // (i.e. the strip hit is shifted to be aligned in the line of sight from interaction point to pixel hit of PS module guaranteeing rt ordering)
            // But I still placed this check for safety. (TODO: After cheking explicitly if not needed remove later?)
            // setDeltaPhiChange(lowerHit.rt() < upperHitMod.rt() ? lowerHit.deltaPhiChange(upperHitMod) : upperHitMod.deltaPhiChange(lowerHit));
            setDeltaPhiChange(((lowerHit.rt() < getUpperShiftedHit().rt()) ? lowerHit.deltaPhiChange(getUpperShiftedHit()) : getUpperShiftedHit().deltaPhiChange(lowerHit)));
            setDeltaPhiChangeNoShift(((lowerHit.rt() < upperHit.rt()) ? lowerHit.deltaPhiChange(upperHit) : upperHit.deltaPhiChange(lowerHit)));
        }
        else
        {
            // SDL::Hit lowerHitMod(lowerHit);
            // lowerHitMod.setXYZ(xn, yn, lowerHit.z());
            // dPhi Change should be calculated so that the upper hit has higher rt.
            // In principle, this kind of check rt_lower < rt_upper should not be necessary because the hit shifting should have taken care of this.
            // (i.e. the strip hit is shifted to be aligned in the line of sight from interaction point to pixel hit of PS module guaranteeing rt ordering)
            // But I still placed this check for safety. (TODO: After cheking explicitly if not needed remove later?)
            // setDeltaPhiChange(lowerHitMod.rt() < upperHit.rt() ? lowerHitMod.deltaPhiChange(upperHit) : upperHit.deltaPhiChange(lowerHitMod));
            setDeltaPhiChange(((getLowerShiftedHit().rt() < upperHit.rt()) ? getLowerShiftedHit().deltaPhiChange(upperHit) : upperHit.deltaPhiChange(getLowerShiftedHit())));
            setDeltaPhiChangeNoShift(((lowerHit.rt() < upperHit.rt()) ? lowerHit.deltaPhiChange(upperHit) : upperHit.deltaPhiChange(lowerHit)));
        }
    }
    else
    {
        // When it is flat lying module, whichever is the lowerSide will always have rt lower
        setDeltaPhiChange(lowerHit.deltaPhiChange(upperHit));
        setDeltaPhiChangeNoShift(lowerHit.deltaPhiChange(upperHit));
    }

    if (not (std::abs(getDeltaPhiChange()) < miniCut)) // If cut fails continue
    {
        /*if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << lowerModule << std::endl;
            SDL::cout << "Debug: " << __FUNCTION__ << "()" << std::endl;
            SDL::cout << "upperHit: " << upperHit << std::endl;
            SDL::cout << "lowerHit: " << lowerHit << std::endl;
            SDL::cout << "fabsdPhiChange : " << getDeltaPhiChange() << std::endl;
            SDL::cout << "miniCut : " << miniCut << std::endl;
        }*/

        // did not pass default algo
        passAlgo_ &= (0 << SDL::Default_MDAlgo);
        return;
    }
    /*else
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << lowerModule << std::endl;
            SDL::cout << "Debug: " << __FUNCTION__ << "()" << std::endl;
            SDL::cout << "upperHit: " << upperHit << std::endl;
            SDL::cout << "lowerHit: " << lowerHit << std::endl;
            SDL::cout << "fabsdPhiChange : " << getDeltaPhiChange() << std::endl;
            SDL::cout << "miniCut : " << miniCut << std::endl;
        }
    }*/

    // If all cut passed this pair is good, and make and add the mini-doublet
    passAlgo_ |= (1 << SDL::Default_MDAlgo);
    return;

}

void SDL::MiniDoublet::runMiniDoubletDefaultAlgoEndcap(SDL::LogLevel logLevel)
{
    // First get the object that the pointer points to
    const SDL::Hit& lowerHit = (*lowerHitPtr_);
    const SDL::Hit& upperHit = (*upperHitPtr_);

    // Retreived the lower module object
    const SDL::Module& lowerModule = lowerHitPtr_->getModule();

    setMiniCut(-999);

    // There are series of cuts that applies to mini-doublet in a "endcap" region

    // Cut #1 : dz cut. The dz difference can't be larger than 1cm. (max separation is 4mm for modules in the endcap)
    // Ref to original code: https://github.com/slava77/cms-tkph2-ntuple/blob/184d2325147e6930030d3d1f780136bc2dd29ce6/doubletAnalysis.C#L3093
    // For PS module in case when it is tilted a different dz (after the strip hit shift) is calculated later.
    // This is because the 10.f cut is meant more for sanity check (most will pass this cut anyway) (TODO: Maybe revisit this cut later?)

    setDz(lowerHit.z() - upperHit.z());
    float dz = getDz(); // Not const since later it might change depending on the type of module

    const float dzCut = ((lowerModule.side() == SDL::Module::Endcap) ?  1.f : 10.f);
    if (not (abs(dz) < dzCut)) // If cut fails continue
    {
/*        if (logLevel >= SDL::Log_Debug2)
        {
            SDL::cout << lowerModule << std::endl;
            SDL::cout << "Debug: " << __FUNCTION__ << "()" << std::endl;
            SDL::cout << "upperHit: " << upperHit << std::endl;
            SDL::cout << "lowerHit: " << lowerHit << std::endl;
            SDL::cout << "dz : " << dz << std::endl;
            SDL::cout << "dzCut : " << dzCut << std::endl;
        }*/

        // did not pass default algo
        passAlgo_ &= (0 << SDL::Default_MDAlgo);
        return;
    }
    /*else
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << lowerModule << std::endl;
            SDL::cout << "Debug: " << __FUNCTION__ << "()" << std::endl;
            SDL::cout << "upperHit: " << upperHit << std::endl;
            SDL::cout << "lowerHit: " << lowerHit << std::endl;
            SDL::cout << "dz : " << dz << std::endl;
            SDL::cout << "dzCut : " << dzCut << std::endl;
        }
    }*/

    // Cut #2 : drt cut. The dz difference can't be larger than 1cm. (max separation is 4mm for modules in the endcap)
    // Ref to original code: https://github.com/slava77/cms-tkph2-ntuple/blob/184d2325147e6930030d3d1f780136bc2dd29ce6/doubletAnalysis.C#L3100
    const float drtCut = lowerModule.moduleType() == SDL::Module::PS ? 2.f : 10.f;
    float drt = abs(lowerHit.rt() - upperHit.rt());
    if (not (drt < drtCut)) // If cut fails continue
    {
       /* if (logLevel >= SDL::Log_Debug2)
        {
            SDL::cout << lowerModule << std::endl;
            SDL::cout << "Debug: " << __FUNCTION__ << "()" << std::endl;
            SDL::cout << "upperHit: " << upperHit << std::endl;
            SDL::cout << "lowerHit: " << lowerHit << std::endl;
            SDL::cout << "drt : " << drt << std::endl;
            SDL::cout << "drtCut : " << drtCut << std::endl;
        }*/

        // did not pass default algo
        passAlgo_ &= (0 << SDL::Default_MDAlgo);
        return;
    }
    /*else
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << lowerModule << std::endl;
            SDL::cout << "Debug: " << __FUNCTION__ << "()" << std::endl;
            SDL::cout << "upperHit: " << upperHit << std::endl;
            SDL::cout << "lowerHit: " << lowerHit << std::endl;
            SDL::cout << "drt : " << drt << std::endl;
            SDL::cout << "drtCut : " << drtCut << std::endl;
        }
    }*/

    // Calculate the cut thresholds for the selection
    
    // Cut #3: dphi difference
    // Ref to original code: https://github.com/slava77/cms-tkph2-ntuple/blob/184d2325147e6930030d3d1f780136bc2dd29ce6/doubletAnalysis.C#L3111
    // // Old comments ----
    // // Old comments Slava, 6:17 PM
    // // Old comments here for the code to work you would need to slide (line extrapolate) the lower or the upper  hit along the strip direction to the radius of the other
    // // Old comments you'll get it almost right by assuming radial strips and just add the d_rt*(cosPhi, sinPhi)
    // // Old comments ----
    // // Old comments The algorithm assumed that the radial position is ~close according to Slava.
    // // Old comments However, for PS modules, it is not the case.
    // // Old comments So we'd have to move the hits to be in same position as the other.
    // // Old comments We'll move the pixel along the radial direction (assuming the radial direction is more or less same as the strip direction)
    // ----
    // The new scheme shifts strip hits to be "aligned" along the line of sight from interaction point to the pixel hit (if it is PS modules)
    float xn = 0, yn = 0, zn = 0;
    // Shift the hits and calculate new xn, yn position
    float shiftedCoords[3];
    shiftStripHits(lowerHit, upperHit, lowerModule, shiftedCoords, logLevel);
    xn = shiftedCoords[0];
    yn = shiftedCoords[1];
    zn = shiftedCoords[2];

    if (lowerModule.moduleType() == SDL::Module::PS)
    {
        // Appropriate lower or upper hit is modified after checking which one was actually shifted
        if (lowerModule.moduleLayerType() == SDL::Module::Pixel)
        {
            // SDL::Hit upperHitMod(upperHit);
            // upperHitMod.setXYZ(xn, yn, upperHit.z());
            // setDeltaPhi(lowerHit.deltaPhi(upperHitMod));
            setUpperShiftedHit(xn, yn, upperHit.z());
            setDeltaPhi(lowerHit.deltaPhi(getUpperShiftedHit()));
            setDeltaPhiNoShift(lowerHit.deltaPhi(upperHit));
        }
        else
        {
            // SDL::Hit lowerHitMod(lowerHit);
            // lowerHitMod.setXYZ(xn, yn, lowerHit.z());
            // setDeltaPhi(lowerHitMod.deltaPhi(upperHit));
            setLowerShiftedHit(xn, yn, lowerHit.z());
            setDeltaPhi(getLowerShiftedHit().deltaPhi(upperHit));
            setDeltaPhiNoShift(lowerHit.deltaPhi(upperHit));
        }
    }
    else
    {
        // SDL::Hit upperHitMod(upperHit);
        // upperHitMod.setXYZ(xn, yn, upperHit.z());
        // setDeltaPhi(lowerHit.deltaPhi(upperHitMod));
        setUpperShiftedHit(xn, yn, upperHit.z());
        setDeltaPhi(lowerHit.deltaPhi(getUpperShiftedHit()));
        setDeltaPhiNoShift(lowerHit.deltaPhi(upperHit));
    }

    // dz needs to change if it is a PS module where the strip hits are shifted in order to properly account for the case when a tilted module falls under "endcap logic"
    // if it was an endcap it will have zero effect
    if (lowerModule.moduleType() == SDL::Module::PS)
    {
        if (lowerModule.moduleLayerType() == SDL::Module::Pixel)
        {
            setShiftedDz(lowerHit.z() - zn);
            dz = getShiftedDz();
        }
        else
        {
            setShiftedDz(upperHit.z() - zn);
            dz = getShiftedDz();
        }
    }

    float miniCut = 0;
    if (lowerModule.moduleLayerType() == SDL::Module::Pixel)
        miniCut = dPhiThreshold(lowerHit, lowerModule, getDeltaPhi(), dz);
    else
        miniCut = dPhiThreshold(upperHit, lowerModule, getDeltaPhi(), dz);

    setMiniCut(miniCut);

    if (not (std::abs(getDeltaPhi()) < miniCut)) // If cut fails continue
    {
       /* if (logLevel >= SDL::Log_Debug2)
        {
            SDL::cout << lowerModule << std::endl;
            SDL::cout << "Debug: " << __FUNCTION__ << "()" << std::endl;
            SDL::cout << "upperHit: " << upperHit << std::endl;
            SDL::cout << "lowerHit: " << lowerHit << std::endl;
            SDL::cout << "fabsdPhi : " << getDeltaPhi() << std::endl;
            SDL::cout << "miniCut : " << miniCut << std::endl;
        }*/

        // did not pass default algo
        passAlgo_ &= (0 << SDL::Default_MDAlgo);
        return;
    }
/*    else
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << lowerModule << std::endl;
            SDL::cout << "Debug: " << __FUNCTION__ << "()" << std::endl;
            SDL::cout << "upperHit: " << upperHit << std::endl;
            SDL::cout << "lowerHit: " << lowerHit << std::endl;
            SDL::cout << "fabsdPhi : " << getDeltaPhi() << std::endl;
            SDL::cout << "miniCut : " << miniCut << std::endl;
        }
    }*/

    // Cut #4: Another cut on the dphi after some modification
    // Ref to original code: https://github.com/slava77/cms-tkph2-ntuple/blob/184d2325147e6930030d3d1f780136bc2dd29ce6/doubletAnalysis.C#L3119-L3124

    
    float dzFrac = abs(dz) / fabs(lowerHit.z());
    setDeltaPhiChange(getDeltaPhi() / dzFrac * (1.f + dzFrac));
    setDeltaPhiChangeNoShift(getDeltaPhiNoShift() / dzFrac * (1.f + dzFrac));
    if (not (abs(getDeltaPhiChange()) < miniCut)) // If cut fails continue
    {
        /*if (logLevel >= SDL::Log_Debug2)
        {
            SDL::cout << lowerModule << std::endl;
            SDL::cout << "Debug: " << __FUNCTION__ << "()" << std::endl;
            SDL::cout << "upperHit: " << upperHit << std::endl;
            SDL::cout << "lowerHit: " << lowerHit << std::endl;
            SDL::cout << "dz : " << dz << std::endl;
            SDL::cout << "dzFrac : " << dzFrac << std::endl;
            SDL::cout << "fabsdPhi : " << std::abs(getDeltaPhi()) << std::endl;
            SDL::cout << "fabsdPhiMod : " << getDeltaPhiChange() << std::endl;
            SDL::cout << "miniCut : " << miniCut << std::endl;
        }*/

        // did not pass default algo
        passAlgo_ &= (0 << SDL::Default_MDAlgo);
        return;
    }
/*  else
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << lowerModule << std::endl;
            SDL::cout << "Debug: " << __FUNCTION__ << "()" << std::endl;
            SDL::cout << "upperHit: " << upperHit << std::endl;
            SDL::cout << "lowerHit: " << lowerHit << std::endl;
            SDL::cout << "dz : " << dz << std::endl;
            SDL::cout << "dzFrac : " << dzFrac << std::endl;
            SDL::cout << "fabsdPhi : " << std::abs(getDeltaPhi()) << std::endl;
            SDL::cout << "fabsdPhiMod : " << getDeltaPhiChange() << std::endl;
            SDL::cout << "miniCut : " << miniCut << std::endl;
        }
    }*/

    // If all cut passed this pair is good, and make and add the mini-doublet
    passAlgo_ |= (1 << SDL::Default_MDAlgo);
    return;
}

bool SDL::MiniDoublet::isIdxMatched(const MiniDoublet& md) const
{
    if (not (lowerHitPtr()->isIdxMatched(*(md.lowerHitPtr()))))
        return false;
    if (not (upperHitPtr()->isIdxMatched(*(md.upperHitPtr()))))
        return false;
    return true;
}

bool SDL::MiniDoublet::isAnchorHitIdxMatched(const MiniDoublet& md) const
{
    if (not anchorHitPtr_->isIdxMatched(*(md.anchorHitPtr())))
        return false;
    return true;
}

namespace SDL
{
    std::ostream& operator<<(std::ostream& out, const MiniDoublet& md)
    {
        out << "dz " << md.getDz() << std::endl;
        out << "shiftedDz " << md.getShiftedDz() << std::endl;
        out << "dphi " << md.getDeltaPhi() << std::endl;
        out << "dphinoshift " << md.getDeltaPhiNoShift() << std::endl;
        out << "dphichange " << md.getDeltaPhiChange() << std::endl;
        out << "dphichangenoshift " << md.getDeltaPhiChangeNoShift() << std::endl;
        out << "ptestimate " << SDL::MathUtil::ptEstimateFromDeltaPhiChangeAndRt(md.getDeltaPhiChange(), md.lowerHitPtr_->rt()) << std::endl;
        out << "ptestimate from noshift " << SDL::MathUtil::ptEstimateFromDeltaPhiChangeAndRt(md.getDeltaPhiChangeNoShift(), md.lowerHitPtr_->rt()) << std::endl;
        out << std::endl;
        out << "Lower Hit " << std::endl;
        out << "------------------------------" << std::endl;
        {
            IndentingOStreambuf indent(out);
            out << "Lower         " << md.lowerHitPtr_ << std::endl;
            out << "Lower Shifted " << md.lowerShiftedHit_ << std::endl;
            out << md.lowerHitPtr_->getModule() << std::endl;
        }
        out << "Upper Hit " << std::endl;
        out << "------------------------------" << std::endl;
        {
            IndentingOStreambuf indent(out);
            out << "Upper         " << md.upperHitPtr_ << std::endl;
            out << "Upper Shifted " << md.upperShiftedHit_ << std::endl;
            out << md.upperHitPtr_->getModule();
        }
        return out;
    }

    std::ostream& operator<<(std::ostream& out, const MiniDoublet* md)
    {
        out << *md;
        return out;
    }
}

float SDL::MiniDoublet::dPhiThreshold(const SDL::Hit& lowerHit, const SDL::Module& module,const float dPhi, const float dz)
{
    // =================================================================
    // Various constants
    // =================================================================
    const float kRinv1GeVf = (2.99792458e-3 * 3.8);
    const float k2Rinv1GeVf = kRinv1GeVf / 2.;
    float ptCut = 1;
    // if (module.layer() == 6 or module.layer() == 5)
    // {
    //     ptCut = 0.96;
    // }
    float sinAlphaMax = 0.95;
    // p2Sim.directionT-r2Sim.directionT smearing around the mean computed with ptSim,rSim
    // (1 sigma based on 95.45% = 2sigma at 2 GeV)
    float miniMulsPtScaleBarrel[] = {0.0052, 0.0038, 0.0034, 0.0034, 0.0032, 0.0034};
    float miniMulsPtScaleEndcap[] =  {0.006, 0.006, 0.006, 0.006, 0.006}; //inter/extra-polated from L11 and L13 both roughly 0.006 [larger R have smaller value by ~50%]
    //mean of the horizontal layer position in y; treat this as R below
    float miniRminMeanBarrel[] = {21.8, 34.6, 49.6, 67.4, 87.6, 106.8}; // TODO: Update this with newest geometry
    float miniRminMeanEndcap[] = {131.4, 156.2, 185.6, 220.3, 261.5};// use z for endcaps // TODO: Update this with newest geometry

    // =================================================================
    // Computing some components that make up the cut threshold
    // =================================================================
    float rt = lowerHit.rt();
    unsigned int iL = module.layer() - 1;
    const float miniSlope = asin(min(rt * k2Rinv1GeVf / ptCut, sinAlphaMax));
    const float rLayNominal = ((module.subdet() == SDL::Module::Barrel) ? miniRminMeanBarrel[iL] : miniRminMeanEndcap[iL]);
    const float miniPVoff = 0.1 / rLayNominal;
    const float miniMuls = ((module.subdet() == SDL::Module::Barrel) ? miniMulsPtScaleBarrel[iL] * 3.f / ptCut : miniMulsPtScaleEndcap[iL] * 3.f / ptCut);
    const bool isTilted = module.subdet() == SDL::Module::Barrel and module.side() != SDL::Module::Center;
    const bool tiltedOT123 = true;
    const float pixelPSZpitch = 0.15;
    const unsigned int detid = ((module.moduleLayerType() == SDL::Module::Pixel) ?  module.partnerDetId() : module.detId());
    const float drdz = isTilted && tiltedOT123 ? drdz_ : 0;
    const float miniTilt = ((isTilted && tiltedOT123) ? 0.5f * pixelPSZpitch * drdz / sqrt(1.f + drdz * drdz) / moduleGapSize(module) : 0);

    // Compute luminous region requirement for endcap
    const float deltaZLum = 15.f;
    const float miniLum = abs(dPhi * deltaZLum/dz); // Balaji's new error
    // const float miniLum = abs(deltaZLum / lowerHit.z()); // Old error


    // =================================================================
    // Return the threshold value
    // =================================================================
    // Following condition is met if the module is central and flatly lying
    if (module.subdet() == SDL::Module::Barrel and module.side() == SDL::Module::Center)
    {
        return miniSlope + sqrt(pow(miniMuls, 2) + pow(miniPVoff, 2));
    }
    // Following condition is met if the module is central and tilted
    else if (module.subdet() == SDL::Module::Barrel and module.side() != SDL::Module::Center) //all types of tilted modules
    {
        return miniSlope + sqrt(pow(miniMuls, 2) + pow(miniPVoff, 2) + pow(miniTilt * miniSlope, 2));
    }
    // If not barrel, it is Endcap
    else
    {
        return miniSlope + sqrt(pow(miniMuls, 2) + pow(miniPVoff, 2) + pow(miniLum, 2));
    }

}

void SDL::MiniDoublet::shiftStripHits(const SDL::Hit& lowerHit, const SDL::Hit& upperHit, const SDL::Module& lowerModule, float* shiftedCoords, SDL::LogLevel logLevel)
{

    // This is the strip shift scheme that is explained in http://uaf-10.t2.ucsd.edu/~phchang/talks/PhilipChang20190607_SDL_Update.pdf (see backup slides)
    // The main feature of this shifting is that the strip hits are shifted to be "aligned" in the line of sight from interaction point to the the pixel hit.
    // (since pixel hit is well defined in 3-d)
    // The strip hit is shifted along the strip detector to be placed in a guessed position where we think they would have actually crossed
    // The size of the radial direction shift due to module separation gap is computed in "radial" size, while the shift is done along the actual strip orientation
    // This means that there may be very very subtle edge effects coming from whether the strip hit is center of the module or the at the edge of the module
    // But this should be relatively minor effect

    // dependent variables for this if statement
    // lowerModule
    // lowerHit
    // upperHit
    // SDL::endcapGeometry
    // SDL::tiltedGeometry

    // Some variables relevant to the function
    float xp; // pixel x (pixel hit x)
    float yp; // pixel y (pixel hit y)
    float xa; // "anchor" x (the anchor position on the strip module plane from pixel hit)
    float ya; // "anchor" y (the anchor position on the strip module plane from pixel hit)
    float xo; // old x (before the strip hit is moved up or down)
    float yo; // old y (before the strip hit is moved up or down)
    float xn; // new x (after the strip hit is moved up or down)
    float yn; // new y (after the strip hit is moved up or down)
    float abszn; // new z in absolute value
    float zn; // new z with the sign (+/-) accounted
    float angleA; // in r-z plane the theta of the pixel hit in polar coordinate is the angleA
    float angleB; // this is the angle of tilted module in r-z plane ("drdz"), for endcap this is 90 degrees
    unsigned int detid; // The detId of the strip module in the PS pair. Needed to access geometry information
    bool isEndcap; // If endcap, drdz = infinity
    const SDL::Hit* pixelHitPtr; // Pointer to the pixel hit
    const SDL::Hit* stripHitPtr; // Pointer to the strip hit
    float moduleSeparation;
    float drprime; // The radial shift size in x-y plane projection
    float drprime_x; // x-component of drprime
    float drprime_y; // y-component of drprime
    float slope; // The slope of the possible strip hits for a given module in x-y plane
    float absArctanSlope;
    float angleM; // the angle M is the angle of rotation of the module in x-y plane if the possible strip hits are along the x-axis, then angleM = 0, and if the possible strip hits are along y-axis angleM = 90 degrees
    float absdzprime; // The distance between the two points after shifting

    // Assign hit pointers based on their hit type
    if (lowerModule.moduleType() == SDL::Module::PS)
    {
        if (lowerModule.moduleLayerType() == SDL::Module::Pixel)
        {
            pixelHitPtr = &lowerHit;
            stripHitPtr = &upperHit;
            detid = lowerModule.partnerDetId(); // partnerDetId returns the partner detId (since lower Module == pixel, get partner ID to access strip one)
        }
        else
        {
            pixelHitPtr = &upperHit;
            stripHitPtr = &lowerHit;
            detid = lowerModule.detId(); // Since the lower module is not pixel, the lower module is the strip
        }
    }
    else // if (lowerModule.moduleType() == SDL::Module::TwoS) // If it is a TwoS module (if this is called likely an endcap module) then anchor the inner hit and shift the outer hit
    {
        pixelHitPtr = &lowerHit; // Even though in this case the "pixelHitPtr" is really just a strip hit, we pretend it is the anchoring pixel hit
        stripHitPtr = &upperHit;
        detid = lowerModule.detId(); // partnerDetId returns the partner detId (since lower Module == pixel, get partner ID to access strip one)
    }

    // If it is endcap some of the math gets simplified (and also computers don't like infinities)
    isEndcap = lowerModule.subdet() == SDL::Module::Endcap;

    // NOTE: TODO: Keep in mind that the sin(atan) function can be simplifed to something like x / sqrt(1 + x^2) and similar for cos
    // I am not sure how slow sin, atan, cos, functions are in c++. If x / sqrt(1 + x^2) are faster change this later to reduce arithmetic computation time

    // The pixel hit is used to compute the angleA which is the theta in polar coordinate
    // angleA = std::atan(pixelHitPtr->rt() / pixelHitPtr->z() + (pixelHitPtr->z() < 0 ? M_PI : 0)); // Shift by pi if the z is negative so that the value of the angleA stays between 0 to pi and not -pi/2 to pi/2
    angleA = fabs(atan(pixelHitPtr->rt() / pixelHitPtr->z())); // Shift by pi if the z is negative so that the value of the angleA stays between 0 to pi and not -pi/2 to pi/2

    // angleB = isEndcap ? M_PI / 2. : -std::atan(tiltedGeometry.getDrDz(detid) * (lowerModule.side() == SDL::Module::PosZ ? -1 : 1)); // The tilt module on the postive z-axis has negative drdz slope in r-z plane and vice versa
    angleB = ((isEndcap) ? M_PI / 2. : atan(drdz_)); // The tilt module on the postive z-axis has negative drdz slope in r-z plane and vice versa
    // https://iopscience.iop.org/article/10.1088/1748-0221/12/02/C02049/pdf says the PS modules have distances of 1.6mm 2.6mm or 4mm
    moduleSeparation = moduleGapSize(lowerModule);

    // Sign flips if the pixel is later layer
    if (lowerModule.moduleType() == SDL::Module::PS and lowerModule.moduleLayerType() != SDL::Module::Pixel)
    {
        moduleSeparation *= -1;
    }

    drprime = (moduleSeparation / sin(angleA + angleB)) * sin(angleA);

    if (lowerModule.subdet() == SDL::Module::Endcap)
    {
        slope = slopeForHitShifting_;//SDL::endcapGeometry.getSlopeLower(detid); // Only need one slope
    }
    if (lowerModule.subdet() == SDL::Module::Barrel)
    {
        slope = slopeForHitShifting_;//SDL::tiltedGeometry.getSlope(detid);
    }

    // Compute arctan of the slope and take care of the slope = infinity case
    absArctanSlope = ((slope != SDL_INF) ? fabs(atan(slope)) : M_PI / 2); // Since C++ can't represent infinity, SDL_INF = 123456789 was used to represent infinity in the data table

    // The pixel hit position
    xp = pixelHitPtr->x();
    yp = pixelHitPtr->y();

    // Depending on which quadrant the pixel hit lies, we define the angleM by shifting them slightly differently
    if (xp > 0 and yp > 0)
    {
        angleM = absArctanSlope;
    }
    else if (xp > 0 and yp < 0)
    {
        angleM = M_PI - absArctanSlope;
    }
    else if (xp < 0 and yp < 0)
    {
        angleM = M_PI + absArctanSlope;
    }
    else // if (xp < 0 and yp > 0)
    {
        angleM = 2 * M_PI - absArctanSlope;
    }

    // Then since the angleM sign is taken care of properly
    drprime_x = drprime * sin(angleM);
    drprime_y = drprime * cos(angleM);

    // The new anchor position is
    xa = xp + drprime_x;
    ya = yp + drprime_y;

    // The original strip hit position
    xo = stripHitPtr->x();
    yo = stripHitPtr->y();

    // Compute the new strip hit position (if the slope vaule is in special condition take care of the exceptions)
    if (slope == SDL_INF) // Special value designated for tilted module when the slope is exactly infinity (module lying along y-axis)
    {
        xn = xa; // New x point is simply where the anchor is
        yn = yo; // No shift in y
    }
    else if (slope == 0)
    {
        xn = xo; // New x point is simply where the anchor is
        yn = ya; // No shift in y
    }
    else
    {
        xn = (slope * xa + (1.f / slope) * xo - ya + yo) / (slope + (1.f / slope)); // new xn
        yn = (xn - xa) * slope + ya; // new yn
    }

    // Computing new Z position
    absdzprime = fabs(moduleSeparation / sin(angleA + angleB) * cos(angleA)); // module separation sign is for shifting in radial direction for z-axis direction take care of the sign later

    // Depending on which one as closer to the interactin point compute the new z wrt to the pixel properly
    if (lowerModule.moduleLayerType() == SDL::Module::Pixel)
    {
        abszn = fabs(pixelHitPtr->z()) + absdzprime;
    }
    else
    {
        abszn = fabs(pixelHitPtr->z()) - absdzprime;
    }

    zn = abszn * ((pixelHitPtr->z() > 0) ? 1 : -1); // Apply the sign of the zn


    shiftedCoords[0] = xn;
    shiftedCoords[1] = yn;
    shiftedCoords[2] = zn;

}


bool SDL::MiniDoublet::isTighterTiltedModules(const SDL::Module& lowerModule)
{
    // The "tighter" tilted modules are the subset of tilted modules that have smaller spacing
    // This is the same as what was previously considered as"isNormalTiltedModules"
    // See Figure 9.1 of https://cds.cern.ch/record/2272264/files/CMS-TDR-014.pdf
    if (
           (lowerModule.subdet() == SDL::Module::Barrel and lowerModule.side() != SDL::Module::Center and lowerModule.layer() == 3)
           or (lowerModule.subdet() == SDL::Module::Barrel and lowerModule.side() == SDL::Module::NegZ and lowerModule.layer() == 2 and lowerModule.rod() > 5)
           or (lowerModule.subdet() == SDL::Module::Barrel and lowerModule.side() == SDL::Module::PosZ and lowerModule.layer() == 2 and lowerModule.rod() < 8)
           or (lowerModule.subdet() == SDL::Module::Barrel and lowerModule.side() == SDL::Module::NegZ and lowerModule.layer() == 1 and lowerModule.rod() > 9)
           or (lowerModule.subdet() == SDL::Module::Barrel and lowerModule.side() == SDL::Module::PosZ and lowerModule.layer() == 1 and lowerModule.rod() < 4)
       )
        return true;
    else
        return false;
}

// The function to determine gap
float SDL::MiniDoublet::moduleGapSize(const Module& lowerModule)
{
    float miniDeltaTilted[] = {0.26, 0.26, 0.26};
    float miniDeltaLooseTilted[] =  {0.4,0.4,0.4};
    //std::array<float, 6> miniDeltaEndcap {0.4, 0.4, 0.4, 0.18, 0.18, 0.18};
    float miniDeltaFlat[] =  {0.26, 0.16, 0.16, 0.18, 0.18, 0.18};
    float miniDeltaEndcap[5][15];

    for (size_t i = 0; i < 5; i++)
    {
        for (size_t j = 0; j < 15; j++)
        {
            if (i == 0 || i == 1)
            {
                if (j < 10)
                {
                    miniDeltaEndcap[i][j] = 0.4;
                }
                else
                {
                    miniDeltaEndcap[i][j] = 0.18;
                }
            }
            else if (i == 2 || i == 3)
            {
                if (j < 8)
                {
                    miniDeltaEndcap[i][j] = 0.4;
                }
                else
                {
                    miniDeltaEndcap[i][j]  = 0.18;
                }
            }
            else
            {
                if (j < 9)
                {
                    miniDeltaEndcap[i][j] = 0.4;
                }
                else
                {
                    miniDeltaEndcap[i][j] = 0.18;
                }
            }
        }
    }

    unsigned int iL = lowerModule.layer() - 1;
    int iR = lowerModule.subdet() == SDL::Module::Endcap ? lowerModule.ring() - 1 : -1;

    float moduleSeparation = 0;

    if (lowerModule.subdet() == SDL::Module::Barrel and lowerModule.side() == SDL::Module::Center)
    {
        moduleSeparation = miniDeltaFlat[iL];
    }
    else if (isTighterTiltedModules(lowerModule))
    {
        moduleSeparation = miniDeltaTilted[iL];
    }
    else if (lowerModule.subdet() == SDL::Module::Endcap)
    {
        moduleSeparation = miniDeltaEndcap[iL][iR];
    }
    else //Loose tilted modules
    {
        moduleSeparation = miniDeltaLooseTilted[iL];
    }

    return moduleSeparation;

}

