#ifndef Study_h
#define Study_h

#include <vector>
#include <tuple>

#include "SDL/Event.cuh"

class Study
{

public:

    virtual void bookStudy();
    virtual void doStudy(SDL::Event& recoevent, std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents);

};


#endif
