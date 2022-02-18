/*

Copyright (c) 2005-2022, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "OffLatticeSimulationWithMonoclonalStoppingEvent.hpp"
#include "Debug.hpp"

bool OffLatticeSimulationWithMonoclonalStoppingEvent::StoppingEventHasOccurred()
{
    if(SimulationTime::Instance()->GetTime() > mStartTime)
    {
        std::set<unsigned> clones;
        for (typename AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
            cell_iter != mrCellPopulation.End();
            ++cell_iter)
        {
            unsigned ancestor = cell_iter->GetAncestor();

            clones.insert(ancestor);
        }

        unsigned num_clones = clones.size();
    
        if(num_clones == 1 )
        {
            PRINT_2_VARIABLES(SimulationTime::Instance()->GetTime(),num_clones);
            return true;
        }
        else
        {
            return false;
        }
    }
    else 
    {
      return false;
    }
}

OffLatticeSimulationWithMonoclonalStoppingEvent::OffLatticeSimulationWithMonoclonalStoppingEvent(
            AbstractCellPopulation<2>& rCellPopulation,
            bool deleteCellPopulationInDestructor,
            bool initialiseCells)
    : OffLatticeSimulation<2>(rCellPopulation,
    deleteCellPopulationInDestructor,
    initialiseCells)
{
    mStartTime = 10;
}

void OffLatticeSimulationWithMonoclonalStoppingEvent::SetStartTime(double startTime)
{
    mStartTime = startTime;
}



// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(OffLatticeSimulationWithMonoclonalStoppingEvent)
