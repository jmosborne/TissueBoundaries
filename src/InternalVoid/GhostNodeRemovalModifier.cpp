/*

Copyright (c) 2005-2021, University of Oxford.
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

#include "GhostNodeRemovalModifier.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"

template<unsigned DIM>
GhostNodeRemovalModifier<DIM>::GhostNodeRemovalModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
    mVolumeThreshold(0.4)
{
}

template<unsigned DIM>
GhostNodeRemovalModifier<DIM>::~GhostNodeRemovalModifier()
{
}

template<unsigned DIM>
void GhostNodeRemovalModifier<DIM>::SetVolumeThreshold(double volumeThreshold)
{
	mVolumeThreshold = volumeThreshold;
}

template<unsigned DIM>
void GhostNodeRemovalModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateGhostNodes(rCellPopulation);
}

template<unsigned DIM>
void GhostNodeRemovalModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
}

template<unsigned DIM>
void GhostNodeRemovalModifier<DIM>::UpdateGhostNodes(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    assert((dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation)));
    
    MeshBasedCellPopulationWithGhostNodes<DIM>* p_static_cast_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

    p_static_cast_cell_population->CreateVoronoiTessellation();

    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = p_static_cast_cell_population->rGetMesh().GetNodeIteratorBegin();
         node_iter != p_static_cast_cell_population->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();

        if (p_static_cast_cell_population->IsGhostNode(node_index))
        {
            double element_volume = p_static_cast_cell_population->GetVoronoiTessellation()->GetVolumeOfElement(node_index);   

            double threshold_element_volume = mVolumeThreshold;

            if (element_volume < threshold_element_volume)
            {
                p_static_cast_cell_population->RemoveGhostNode(node_index);
            }
            else 
            {
                // If no ghost neighbours remove the ghost node
                std::set<unsigned> neighbour_indices = p_static_cast_cell_population->GetNeighbouringNodeIndices(node_index);

                for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                    iter != neighbour_indices.end();)
                {
                    if (!p_static_cast_cell_population->IsGhostNode(*iter))
                    {
                        neighbour_indices.erase(iter++);
                    }
                    else
                    {
                        ++iter;
                    }
                }
                if (neighbour_indices.size()==0)
                {
                    p_static_cast_cell_population->RemoveGhostNode(node_index);
                }
            }
        }
    }

    // This calls ReMesh() and updates all the Cell locations. 
    p_static_cast_cell_population->Update();


}

template<unsigned DIM>
void GhostNodeRemovalModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class GhostNodeRemovalModifier<1>;
template class GhostNodeRemovalModifier<2>;
template class GhostNodeRemovalModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(GhostNodeRemovalModifier)
