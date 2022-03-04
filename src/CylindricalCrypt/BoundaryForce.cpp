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

#include "BoundaryForce.hpp"
#include "MathsCustomFunctions.hpp"
#include "Exception.hpp"

template<unsigned DIM>
BoundaryForce<DIM>::BoundaryForce()
   : AbstractForce<DIM>(),
     mForceStrength(100.0),
     mCutOffHeight(1.0)
{
}

template<unsigned DIM>
BoundaryForce<DIM>::~BoundaryForce()
{
}

template<unsigned DIM>
void BoundaryForce<DIM>::SetForceStrength(double forceStrength)
{
    assert(forceStrength > 0.0);
    mForceStrength = forceStrength;
}

template<unsigned DIM>
double BoundaryForce<DIM>::GetForceStrength()
{
    return mForceStrength;
}

template<unsigned DIM>
void BoundaryForce<DIM>::SetCutOffHeight(double cutOffHeight)
{
    assert(cutOffHeight > 0.0);
    mCutOffHeight = cutOffHeight;
}

template<unsigned DIM>
double BoundaryForce<DIM>::GetCutOffHeight()
{
    return mCutOffHeight;
}

template<unsigned DIM>
void BoundaryForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    if( (dynamic_cast<AbstractCentreBasedCellPopulation<DIM,DIM>*>(&rCellPopulation) == nullptr)
            && (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr) )
    {
        EXCEPTION("BoundaryForce is to be used with CentreBasedCellPopulations or VertexBasedCellPopulations only");
    }

    if (bool(dynamic_cast<AbstractCentreBasedCellPopulation<DIM,DIM>*>(&rCellPopulation)))
    {
        // Helper variable that is a static cast of the cell population
        AbstractCentreBasedCellPopulation<DIM,DIM>* p_cell_population = static_cast<AbstractCentreBasedCellPopulation<DIM,DIM>*>(&rCellPopulation);
        // Iterate over the nodes
        for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = p_cell_population->rGetMesh().GetNodeIteratorBegin();
            node_iter != p_cell_population->rGetMesh().GetNodeIteratorEnd();
            ++node_iter)
        {
            // Get the node location
            unsigned node_index = node_iter->GetIndex();
            double y = node_iter->rGetLocation()[1]; // y-coordinate of node

            // If the node lies below the line y=mCutOffHeight and is not a ghost, then add the boundary force contribution
            if ((y < mCutOffHeight) && !(p_cell_population->IsGhostNode(node_index)))
            {
                c_vector<double, DIM> boundary_force = zero_vector<double>(DIM);
                // boundary_force[1] = mForceStrength*1.0/SmallPow(y+1.0, 2);
                boundary_force[1] = mForceStrength*exp(-5*y);
                node_iter->AddAppliedForceContribution(boundary_force);
            }
        }
    }
    else
    {
        // Helper variable that is a static cast of the cell population
        VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));

        // Iterate over elements
        for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
            elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
            ++elem_iter)
        {
            unsigned elem_index = elem_iter->GetIndex();
            c_vector<double, DIM> centroid_location = p_cell_population->rGetMesh().GetCentroidOfElement(elem_index);
            double y_centroid = centroid_location[1];

            // If the centroid lies below the line y=mCutOffHeight, then add the boundary force contribution to the node forces
            if (y_centroid < mCutOffHeight)
            {
                unsigned num_nodes = elem_iter->GetNumNodes();
                for (unsigned node_local_index = 0; node_local_index < num_nodes; node_local_index++)
                {
                    unsigned node_global_index = elem_iter->GetNodeGlobalIndex(node_local_index);
                    // double y = p_cell_population->GetNode(node_global_index)->rGetLocation()[1]; // y-coordinate of node
                    c_vector<double, DIM> boundary_force = zero_vector<double>(DIM);
                    // boundary_force[1] = mForceStrength*1.0/SmallPow(y_centroid, 2);
                    boundary_force[1] = mForceStrength*exp(-5*y_centroid);
                    p_cell_population->GetNode(node_global_index)->AddAppliedForceContribution(boundary_force);
                }
            }
            // Add some noise to the bottom nodes
            bool add_noise_to_boundary_nodes = false;
            if(y_centroid < 2.0*mCutOffHeight && add_noise_to_boundary_nodes)
            {
                unsigned num_nodes = elem_iter->GetNumNodes();
                for (unsigned node_local_index = 0; node_local_index < num_nodes; node_local_index++)
                {
                    unsigned node_global_index = elem_iter->GetNodeGlobalIndex(node_local_index);
                    Node<DIM>* p_node = p_mesh->GetNode(node_global_index);
                    if(p_node->IsBoundaryNode())
                    {
                        c_vector<double, DIM> r_node = p_node->rGetLocation();
                        // if(r_node[1] < mCutOffHeight)
                        // {
                            // double y = p_cell_population->GetNode(node_global_index)->rGetLocation()[1]; // y-coordinate of node
                            c_vector<double, DIM> boundary_force = zero_vector<double>(DIM);
                            // boundary_force[1] = mForceStrength*1.0/SmallPow(y_centroid, 2);
                            // boundary_force[1] = 0.1*(RandomNumberGenerator::Instance()->ranf())*mForceStrength*exp(-5*y_centroid);
                            boundary_force[1] = 0.1*(RandomNumberGenerator::Instance()->ranf())*mForceStrength*(mCutOffHeight - r_node[1]);
                            p_cell_population->GetNode(node_global_index)->AddAppliedForceContribution(boundary_force);
                        // }
                        
                    }
                }

            }
        }
    }
    
}

template<unsigned DIM>
void BoundaryForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ForceStrength>" << mForceStrength << "</ForceStrength>\n";
    *rParamsFile << "\t\t\t<CutOffHeight>" << mCutOffHeight << "</CutOffHeight>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class BoundaryForce<1>;
template class BoundaryForce<2>;
template class BoundaryForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BoundaryForce)
