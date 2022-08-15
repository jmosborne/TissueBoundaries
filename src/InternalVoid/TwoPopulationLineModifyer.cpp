/*

Copyright (c) 2005-2019, University of Oxford.
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

#include "TwoPopulationLineModifyer.hpp"

#include "NodeBasedCellPopulation.hpp"

#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"

#include "VertexBasedCellPopulation.hpp"

#include "CellLabel.hpp"


#include "Debug.hpp"



template<unsigned DIM>
TwoPopulationLineModifyer<DIM>::TwoPopulationLineModifyer()
    : AbstractCellBasedSimulationModifier<DIM>(),
    mOutputDirectory(""),
    mCutoff(1.5)
{
}

template<unsigned DIM>
void TwoPopulationLineModifyer<DIM>::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;

    OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
    out_stream locationFile = output_file_handler.OpenOutputFile("TwoPopulationsLine.dat");
    locationFile->close();

}

template<unsigned DIM>
void TwoPopulationLineModifyer<DIM>::SetCutoff(double cutoff)
{
	mCutoff = cutoff;
}

template<unsigned DIM>
TwoPopulationLineModifyer<DIM>::~TwoPopulationLineModifyer()
{
}

template<unsigned DIM>
void TwoPopulationLineModifyer<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
}

template<unsigned DIM>
void TwoPopulationLineModifyer<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
    if(DIM == 2)
    {
        OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
        out_stream locationFile = output_file_handler.OpenOutputFile("TwoPopulationsLine.dat", std::ios::app);
        SimulationTime* p_time = SimulationTime::Instance();
        *locationFile << p_time->GetTime() << "\t";
        
        if (bool(dynamic_cast<MeshBasedCellPopulation<DIM,DIM>*>(&rCellPopulation)))
        {
            MeshBasedCellPopulation<DIM,DIM>* p_static_cast_cell_population = static_cast<MeshBasedCellPopulation<DIM,DIM>*>(&rCellPopulation);
            MutableMesh<DIM,DIM>* p_mesh = static_cast<MutableMesh<DIM,DIM>*>(&(p_static_cast_cell_population->rGetMesh()));
            // Iterate over all springs and add force contributions
            for (typename MeshBasedCellPopulation<DIM,DIM>::SpringIterator spring_iterator = p_static_cast_cell_population->SpringsBegin();
                spring_iterator != p_static_cast_cell_population->SpringsEnd();
                ++spring_iterator)
            {
                unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
                unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

                CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeA_global_index);
                CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeB_global_index);

                if( (p_cell_A->HasCellProperty<CellLabel>()) && !(p_cell_B->HasCellProperty<CellLabel>()) )
                {
                    c_vector<double, DIM> cell_location_A = rCellPopulation.GetLocationOfCellCentre(p_cell_A);
                    c_vector<double, DIM> cell_location_B = rCellPopulation.GetLocationOfCellCentre(p_cell_B);

                    c_vector<double,DIM> edge = p_mesh->GetVectorFromAtoB(cell_location_A,cell_location_B);
                    
                    if(norm_2(edge) <= mCutoff)
                    {
                        c_vector<double,DIM> interface_location = cell_location_A + 0.5*edge;

                        *locationFile << interface_location[0] << "\t";
                        *locationFile << interface_location[1] << "\t";
                    }
                        
                }
            }
        }
        else if (bool(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation)))
        {
            AbstractCentreBasedCellPopulation<DIM,DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<DIM,DIM>*>(&rCellPopulation);
            NodesOnlyMesh<DIM>* p_mesh = static_cast<NodesOnlyMesh<DIM>*>(&(p_static_cast_cell_population->rGetMesh()));

            std::vector< std::pair<Node<DIM>*, Node<DIM>* > >& r_node_pairs = p_static_cast_cell_population->rGetNodePairs();

            for (typename std::vector< std::pair<Node<DIM>*, Node<DIM>* > >::iterator iter = r_node_pairs.begin();
                iter != r_node_pairs.end();
                iter++)
            {
                std::pair<Node<DIM>*, Node<DIM>* > pair = *iter;

                unsigned node_a_index = pair.first->GetIndex();
                unsigned node_b_index = pair.second->GetIndex();

                CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(node_a_index);
                CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(node_b_index);

                if( (p_cell_A->HasCellProperty<CellLabel>()) && !(p_cell_B->HasCellProperty<CellLabel>()) )
                {
                    c_vector<double, DIM> cell_location_A = rCellPopulation.GetLocationOfCellCentre(p_cell_A);
                    c_vector<double, DIM> cell_location_B = rCellPopulation.GetLocationOfCellCentre(p_cell_B);

                    c_vector<double,DIM> edge = p_mesh->GetVectorFromAtoB(cell_location_A,cell_location_B);

                    if(norm_2(edge) <= mCutoff)
                    {
                        c_vector<double,DIM> interface_location = cell_location_A + 0.5*edge;

                        *locationFile << interface_location[0] << "\t";
                        *locationFile << interface_location[1] << "\t";
                    }
                    
                }

            }
        }
        else if (bool(dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation)))
        {
            VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
            MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));
            
            unsigned num_nodes = p_cell_population->GetNumNodes();

            for (unsigned node_index=0; node_index<num_nodes; node_index++)
            {
                Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
                // Find the indices of the elements owned by this node
                std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

                // Iterate over these elements
                for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
                    iter != containing_elem_indices.end();
                    ++iter)
                {
                    // Get this element, its index and its number of nodes
                    VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
                    unsigned elem_index = p_element->GetIndex();
                    unsigned num_nodes_elem = p_element->GetNumNodes();

                    // Find the local index of this node in this element
                    unsigned local_index = p_element->GetNodeLocalIndex(node_index);

                    unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
                    Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);

                    // Find the indices of the elements owned by each node
                    std::set<unsigned> elements_containing_nodeA = p_this_node->rGetContainingElementIndices();
                    std::set<unsigned> elements_containing_nodeB = p_next_node->rGetContainingElementIndices();

                    // Find common elements
                    std::set<unsigned> shared_elements;
                    std::set_intersection(elements_containing_nodeA.begin(),
                                        elements_containing_nodeA.end(),
                                        elements_containing_nodeB.begin(),
                                        elements_containing_nodeB.end(),
                                        std::inserter(shared_elements, shared_elements.begin()));
                    if (shared_elements.size() > 1)
                    {
                        {
                            // Work out the number of labelled cells: 0,1 or 2
                            unsigned num_labelled_cells = 0;
                            for (std::set<unsigned>::iterator iter = shared_elements.begin();
                                iter != shared_elements.end();
                                ++iter)
                            {
                                unsigned element_index = *(iter);

                                // Get cell associated with this element
                                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(element_index);

                                if(p_cell->HasCellProperty<CellLabel>())
                                {
                                    num_labelled_cells++;
                                }
                            }

                            if (num_labelled_cells == 1)
                            {
                                c_vector<double, DIM> cell_location_A = p_this_node->rGetLocation();
                                c_vector<double, DIM> cell_location_B = p_next_node->rGetLocation();

                                c_vector<double,DIM> edge = p_mesh->GetVectorFromAtoB(cell_location_A,cell_location_B);
                                
                                c_vector<double,DIM> interface_location = cell_location_A + 0.5*edge;

                                *locationFile << interface_location[0] << "\t";
                                *locationFile << interface_location[1] << "\t";
                                
                            }
                            
                        }
                    }


                }

            }   

        }











        *locationFile << "\n";
        locationFile->close();
        
    }

}



template<unsigned DIM>
void TwoPopulationLineModifyer<DIM>::SetupSolve(AbstractCellPopulation<DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateAtEndOfOutputTimeStep(rCellPopulation);


}


template<unsigned DIM>
void TwoPopulationLineModifyer<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class TwoPopulationLineModifyer<1>;
template class TwoPopulationLineModifyer<2>;
template class TwoPopulationLineModifyer<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TwoPopulationLineModifyer)

