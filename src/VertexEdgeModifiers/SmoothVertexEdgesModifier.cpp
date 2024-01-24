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

#include "SmoothVertexEdgesModifier.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "Debug.hpp"

template<unsigned DIM>
SmoothVertexEdgesModifier<DIM>::SmoothVertexEdgesModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
SmoothVertexEdgesModifier<DIM>::~SmoothVertexEdgesModifier()
{
}

template<unsigned DIM>
void SmoothVertexEdgesModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    SmoothEdges(rCellPopulation);
}

template<unsigned DIM>
void SmoothVertexEdgesModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We also smooth at the start.
     */
    SmoothEdges(rCellPopulation);
}

template<unsigned DIM>
void SmoothVertexEdgesModifier<DIM>::SmoothEdges(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("SmoothVertexEdgesModifier is to be used with a VertexBasedCellPopulation only");
    }

    // unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    // std::stringstream time;
    // time << num_timesteps;
    // PRINT_VARIABLE(num_timesteps);

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));

    for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
        node_iter != p_mesh->GetNodeIteratorEnd();
        ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
        if (containing_element_indices.size() == 1)
        {
            // Get this element
            unsigned elem_index = (*containing_element_indices.begin());

            VertexElement<DIM,DIM>* p_element = p_mesh->GetElement(elem_index);

            // Remove node from this element and delete the node
            p_element->DeleteNode(p_element->GetNodeLocalIndex(node_index));
            p_mesh->DeleteNodePriorToReMesh(node_index);
        }
    }
    p_mesh->ReMesh();

    
    
    double distanceBetweenVerteciesThreshold = 0.025;
    double mThetaThreshold = 0.0;

    bool performed_edge_modifier = true;
    while(performed_edge_modifier)
    {
        performed_edge_modifier = false;

        for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
            node_iter != p_mesh->GetNodeIteratorEnd();
            ++node_iter)
        {
            if (node_iter->IsBoundaryNode())
            {
                unsigned node_index = node_iter->GetIndex();
                std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);
                std::vector<unsigned> boundary_neighbours;
                for (std::set<unsigned>::iterator neighbour_iter = node_neighbours.begin();
                    neighbour_iter != node_neighbours.end();
                    ++neighbour_iter)
                {
                    unsigned neighbour_index = *neighbour_iter;
                    Node<DIM>* p_neighbour_node = p_mesh->GetNode(neighbour_index);

                    if(p_neighbour_node->IsBoundaryNode() )
                    {        
                        boundary_neighbours.push_back(neighbour_index);
                    }
                }
                if(boundary_neighbours.size()==2)
                {
                    Node<DIM>* p_neighbour_1 = p_mesh->GetNode(boundary_neighbours[0]);
                    Node<DIM>* p_neighbour_2 = p_mesh->GetNode(boundary_neighbours[1]);
                    Node<DIM>* p_node = p_mesh->GetNode(node_index);

                    c_vector<double, DIM> r_node = p_node->rGetLocation();
                    c_vector<double, DIM> r_neighbour_1 = p_neighbour_1->rGetLocation();
                    c_vector<double, DIM> r_neighbour_2 = p_neighbour_2->rGetLocation();

                    c_vector<double, DIM> r_node_to_1 = p_mesh->GetVectorFromAtoB(r_neighbour_1,r_node);
                    c_vector<double, DIM> r_node_to_2 = p_mesh->GetVectorFromAtoB(r_neighbour_2,r_node);

                    double cos_theta = ((r_node_to_1[0]*r_node_to_2[0] + r_node_to_1[1]*r_node_to_2[1])/(norm_2(r_node_to_1)*norm_2(r_node_to_2)));

                    if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_2)) < distanceBetweenVerteciesThreshold || acos(cos_theta)< mThetaThreshold)  
                    {
                        // TRACE("Doing thins!!");

                        std::set<unsigned> containing_element_0 = p_node->rGetContainingElementIndices();
                        std::set<unsigned> containing_element_1 = p_neighbour_1->rGetContainingElementIndices();
                        std::set<unsigned> containing_element_2 = p_neighbour_2->rGetContainingElementIndices();

                        // Find common elements
                        std::set<unsigned> shared_elements_1;
                        std::set_intersection(containing_element_0.begin(),
                                                containing_element_0.end(),
                                                containing_element_1.begin(),
                                                containing_element_1.end(),
                                                std::inserter(shared_elements_1, shared_elements_1.begin()));

                        std::set<unsigned> shared_elements_2;
                        std::set_intersection(containing_element_0.begin(),
                                                containing_element_0.end(),
                                                containing_element_2.begin(),
                                                containing_element_2.end(),
                                                std::inserter(shared_elements_2, shared_elements_2.begin()));

                        unsigned elem_index_1 = (*shared_elements_1.begin());
                        VertexElement<DIM,DIM>* p_element_1 = p_mesh->GetElement(elem_index_1);
                        p_element_1->DeleteNode(p_element_1->GetNodeLocalIndex(node_index));

                        unsigned elem_index_2 = (*shared_elements_2.begin());
                        VertexElement<DIM,DIM>* p_element_2 = p_mesh->GetElement(elem_index_2);
                        p_element_2->DeleteNode(p_element_2->GetNodeLocalIndex(node_index));
                        p_mesh->DeleteNodePriorToReMesh(node_index);

                        p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);
                        p_neighbour_1->SetAsBoundaryNode(true);
                        
                        // PRINT_VARIABLE(SimulationTime::Instance()->GetTime());
                        // PRINT_VECTOR(r_neighbour_1);
                        
                        performed_edge_modifier = true;
                        break;
                    }
                }
            }
        }
        p_mesh->ReMesh();

    }


    performed_edge_modifier = true;
    double mdistanceBetweenVerteciesThreshold = 0.025;
    while(performed_edge_modifier)
    {
        performed_edge_modifier = false;

        for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
            node_iter != p_mesh->GetNodeIteratorEnd();
            ++node_iter)
        {
            if (node_iter->IsBoundaryNode())
            {
                unsigned node_index = node_iter->GetIndex();
                std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);
                for (std::set<unsigned>::iterator neighbour_iter = node_neighbours.begin();
                    neighbour_iter != node_neighbours.end();
                    ++neighbour_iter)
                {
                    unsigned neighbour_index = *neighbour_iter;
                    Node<DIM>* p_neighbour_node = p_mesh->GetNode(neighbour_index);

                    if(p_neighbour_node->IsBoundaryNode() )
                    {        
                        Node<DIM>* p_node = p_mesh->GetNode(node_index);

                        c_vector<double, DIM> r_node = p_node->rGetLocation();
                        c_vector<double, DIM> r_neighbour_1 = p_neighbour_node->rGetLocation();

                        if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_node)) < mdistanceBetweenVerteciesThreshold )  
                        {
                            // TRACE("Node Merge");

                            p_mesh->PerformNodeMerge(p_neighbour_node,p_node);
                            p_neighbour_node->SetAsBoundaryNode(true);
                        
                            // PRINT_VARIABLE(SimulationTime::Instance()->GetTime());
                            // PRINT_VECTOR(r_neighbour_1);
                        
                            performed_edge_modifier = true;
                            break;
                        }
                    }
                    if(performed_edge_modifier)
                    {
                        break;
                    }
                }
            }
        }
        p_mesh->ReMesh();

    }

    double distanceBetweenVerteciesThreshold_2 = 0.05; //0.075


        if(true)
        {
            for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter_1 = p_mesh->GetNodeIteratorBegin();
                node_iter_1 != p_mesh->GetNodeIteratorEnd();
                ++node_iter_1)
            {
                if (node_iter_1->IsBoundaryNode())
                {
                    std::set<unsigned> containing_element_indices_1 = node_iter_1->rGetContainingElementIndices();
                    
                    if(containing_element_indices_1.size() == 2)
                    {
                        for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter_2 = p_mesh->GetNodeIteratorBegin();
                        node_iter_2 != p_mesh->GetNodeIteratorEnd();
                        ++node_iter_2)
                        {
                            unsigned node_index_1 = node_iter_1->GetIndex();
                            unsigned node_index_2 = node_iter_2->GetIndex();

                            if (node_iter_2->IsBoundaryNode() && node_index_1 != node_index_2 )
                            {
                                std::set<unsigned> containing_element_indices_2 = node_iter_2->rGetContainingElementIndices();
                    
                                if(containing_element_indices_2.size() == 2)
                                {
                                    c_vector<double, DIM> r_node_1 = node_iter_1->rGetLocation();
                                    c_vector<double, DIM> r_node_2 = node_iter_2->rGetLocation();

                                    if(norm_2(p_mesh->GetVectorFromAtoB(r_node_1,r_node_2)) < distanceBetweenVerteciesThreshold_2)
                                    {
                                        // PRINT_VECTOR(r_node_1);
                                        // PRINT_VECTOR(r_node_2);
                                        Node<DIM>* p_node_1 = p_mesh->GetNode(node_index_1);
                                        Node<DIM>* p_node_2 = p_mesh->GetNode(node_index_2);

                                        c_vector<double, DIM> r_node_12 = p_mesh->GetVectorFromAtoB(r_node_1,r_node_2);

                                        p_node_1->rGetModifiableLocation()[0] = r_node_1[0] + 0.05*r_node_12[0];
                                        p_node_1->rGetModifiableLocation()[1] = r_node_1[1] - 0.05*r_node_12[1];

                                        p_node_2->rGetModifiableLocation()[0] = r_node_2[0] - 0.05*r_node_12[0];
                                        p_node_2->rGetModifiableLocation()[1] = r_node_2[1] + 0.05*r_node_12[1];

                                        // p_mesh->PerformNodeMerge(p_node_1,p_node_2);

                                        // p_node_1->SetAsBoundaryNode(true);
                                        performed_edge_modifier = true;
                                        break;

                                    }


                                }

                            }
                        }



                    }
                }
            }
        }
        if(performed_edge_modifier)
        {
            // TRACE("Here");
            p_mesh->ReMesh();
            // TRACE("11_mesh");
            // unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
            // std::stringstream time;
            // time << num_timesteps;
            // PRINT_VARIABLE(num_timesteps);
            // VertexMeshWriter<DIM,DIM> vertexmesh_writer("tmp", "11_mesh", false);
            // vertexmesh_writer.WriteVtkUsingMesh(*p_mesh, time.str());

            // PRINT_VARIABLE(SimulationTime::Instance()->GetTime());
        }
        performed_edge_modifier = false;


    // Might have this sometimes due to boundary node misslabelling... 
        /*  
        *           
        *   --------o---------
        *           |
        *  (Cell)   o   (Cell)
        *           |
        *           |
        *   --------o----------
        *           |
        */
        // i.e. internal nodes...
        // TRACE("Checking Internal Nods");
 
        performed_edge_modifier = true;
        while(performed_edge_modifier)
        {
            performed_edge_modifier = false;
            for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
                node_iter != p_mesh->GetNodeIteratorEnd();
                ++node_iter)
            {

                    std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
                    if(containing_element_indices.size() == 2)
                    {
                        unsigned node_index = node_iter->GetIndex();
                        std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);
                        
                        if(node_neighbours.size()==2)
                        {
                            std::set<unsigned>::const_iterator neigh_it = node_neighbours.begin();
                            unsigned node_neighbour_1 = (*neigh_it);
                            neigh_it++;
                            unsigned node_neighbour_2 = (*neigh_it);

                            Node<DIM>* p_neighbour_node_1 = p_mesh->GetNode(node_neighbour_1);
                            Node<DIM>* p_neighbour_node_2 = p_mesh->GetNode(node_neighbour_2);

                            std::set<unsigned> containing_element_indices_1 = p_neighbour_node_1->rGetContainingElementIndices();
                            std::set<unsigned> containing_element_indices_2 = p_neighbour_node_2->rGetContainingElementIndices();

                            // Find common elements
                            std::set<unsigned> shared_elements;
                            std::set_intersection(containing_element_indices_1.begin(),
                                                containing_element_indices_1.end(),
                                                containing_element_indices_2.begin(),
                                                containing_element_indices_2.end(),
                                                std::inserter(shared_elements, shared_elements.begin()));

                            if(containing_element_indices_1.size() >= 2 && containing_element_indices_2.size() >= 2 && shared_elements.size() >= 2)
                            {
                                // Node<DIM>* p_node = p_mesh->GetNode(node_index);

                                std::set<unsigned>::const_iterator elem_it = containing_element_indices.begin();

                                unsigned elem_index_1 = (*elem_it);
                                VertexElement<DIM,DIM>* p_element_1 = p_mesh->GetElement(elem_index_1);
                                p_element_1->DeleteNode(p_element_1->GetNodeLocalIndex(node_index));

                                elem_it++;
                                unsigned elem_index_2 = (*elem_it);
                                VertexElement<DIM,DIM>* p_element_2 = p_mesh->GetElement(elem_index_2);
                                p_element_2->DeleteNode(p_element_2->GetNodeLocalIndex(node_index));

                                p_mesh->DeleteNodePriorToReMesh(node_index);

                                performed_edge_modifier = true;
                                break;
                            }
                        }
                    }
            }
        }
        

        // Clean any rouge nodes that have incorrectly been labeled a boundary node.
        if(true)
        {
            for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
                node_iter != p_mesh->GetNodeIteratorEnd();
                ++node_iter)
            {
                bool is_boundary = false;
                // bool is_boundary_tmp = false;
                if (node_iter->IsBoundaryNode())
                {
                    std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
                    if(containing_element_indices.size() == 1)
                    {
                        is_boundary = true;
                    }
                    else
                    {
                        unsigned node_index = node_iter->GetIndex();
                        std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);
                        unsigned number_of_boundary_neighbours = 0;

                        for (std::set<unsigned>::iterator neighbour_iter = node_neighbours.begin();
                            neighbour_iter != node_neighbours.end();
                            ++neighbour_iter)
                        {
                            unsigned neighbour_index = *neighbour_iter;
                            Node<DIM>* p_neighbour_node = p_mesh->GetNode(neighbour_index);

                            if(p_neighbour_node->IsBoundaryNode())
                            {
                                number_of_boundary_neighbours++;
                            }
                        }
                        if(number_of_boundary_neighbours >= 2)
                        {
                            is_boundary = true;
                        }
                    }

                    
                }
                else
                {
                    std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
                    if(containing_element_indices.size() == 1)
                    {
                        is_boundary = true;
                    }
                    else if(containing_element_indices.size() > 1)
                    {
                        unsigned node_index = node_iter->GetIndex();
                        std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);

                        for (std::set<unsigned>::iterator neighbour_iter = node_neighbours.begin();
                            neighbour_iter != node_neighbours.end();
                            ++neighbour_iter)
                        {
                            unsigned neighbour_index = *neighbour_iter;
                            Node<DIM>* p_neighbour_node = p_mesh->GetNode(neighbour_index);
                            std::set<unsigned> element_index_set_n = p_neighbour_node->rGetContainingElementIndices();

                            if(p_neighbour_node->IsBoundaryNode() && (element_index_set_n.size()==1) )
                            {
                                is_boundary = true;
                            }
                        }
                    }
                }

                node_iter->SetAsBoundaryNode(is_boundary);
            }
        }
        p_mesh->ReMesh();
       
    
}

template<unsigned DIM>
void SmoothVertexEdgesModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{

    // Call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class SmoothVertexEdgesModifier<1>;
template class SmoothVertexEdgesModifier<2>;
template class SmoothVertexEdgesModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SmoothVertexEdgesModifier)

