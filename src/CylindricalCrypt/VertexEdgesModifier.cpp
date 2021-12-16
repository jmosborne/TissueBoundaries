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

#include "VertexEdgesModifier.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned DIM>
VertexEdgesModifier<DIM>::VertexEdgesModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
VertexEdgesModifier<DIM>::~VertexEdgesModifier()
{
}

template<unsigned DIM>
void VertexEdgesModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    SmoothEdges(rCellPopulation);
}

template<unsigned DIM>
void VertexEdgesModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We also smooth at the start.
     */
    SmoothEdges(rCellPopulation);
}

template<unsigned DIM>
void VertexEdgesModifier<DIM>::SmoothEdges(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("VertexEdgesModifier is to be used with a VertexBasedCellPopulation only");
    }

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));

    /*
    *        |                              |
    *        |                              |
    *        oB                             |
    *       /  \                            |
    *      /    \           ---->           |
    *   Ao|--rm--|oC                        oB
    *   /          \                        /\
    *  /            \                      /  \
    *
    *   rm = max seperation between vertex_A and vertex_C.
    *   if |R_vertex_A - R_vertex_C| <= rm, then we 
    *   delete both vertex_A and Vertex_C and move vertex_B
    *   to their midpoint. 
    */
    double distanceBetweenVerteciesThreshold = 0.075;
    double distanceToCommonVertexThreshold = 0.20;

    for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
        node_iter != p_mesh->GetNodeIteratorEnd();
        ++node_iter)
    {
        if (node_iter->IsBoundaryNode())
        {
            std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
            if(containing_element_indices.size() == 2)
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
                    std::set<unsigned> element_index_set_n = p_neighbour_node->rGetContainingElementIndices();

                    if(p_neighbour_node->IsBoundaryNode() && (element_index_set_n.size()==1) )
                    {
                        c_vector<double, DIM> r_node_iter = node_iter->rGetLocation();
                        c_vector<double, DIM> r_neighbour = p_neighbour_node->rGetLocation();
                        if(norm_2(p_mesh->GetVectorFromAtoB(r_node_iter,r_neighbour)) < distanceToCommonVertexThreshold) 
                        {
                            boundary_neighbours.push_back(neighbour_index);
                        }

                    }
                }

                if(boundary_neighbours.size() == 2)
                {
                    Node<DIM>* p_neighbour_1 = p_mesh->GetNode(boundary_neighbours[0]);
                    Node<DIM>* p_neighbour_2 = p_mesh->GetNode(boundary_neighbours[1]);
                    Node<DIM>* p_node = p_mesh->GetNode(node_index);

                    c_vector<double, DIM> r_neighbour_1 = p_neighbour_1->rGetLocation();
                    c_vector<double, DIM> r_neighbour_2 = p_neighbour_2->rGetLocation();

                    if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_2)) < distanceBetweenVerteciesThreshold)  
                    {
                        // Delete neighbour 1
                        std::set<unsigned> element_index_set_1 = p_neighbour_1->rGetContainingElementIndices();
                        unsigned elem_index_1 = (*element_index_set_1.begin());
                        VertexElement<DIM,DIM>* p_element_1 = p_mesh->GetElement(elem_index_1);
                        p_element_1->DeleteNode(p_element_1->GetNodeLocalIndex(boundary_neighbours[0]));
                        p_mesh->DeleteNodePriorToReMesh(boundary_neighbours[0]);

                        // Delete neighbour 2
                        std::set<unsigned> element_index_set_2 = p_neighbour_2->rGetContainingElementIndices();
                        unsigned elem_index_2 = (*element_index_set_2.begin());
                        VertexElement<DIM,DIM>* p_element_2 = p_mesh->GetElement(elem_index_2);
                        p_element_2->DeleteNode(p_element_2->GetNodeLocalIndex(boundary_neighbours[1]));
                        p_mesh->DeleteNodePriorToReMesh(boundary_neighbours[1]);

                        // Move node to where the other 2 used to be:
                        p_node->rGetModifiableLocation() = r_neighbour_1 + 0.5 * p_mesh->GetVectorFromAtoB(r_neighbour_1, r_neighbour_2);
                    }
                }
            }
        }
    }

    /*
    *        |                           |
    *        |                           |
    *        oB                          |
    *       /|                           oB
    *     Ao |           ---->          /|
    *    /   |                         / |
    *   /    |                        /  |
    *  /     oC                      /   oC
    *
    *   if vertex_A intersects edge defined by vertex_B-vertex_C
    *   then we delete vertex_A and move vertex_B to the midpoint
    *   of vertex_A-vertex_B
    */
    double distanceToEdgeThreshold = 0.05;

    for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
        node_iter != p_mesh->GetNodeIteratorEnd();
        ++node_iter)
    {
        if (node_iter->IsBoundaryNode())
        {
            std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
            if(containing_element_indices.size() == 2)
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
                    std::set<unsigned> element_index_set_n = p_neighbour_node->rGetContainingElementIndices();

                    if(p_neighbour_node->IsBoundaryNode() && (element_index_set_n.size()==1) )
                    {
                        c_vector<double, DIM> r_node_iter = node_iter->rGetLocation();
                        c_vector<double, DIM> r_neighbour = p_neighbour_node->rGetLocation();
                        
                        boundary_neighbours.push_back(neighbour_index);
                    }
                }

                if(boundary_neighbours.size() == 2)
                {
                    Node<DIM>* p_neighbour_1 = p_mesh->GetNode(boundary_neighbours[0]);
                    Node<DIM>* p_neighbour_2 = p_mesh->GetNode(boundary_neighbours[1]);
                    Node<DIM>* p_node = p_mesh->GetNode(node_index);

                    c_vector<double, DIM> r_node = p_node->rGetLocation();
                    c_vector<double, DIM> r_neighbour_1 = p_neighbour_1->rGetLocation();
                    c_vector<double, DIM> r_neighbour_2 = p_neighbour_2->rGetLocation();
                    
                    if(norm_2(p_mesh->GetVectorFromAtoB(r_node, r_neighbour_2)) < norm_2(p_mesh->GetVectorFromAtoB(r_node, r_neighbour_1)) )
                    {
                        // p_neighbour_2 is closer to p_node, so check if p_neighbour_2 intersects the edge p_node-p_neighbour_1

                        double numerator_signed = (r_neighbour_1[0] - r_node[0])*(r_node[1] - r_neighbour_2[1]) - (r_node[0] - r_neighbour_2[0])*(r_neighbour_1[1] - r_node[1]);
                        double numerator = sqrt(numerator_signed*numerator_signed);
                        double denominator = sqrt((r_neighbour_1[0] - r_node[0])*(r_neighbour_1[0] - r_node[0]) + (r_neighbour_1[1] - r_node[1])*(r_neighbour_1[1] - r_node[1]));
                        
                        double distance_p_2_to_edge = numerator/denominator;

                        if(distance_p_2_to_edge < distanceToEdgeThreshold)
                        {
                            std::set<unsigned> element_index_set_2 = p_neighbour_2->rGetContainingElementIndices();
                            unsigned elem_index_2 = (*element_index_set_2.begin());
                            VertexElement<DIM,DIM>* p_element_2 = p_mesh->GetElement(elem_index_2);
                            p_element_2->DeleteNode(p_element_2->GetNodeLocalIndex(boundary_neighbours[1]));
                            p_mesh->DeleteNodePriorToReMesh(boundary_neighbours[1]);

                            p_node->rGetModifiableLocation() = r_neighbour_2;

                        }
                    }
                    
                    else if(norm_2(p_mesh->GetVectorFromAtoB(r_node, r_neighbour_1)) < norm_2(p_mesh->GetVectorFromAtoB(r_node, r_neighbour_2)))
                    {
                        // p_neighbour_1 is closer to p_node, so check if it intersects the edge p_node-p_neighbour_2

                        double numerator_signed = (r_neighbour_2[0] - r_node[0])*(r_node[1] - r_neighbour_1[1]) - (r_node[0] - r_neighbour_1[0])*(r_neighbour_2[1] - r_node[1]);
                        double numerator = sqrt(numerator_signed*numerator_signed);
                        double denominator = sqrt((r_neighbour_2[0] - r_node[0])*(r_neighbour_2[0] - r_node[0]) + (r_neighbour_2[1] - r_node[1])*(r_neighbour_2[1] - r_node[1]));
                        
                        double distance_p_1_to_edge = numerator/denominator;

                        if(distance_p_1_to_edge < distanceToEdgeThreshold)
                        {
                            std::set<unsigned> element_index_set_1 = p_neighbour_1->rGetContainingElementIndices();
                            unsigned elem_index_1 = (*element_index_set_1.begin());
                            VertexElement<DIM,DIM>* p_element_1 = p_mesh->GetElement(elem_index_1);
                            p_element_1->DeleteNode(p_element_1->GetNodeLocalIndex(boundary_neighbours[0]));
                            p_mesh->DeleteNodePriorToReMesh(boundary_neighbours[0]);

                            p_node->rGetModifiableLocation() = r_neighbour_1;
                        }
                    }
                }
            }
        }
    }

    p_mesh->ReMesh();
    
}

template<unsigned DIM>
void VertexEdgesModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{

    // Call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class VertexEdgesModifier<1>;
template class VertexEdgesModifier<2>;
template class VertexEdgesModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexEdgesModifier)

