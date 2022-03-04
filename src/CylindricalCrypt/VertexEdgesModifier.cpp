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
#include "UblasCustomFunctions.hpp"

#include "VertexMeshWriter.hpp"
#include "MutableMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "Debug.hpp"

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
double VertexEdgesModifier<DIM>::DistanceToEdgeFromPoint( c_vector<double, DIM> P1,  c_vector<double, DIM> P2,  c_vector<double, DIM> P0)
{
    // Line is P1-------P2
    // Point is    P0
    if(DIM == 2)
    {
        double numerator = fabs((P2[0] - P1[0])*(P1[1] - P0[1]) - (P1[0] - P0[0])*(P2[1] - P1[1]));
        double denominator = sqrt((P2[0] - P1[0])*(P2[0] - P1[0]) + (P2[1] - P1[1])*(P2[1] - P1[1]));
                                
        return numerator/denominator;
    }
    else
    {
        return 0.0;
    }
    
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

    bool performed_edge_modifier = false;
    bool performed_edge_modifier_2 = false;

    // PRINT_VARIABLE(SimulationTime::Instance()->GetTime());


    /*
    *        |                              |
    *        |                              |
    *  Cell  oB  Cell                 Cell  |  Cell
    *       /  \                            |
    *      /    \           ---->           |
    *   Ao|--rm--|oC                        oB
    *   /          \                        /\
    *  /            \                      /  \
    *
    *   rm = max seperation between vertex_A and vertex_C -> distanceBetweenVerteciesThreshold.
    *   if |R_vertex_A - R_vertex_C| <= rm, then we 
    *   delete both vertex_A and Vertex_C and move vertex_B
    *   to their midpoint. 
    */
    double distanceBetweenVerteciesThreshold = 0.1; //0.075
    double distanceToCommonVertexThreshold = 0.5; // must be greater than mMaxEdgeLength in VertexBoundaryRefinementModifier.cpp

    if(true)
    {
        for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
            node_iter != p_mesh->GetNodeIteratorEnd();
            ++node_iter)
        {
            if (node_iter->IsBoundaryNode())
            {
                std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
                
                c_vector<double, DIM> r_node = node_iter->rGetLocation();
                // PRINT_VARIABLE(containing_element_indices.size());
                // PRINT_VECTOR(r_node);

                // if(containing_element_indices.size() == 2)
                if(containing_element_indices.size() > 1)
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

                    // PRINT_VECTOR(node_iter->rGetLocation());
                    // PRINT_VARIABLE(boundary_neighbours.size());

                    // If we have two nodes which are closer to the common vertex than thethreshold distanceToCommonVertexThreshold,
                    // we might be able to delete them if they are close enough to eachother.
                    if(boundary_neighbours.size() == 2)
                    {
                        Node<DIM>* p_neighbour_1 = p_mesh->GetNode(boundary_neighbours[0]);
                        Node<DIM>* p_neighbour_2 = p_mesh->GetNode(boundary_neighbours[1]);
                        Node<DIM>* p_node = p_mesh->GetNode(node_index);

                        c_vector<double, DIM> r_neighbour_1 = p_neighbour_1->rGetLocation();
                        c_vector<double, DIM> r_neighbour_2 = p_neighbour_2->rGetLocation();

                        // We have this:
                        /*  
                        *     \ (Cell) /
                        *      \     /
                        *       \   /
                        *        \ /
                        * (Cell)  o  (Cell)
                        *       /   \
                        *      / (v) \
                        * ----o       o -----
                        *
                        */
                        if( containing_element_indices.size() == 3 )
                        {
                            if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_2)) < distanceBetweenVerteciesThreshold)  
                            {
                                p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);
                                p_node->SetAsBoundaryNode(false);

                                performed_edge_modifier = true;
                            }

                        }
                        else
                        {
                            // PRINT_VECTOR( r_neighbour_1 );
                            // PRINT_VECTOR( r_neighbour_2 );
                            // PRINT_VARIABLE( norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_2)) );

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

                                performed_edge_modifier = true;
                                // PRINT_VECTOR(p_node->rGetModifiableLocation());
                                // TRACE("Custom 1.1");

                            }
                        }
                    }

                    /* We might have the weird case as bellow:
                    *                  | (v)/                                  | (v)/      
                    *                  |   /                                   |   /
                    *                  |  o(0)                                 |  o       
                    *                  | /                                     | /
                    *                  |/                                      |/
                    *      (Cell)      o(p)  (Cell)   ---------->      (Cell)  o   (Cell)
                    *                 / |                                      |
                    *                /  |                                      |
                    *               /(v)|                                      |
                    *        ------o    o------                        --------o-------
                    *             (1)  (2)                                  (v)
                    * 
                    *   Here, nodes 0,1,2 can all be boundary nodes which 
                    *   only belong to a single cell(element) each.
                    */
                    else if(boundary_neighbours.size() == 3)
                    {
                        Node<DIM>* p_node = p_mesh->GetNode(node_index);


                        unsigned boundary_neighbour_0 = boundary_neighbours[0];
                        Node<DIM>* p_neighbour_0= p_mesh->GetNode(boundary_neighbour_0);
                        std::set<unsigned> element_index_set_0 = p_neighbour_0->rGetContainingElementIndices();
                        unsigned element_index_0 = *element_index_set_0.begin();
                        c_vector<double, DIM> r_neighbour_0 = p_neighbour_0->rGetLocation();

                        unsigned boundary_neighbour_1 = boundary_neighbours[1];
                        Node<DIM>* p_neighbour_1= p_mesh->GetNode(boundary_neighbour_1);
                        std::set<unsigned> element_index_set_1 = p_neighbour_1->rGetContainingElementIndices();
                        unsigned element_index_1 = *element_index_set_1.begin();
                        c_vector<double, DIM> r_neighbour_1 = p_neighbour_1->rGetLocation();

                        unsigned boundary_neighbour_2 = boundary_neighbours[2];
                        Node<DIM>* p_neighbour_2= p_mesh->GetNode(boundary_neighbour_2);
                        std::set<unsigned> element_index_set_2 = p_neighbour_2->rGetContainingElementIndices();
                        unsigned element_index_2 = *element_index_set_2.begin();
                        c_vector<double, DIM> r_neighbour_2 = p_neighbour_2->rGetLocation();

                        if(element_index_0 == element_index_1)
                        {
                            if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_0,r_neighbour_2)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 0 and 2
                                p_mesh->PerformNodeMerge(p_neighbour_0,p_neighbour_2);

                                // // Move node 0 to the mid-point
                                // p_neighbour_0->rGetModifiableLocation() = r_neighbour_0 + 0.5 * p_mesh->GetVectorFromAtoB(r_neighbour_0, r_neighbour_2);

                                // // Update the elements previously containing node 2 to contain node 0
                                // // Find the local index of node 2 in this element
                                // VertexElement<DIM,DIM>* p_element_2 = p_mesh->GetElement(element_index_2);
                                // unsigned node_2_local_index = p_element_2->GetNodeLocalIndex(boundary_neighbour_2);

                                // /*
                                // * If this element already contains node 0, then just remove node 2.
                                // * Otherwise replace it with node 0 in the element and remove it from mNodes.
                                // */
                                // // if (element_index_set_0.count(element_index_2) != 0)
                                // // {
                                //     TRACE("Deleted");
                                //     p_element_2->DeleteNode(p_element_2->GetNodeLocalIndex(node_2_local_index));
                                //     // p_element_2->DeleteNode(p_element_2->GetNodeLocalIndex(boundary_neighbour_2));
                                //     p_mesh->DeleteNodePriorToReMesh(boundary_neighbour_2);
                                // // }
                                // // else
                                // // {
                                //     TRACE("Replaced");
                                //     // Replace node 2 with node 0 in this element
                                //     p_mesh->GetElement(element_index_2)->AddNode(p_neighbour_0, node_2_local_index);
                                //     // p_mesh->GetElement(element_index_2)->AddNode(p_neighbour_0, boundary_neighbour_2);

                                // // }


                                performed_edge_modifier = true;
                                performed_edge_modifier_2 = true;
                                // TRACE("Custom 1.2.1");

                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_2)) < distanceBetweenVerteciesThreshold)
                            {
                                // Merge 1 and 2
                                p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);

                                performed_edge_modifier = true;
                                // TRACE("Custom 1.2.2");
                                performed_edge_modifier_2 = true;

                            }
                        }
                        else if(element_index_1 == element_index_2)
                        {
                            if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_0)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 0 and 1
                                p_mesh->PerformNodeMerge(p_neighbour_0,p_neighbour_1);

                                performed_edge_modifier = true;
                                // TRACE("Custom 1.2.3");
                                performed_edge_modifier_2 = true;

                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_2,r_neighbour_0)) < distanceBetweenVerteciesThreshold)
                            {
                                // Merge 0 and 2
                                p_mesh->PerformNodeMerge(p_neighbour_0,p_neighbour_2);

                                performed_edge_modifier = true;
                                // TRACE("Custom 1.2.4");
                                performed_edge_modifier_2 = true;
                                
                            }
                        }
                        else if(element_index_2 == element_index_0)
                        {
                            if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_2,r_neighbour_1)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 1 and 2
                                p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);

                                performed_edge_modifier = true;
                                // TRACE("Custom 1.2.5");
                                performed_edge_modifier_2 = true;

                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_0,r_neighbour_1)) < distanceBetweenVerteciesThreshold)
                            {
                                // Merge 0 and 1
                                p_mesh->PerformNodeMerge(p_neighbour_0,p_neighbour_1);

                                performed_edge_modifier = true;
                                // TRACE("Custom 1.2.6");
                                performed_edge_modifier_2 = true;
                                
                            }
                        }

                    }
                }
                
                else if(containing_element_indices.size() == 1)
                {
                    /*  |    |                       |   |
                    *   oD   oB                      oD  oB
                    *    \   |  Cell                   \ |  Cell
                    *     \  |                   Cell   \|
                    * Cell \ |                    _______oA
                    * _____Ao|           ---->           |
                    *        |                           |
                    *        |                           |
                    *        oC                          oC
                    */
                    // push a node into an edge pf a neighbouring cell...

                    for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter_2 = p_mesh->GetNodeIteratorBegin();
                        node_iter_2 != p_mesh->GetNodeIteratorEnd();
                        ++node_iter_2)
                    {
                        // Make sure node_2 is a boundary node
                        if (node_iter_2->IsBoundaryNode())
                        {
                             unsigned elem_index_1 = (*containing_element_indices.begin());
                            // make sure node_2 and node 1 are not in the same element
                            std::set<unsigned> containing_element_indices_2 = node_iter_2->rGetContainingElementIndices();
                            bool node_2_and_node_in_same_element = containing_element_indices_2.find(elem_index_1) != containing_element_indices_2.end();
                            if( !node_2_and_node_in_same_element )
                            {
                                c_vector<double, DIM> r_node_1 = node_iter->rGetLocation();
                                c_vector<double, DIM> r_node_2 = node_iter_2->rGetLocation();

                                // 1 and 2 are close, check if 1 is close to the edge of node 2 and its neighbours
                                if(norm_2(p_mesh->GetVectorFromAtoB(r_node_1,r_node_2)) < 0.05)
                                {
                                    unsigned node_2_index = node_iter_2->GetIndex();
                                    std::set<unsigned> node_2_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_2_index);

                                    std::vector<unsigned> boundary_neighbours_2;
                                    for (std::set<unsigned>::iterator neighbour_iter = node_2_neighbours.begin();
                                        neighbour_iter != node_2_neighbours.end();
                                        ++neighbour_iter)
                                    {
                                        unsigned neighbour_index = *neighbour_iter;
                                        Node<DIM>* p_neighbour_node = p_mesh->GetNode(neighbour_index);
                                        // std::set<unsigned> element_index_set_n = p_neighbour_node->rGetContainingElementIndices();

                                        if(p_neighbour_node->IsBoundaryNode() )
                                        {
                                            c_vector<double, DIM> r_neighbour = p_neighbour_node->rGetLocation();
                                                                                        
                                            double distance_to_edge = DistanceToEdgeFromPoint(r_neighbour,r_node_2,r_node_1);

                                            if(distance_to_edge < 0.05)
                                            {
                                                
                                                unsigned node_index = node_iter->GetIndex();

                                                Node<DIM>* p_neighbour_1= p_mesh->GetNode(node_index);
                                                Node<DIM>* p_neighbour_2= p_mesh->GetNode(node_2_index);
                                                p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);

                                                TRACE("Custom 1.3.1");

                                                p_neighbour_1->SetAsBoundaryNode(true);
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
            }
        }
    }

    /*
    *        |                           |
    *  Cell  |  Cell               Cell  |  Cell
    *        oB                          |
    *       /|                    _______oB
    * ___ Ao |           ---->           |
    *        |                           |
    *        |                           |
    *        oC                          oC
    *
    *   if vertex_A intersects edge defined by vertex_B-vertex_C
    *   then we delete vertex_A and move vertex_B to the midpoint
    *   of vertex_A-vertex_B
    */
    double distanceToEdgeThreshold = 0.05;
    if(true)
    {
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
                        // std::set<unsigned> element_index_set_n = p_neighbour_node->rGetContainingElementIndices();

                        // if(p_neighbour_node->IsBoundaryNode() && (element_index_set_n.size()==1) )
                        if(p_neighbour_node->IsBoundaryNode() )
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

                        if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1, r_neighbour_2)) < 0.25 )
                        {
                            if(norm_2(p_mesh->GetVectorFromAtoB(r_node, r_neighbour_2)) < norm_2(p_mesh->GetVectorFromAtoB(r_node, r_neighbour_1)) )
                            {
                                // p_neighbour_2 is closer to p_node, so check if p_neighbour_2 intersects the edge p_node-p_neighbour_1

                                // double numerator_signed = (r_neighbour_1[0] - r_node[0])*(r_node[1] - r_neighbour_2[1]) - (r_node[0] - r_neighbour_2[0])*(r_neighbour_1[1] - r_node[1]);
                                // double numerator_signed = (r_node[0]-r_neighbour_1[0])*(r_neighbour_1[1]-r_neighbour_2[1]) - (r_neighbour_1[0]-r_neighbour_2[0])*(r_node[1]-r_neighbour_1[1]);
                                // double numerator = sqrt(numerator_signed*numerator_signed);
                                // double denominator = sqrt((r_neighbour_1[0] - r_node[0])*(r_neighbour_1[0] - r_node[0]) + (r_neighbour_1[1] - r_node[1])*(r_neighbour_1[1] - r_node[1]));
                                
                                // double distance_p_2_to_edge = numerator/denominator;
                                double distance_p_2_to_edge = DistanceToEdgeFromPoint(r_node,r_neighbour_1,r_neighbour_2);

                                // TRACE(" ");
                                // PRINT_VECTOR(r_neighbour_2);
                                // PRINT_VECTOR(r_node);
                                // PRINT_VECTOR(r_neighbour_1);
                                // PRINT_VARIABLE(distance_p_2_to_edge);
                                if(distance_p_2_to_edge < distanceToEdgeThreshold)
                                {
                                    std::set<unsigned> element_index_set_2 = p_neighbour_2->rGetContainingElementIndices();
                                    if(element_index_set_2.size() == 1)
                                    {
                                        unsigned elem_index_2 = (*element_index_set_2.begin());
                                        VertexElement<DIM,DIM>* p_element_2 = p_mesh->GetElement(elem_index_2);
                                        p_element_2->DeleteNode(p_element_2->GetNodeLocalIndex(boundary_neighbours[1]));
                                        p_mesh->DeleteNodePriorToReMesh(boundary_neighbours[1]);

                                        p_node->rGetModifiableLocation() = r_neighbour_2;
                                        performed_edge_modifier = true;
                                        // TRACE("Custom 2.1");
                                        // PRINT_VECTOR(p_node->rGetModifiableLocation());
                                    }

                                }
                            }
                            
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_node, r_neighbour_1)) < norm_2(p_mesh->GetVectorFromAtoB(r_node, r_neighbour_2)))
                            {
                                // p_neighbour_1 is closer to p_node, so check if it intersects the edge p_node-p_neighbour_2

                                // double numerator_signed = (r_neighbour_2[0] - r_node[0])*(r_node[1] - r_neighbour_1[1]) - (r_node[0] - r_neighbour_1[0])*(r_neighbour_2[1] - r_node[1]);
                                // double numerator_signed = (r_node[0]-r_neighbour_2[0])*(r_neighbour_2[1]-r_neighbour_1[1]) - (r_neighbour_2[0]-r_neighbour_1[0])*(r_node[1]-r_neighbour_2[1]);
                                // double numerator = sqrt(numerator_signed*numerator_signed);
                                // double denominator = sqrt((r_neighbour_2[0] - r_node[0])*(r_neighbour_2[0] - r_node[0]) + (r_neighbour_2[1] - r_node[1])*(r_neighbour_2[1] - r_node[1]));
                                
                                // double distance_p_1_to_edge = numerator/denominator;
                                double distance_p_1_to_edge = DistanceToEdgeFromPoint(r_node,r_neighbour_2,r_neighbour_1);

                                // TRACE(" ");
                                // PRINT_VECTOR(r_neighbour_1);
                                // PRINT_VECTOR(r_node);
                                // PRINT_VECTOR(r_neighbour_2);
                                // PRINT_VARIABLE(distance_p_1_to_edge);
                                if(distance_p_1_to_edge < distanceToEdgeThreshold)
                                {
                                    std::set<unsigned> element_index_set_1 = p_neighbour_1->rGetContainingElementIndices();
                                    if( element_index_set_1.size() == 1 )
                                    {
                                        unsigned elem_index_1 = (*element_index_set_1.begin());
                                        VertexElement<DIM,DIM>* p_element_1 = p_mesh->GetElement(elem_index_1);
                                        p_element_1->DeleteNode(p_element_1->GetNodeLocalIndex(boundary_neighbours[0]));
                                        p_mesh->DeleteNodePriorToReMesh(boundary_neighbours[0]);

                                        p_node->rGetModifiableLocation() = r_neighbour_1;
                                        performed_edge_modifier = true;
                                        // TRACE("Custom 2.2");
                                        // PRINT_VECTOR(p_node->rGetModifiableLocation());

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Might have this sometimes due to boundary node misslabelling... 
    /*  --------o---------
    *           |
    *  (Cell)   o   (Cell)
    *           |
    *           |
    *   --------o----------
    */
    // i.e. internal nodes...
    if(true)
    {
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
                    
                    if(node_neighbours.size()==2)
                    {
                        std::set<unsigned>::iterator neigh_it = node_neighbours.begin();
                        unsigned node_neighbour_1 = (*neigh_it);
                        std::advance(neigh_it, 1);
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
                            Node<DIM>* p_node = p_mesh->GetNode(node_index);

                            std::set<unsigned>::iterator elem_it = containing_element_indices.begin();


                            unsigned elem_index_1 = (*elem_it);
                            VertexElement<DIM,DIM>* p_element_1 = p_mesh->GetElement(elem_index_1);
                            p_element_1->DeleteNode(p_element_1->GetNodeLocalIndex(node_index));

                            std::advance(elem_it, 1);
                            unsigned elem_index_2 = (*elem_it);
                            VertexElement<DIM,DIM>* p_element_2 = p_mesh->GetElement(elem_index_2);
                            p_element_2->DeleteNode(p_element_2->GetNodeLocalIndex(node_index));

                            p_mesh->DeleteNodePriorToReMesh(node_index);

                            TRACE("Custom 3.1");
                            performed_edge_modifier = true;
                        }
                    }
                }
            }
        }
    }

    double mDistanceForT3SwapChecking = 5.0;
    double mCellRearrangementRatio = 1.5;
    double mCellRearrangementThreshold = 0.1;
    if(false)
    {
        // If checking for T3 swaps, check that no boundary nodes have overlapped any boundary elements
        // First: find all boundary element and calculate their centroid only once
        std::vector<unsigned> boundary_element_indices;
        std::vector< c_vector<double, DIM> > boundary_element_centroids;
        for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_mesh->GetElementIteratorBegin();
                elem_iter != p_mesh->GetElementIteratorEnd();
                ++elem_iter)
        {
            if (elem_iter->IsElementOnBoundary())
            {
                unsigned element_index = elem_iter->GetIndex();
                boundary_element_indices.push_back(element_index);
                // should be a map but I am too lazy to look up the syntax
                boundary_element_centroids.push_back(p_mesh->GetCentroidOfElement(element_index));
            }
        }

        for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
            node_iter != p_mesh->GetNodeIteratorEnd();
            ++node_iter)
        {
            if (node_iter->IsBoundaryNode())
            {
                // index in boundary_element_centroids and boundary_element_indices
                unsigned boundary_element_index = 0;
                for (std::vector<unsigned>::iterator elem_iter = boundary_element_indices.begin();
                        elem_iter != boundary_element_indices.end();
                        ++elem_iter)
                {
                    // Check that the node is not part of this element
                    if (node_iter->rGetContainingElementIndices().count(*elem_iter) == 0)
                    {
                        c_vector<double, DIM> node_location = node_iter->rGetLocation();
                        c_vector<double, DIM> element_centroid = boundary_element_centroids[boundary_element_index];
                        double node_element_distance = norm_2(p_mesh->GetVectorFromAtoB(node_location, element_centroid));
                        
                        if ( node_element_distance < mDistanceForT3SwapChecking )
                        {
                            // PRINT_VARIABLE(node_element_distance);

                            
                            bool element_includes_point = false;
                            if(true) // (p_mesh->ElementIncludesPoint(node_iter->rGetLocation(), *elem_iter))
                            {
                                // unsigned elementIndex = elem_iter->GetIndex();

                                VertexElement<DIM,DIM>* p_element = p_mesh->GetElement(*elem_iter);
                                unsigned num_nodes = p_element->GetNumNodes();

                                // Initialise boolean
                                
                                c_vector<double, DIM> first_vertex = p_element->GetNodeLocation(0);
                                c_vector<double, DIM> test_point = p_mesh->GetVectorFromAtoB(first_vertex, node_iter->rGetLocation());

                                // Loop over edges of the element
                                c_vector<double, DIM> vertexA = zero_vector<double>(DIM);
                                for (unsigned local_index = 0; local_index < num_nodes; local_index++)
                                {
                                    // Check if this edge crosses the ray running out horizontally (increasing x, fixed y) from the test point
                                    c_vector<double, DIM> vector_a_to_point = p_mesh->GetVectorFromAtoB(vertexA, test_point);

                                    // Pathological case - test point coincides with vertexA
                                    // (we will check vertexB next time we go through the for loop)
                                    if (norm_2(vector_a_to_point) < DBL_EPSILON)
                                    {
                                        element_includes_point = false;
                                    }

                                    c_vector<double, DIM> vertexB = p_mesh->GetVectorFromAtoB(first_vertex, p_element->GetNodeLocation((local_index + 1) % num_nodes));
                                    c_vector<double, DIM> vector_b_to_point = p_mesh->GetVectorFromAtoB(vertexB, test_point);
                                    c_vector<double, DIM> vector_a_to_b = p_mesh->GetVectorFromAtoB(vertexA, vertexB);

                                    // Pathological case - ray coincides with horizontal edge
                                    if ((fabs(vector_a_to_b[1]) < DBL_EPSILON) && (fabs(vector_a_to_point[1]) < DBL_EPSILON) && (fabs(vector_b_to_point[1]) < DBL_EPSILON))
                                    {
                                        if ((vector_a_to_point[0] > 0) != (vector_b_to_point[0] > 0))
                                        {
                                            element_includes_point = false;
                                        }
                                    }

                                    // Non-pathological case
                                    // A and B on different sides of the line y = test_point[1]
                                    if ((vertexA[1] > test_point[1]) != (vertexB[1] > test_point[1]))
                                    {
                                        // Intersection of y=test_point[1] and vector_a_to_b is on the right of test_point
                                        if (test_point[0] < vertexA[0] + vector_a_to_b[0] * vector_a_to_point[1] / vector_a_to_b[1])
                                        {
                                            element_includes_point = !element_includes_point;
                                        }
                                    }

                                    vertexA = vertexB;
                                }
                                // element_includes_point;
                            }


                            bool is_point_close_to_edge = false;
                            // // Check to see if the node is close to the edge
                            // if(true)
                            // {

                            // }



                            // if (p_mesh->ElementIncludesPoint(node_iter->rGetLocation(), *elem_iter))
                            if(element_includes_point || is_point_close_to_edge)
                            {
                                TRACE("MY T3 1.2");
                                // this->PerformT3Swap(&(*node_iter), *elem_iter);
                                
                                unsigned node_index = node_iter->GetIndex();
                                Node<DIM>* pNode = p_mesh->GetNode(node_index);
                                // unsigned elementIndex = elem_iter->GetIndex();
                                
                                
                                // Get element
                                VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(*elem_iter);
                                unsigned num_nodes = p_element->GetNumNodes();
                                
                                // Store the index of the elements containing the intersecting node
                                std::set<unsigned> elements_containing_intersecting_node = pNode->rGetContainingElementIndices();

                                // Get the local index of the node in the intersected element after which the new node is to be added
                                // unsigned node_A_local_index = p_mesh->GetLocalIndexForElementEdgeClosestToPoint(pNode->rGetLocation(), *elem_iter);
                                unsigned node_A_local_index;
                                if(true)
                                {
                                    // Get the element
                                    VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(*elem_iter);
                                    unsigned num_nodes = p_element->GetNumNodes();

                                    double min_squared_normal_distance = DBL_MAX;
                                    unsigned min_distance_edge_index = UINT_MAX;

                                    // Loop over edges of the element
                                    for (unsigned local_index = 0; local_index < num_nodes; local_index++)
                                    {
                                        // Get the end points of this edge
                                        c_vector<double, DIM> vertexA = p_element->GetNodeLocation(local_index);
                                        c_vector<double, DIM> vertexB = p_element->GetNodeLocation((local_index + 1) % num_nodes);

                                        c_vector<double, DIM> vector_a_to_point = p_mesh->GetVectorFromAtoB(vertexA, pNode->rGetLocation());
                                        c_vector<double, DIM> vector_a_to_b = p_mesh->GetVectorFromAtoB(vertexA, vertexB);
                                        double distance_a_to_b = norm_2(vector_a_to_b);

                                        c_vector<double, DIM> edge_ab_unit_vector = vector_a_to_b / norm_2(vector_a_to_b);
                                        double distance_parallel_to_edge = inner_prod(vector_a_to_point, edge_ab_unit_vector);

                                        double squared_distance_normal_to_edge = SmallPow(norm_2(vector_a_to_point), 2) - SmallPow(distance_parallel_to_edge, 2);

                                        /*
                                        * If the point lies almost bang on the supporting line of the edge, then snap to the line.
                                        * This allows us to do floating point tie-breaks when line is exactly at a node.
                                        * We adopt a similar approach if the point is at the same position as a point in the
                                        * element.
                                        */
                                        if (squared_distance_normal_to_edge < DBL_EPSILON)
                                        {
                                            squared_distance_normal_to_edge = 0.0;
                                        }

                                        if (fabs(distance_parallel_to_edge) < DBL_EPSILON)
                                        {
                                            distance_parallel_to_edge = 0.0;
                                        }
                                        else if (fabs(distance_parallel_to_edge - distance_a_to_b) < DBL_EPSILON)
                                        {
                                            distance_parallel_to_edge = distance_a_to_b;
                                        }

                                        // Make sure node is within the confines of the edge and is the nearest edge to the node \this breaks for convex elements
                                        if (squared_distance_normal_to_edge < min_squared_normal_distance && distance_parallel_to_edge >= 0 && distance_parallel_to_edge <= distance_a_to_b)
                                        {
                                            min_squared_normal_distance = squared_distance_normal_to_edge;
                                            min_distance_edge_index = local_index;
                                        }
                                    }

                                    node_A_local_index =  min_distance_edge_index;
                                }

                                // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
                                c_vector<double, DIM> node_location;
                                node_location = pNode->rGetModifiableLocation();

                                // Get the nodes at either end of the edge to be divided
                                unsigned vertexA_index = p_element->GetNodeGlobalIndex(node_A_local_index);
                                unsigned vertexB_index = p_element->GetNodeGlobalIndex((node_A_local_index+1)%num_nodes);

                                // Check these nodes are also boundary nodes if this fails then the elements have become concave and you need a smaller timestep
                                // if (!p_mesh->mNodes[vertexA_index]->IsBoundaryNode() || !p_mesh->mNodes[vertexB_index]->IsBoundaryNode())
                                // {
                                //     EXCEPTION("A boundary node has intersected a non-boundary edge; this is because the boundary element has become concave. You need to rerun the simulation with a smaller time step to prevent this.");
                                // }

                                // Get the nodes at either end of the edge to be divided and calculate intersection
                                c_vector<double, DIM> vertexA = p_element->GetNodeLocation(node_A_local_index);
                                c_vector<double, DIM> vertexB = p_element->GetNodeLocation((node_A_local_index+1)%num_nodes);
                                c_vector<double, DIM> vector_a_to_point = p_mesh->GetVectorFromAtoB(vertexA, node_location);

                                c_vector<double, DIM> vector_a_to_b = p_mesh->GetVectorFromAtoB(vertexA, vertexB);

                                c_vector<double, DIM> edge_ab_unit_vector = vector_a_to_b/norm_2(vector_a_to_b);
                                c_vector<double, DIM> intersection = vertexA + edge_ab_unit_vector*inner_prod(vector_a_to_point, edge_ab_unit_vector);

                                

                                /*
                                *  From          To
                                *   ____        _______
                                *                 / \
                                *    /\   ^      /   \
                                *   /  \  |
                                *
                                *  The edge goes from vertexA--vertexB to vertexA--new_node--pNode--vertexB
                                */

                                // Check whether the intersection location fits into the edge and update distances and vertex positions afterwards.
                                // intersection = p_mesh->WidenEdgeOrCorrectIntersectionLocationIfNecessary(vertexA_index, vertexB_index, intersection);
                                edge_ab_unit_vector = p_mesh->GetPreviousEdgeGradientOfElementAtNode(p_element, (node_A_local_index+1)%num_nodes);

                                // Move original node
                                pNode->rGetModifiableLocation() = intersection + 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;

                                // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
                                c_vector<double, DIM> new_node_location;
                                new_node_location = intersection - 0.5*mCellRearrangementRatio*mCellRearrangementThreshold*edge_ab_unit_vector;

                                // Add new node which will always be a boundary node
                                unsigned new_node_global_index = p_mesh->AddNode(new Node<DIM>(0, true, new_node_location[0], new_node_location[1]));
                                Node<DIM>* p_new_node = p_mesh->GetNode(new_node_global_index);

                                // Add the moved and new nodes to the element (this also updates the node)
                                p_mesh->GetElement(*elem_iter)->AddNode(pNode, node_A_local_index);
                                // p_mesh->GetElement(*elem_iter)->AddNode(p_mesh->mNodes[new_node_global_index], node_A_local_index);
                                p_mesh->GetElement(*elem_iter)->AddNode(p_new_node, node_A_local_index);

                                // Add the new node to the original element containing pNode (this also updates the node)
                                unsigned intersecting_element_index = *elements_containing_intersecting_node.begin();
                                // p_mesh->GetElement(intersecting_element_index)->AddNode(p_mesh->mNodes[new_node_global_index], p_mesh->GetElement(intersecting_element_index)->GetNodeLocalIndex(pNode->GetIndex()));
                                p_mesh->GetElement(intersecting_element_index)->AddNode(p_new_node, p_mesh->GetElement(intersecting_element_index)->GetNodeLocalIndex(pNode->GetIndex()));

                                // The nodes must have been updated correctly
                                // assert(pNode->GetNumContainingElements() == 2);
                                // assert(p_mesh->mNodes[new_node_global_index]->GetNumContainingElements() == 2);


                            }
                        }
                    }
                    // increment the boundary element index
                    boundary_element_index +=1u;
                }

            }
        }
    }



    if(performed_edge_modifier)
    {
        p_mesh->ReMesh();
        // PRINT_VARIABLE(SimulationTime::Instance()->GetTime());

    }
    // if(performed_edge_modifier_2)
    // {
    //     PRINT_VARIABLE(SimulationTime::Instance()->GetTime());
    // }
    
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

