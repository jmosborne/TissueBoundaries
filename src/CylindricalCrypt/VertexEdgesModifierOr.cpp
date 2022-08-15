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
bool VertexEdgesModifier<DIM>::IsNodeABNeighbours( unsigned node_a,unsigned node_b, AbstractCellPopulation<DIM,DIM>& rCellPopulation )
{
    bool is_neighbour = false;
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_a);

    for (std::set<unsigned>::iterator neighbour_iter = node_neighbours.begin();
        neighbour_iter != node_neighbours.end();
        ++neighbour_iter)
    {
        unsigned neighbour_index = *neighbour_iter;
        if(neighbour_index == node_b)
        {
            is_neighbour = true;
        }
    }
    return is_neighbour;

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

    PRINT_VARIABLE(SimulationTime::Instance()->GetTime());


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
                        *     \ (Cell) /                 \ (Cell)  /
                        *      \      /                   \       /
                        *       \   /                      \     /
                        *        \ /           ---->        \   /
                        * (Cell)  o  (Cell)           (Cell)  o  (Cell)  
                        *       /   \                         |
                        *      / (v) \                        |
                        * ----o       o -----           ------o-------
                        *
                        */
                        if( containing_element_indices.size() == 3 )
                        {
                            if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_2)) < distanceBetweenVerteciesThreshold)  
                            {
                                p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);
                                p_node->SetAsBoundaryNode(false);

                                performed_edge_modifier = true;
                                TRACE("Performed 3 cell edge stitch");
                            }

                        }
                        else
                        {
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
                                p_node->SetAsBoundaryNode(true);

                                performed_edge_modifier = true;
                                TRACE("Performed 3 cell edge stitch");

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
                        // Node<DIM>* p_node = p_mesh->GetNode(node_index);


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
                                p_neighbour_0->SetAsBoundaryNode(true);

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
                                TRACE("Performed 2 cell edge stitch with void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_2)) < distanceBetweenVerteciesThreshold)
                            {
                                // Merge 1 and 2
                                p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);
                                p_neighbour_1->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with void");
                            }
                        }
                        else if(element_index_1 == element_index_2)
                        {
                            if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_0)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 0 and 1
                                p_mesh->PerformNodeMerge(p_neighbour_0,p_neighbour_1);
                                p_neighbour_0->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_2,r_neighbour_0)) < distanceBetweenVerteciesThreshold)
                            {
                                // Merge 0 and 2
                                p_mesh->PerformNodeMerge(p_neighbour_0,p_neighbour_2);
                                p_neighbour_0->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with void");
                            }
                        }
                        else if(element_index_2 == element_index_0)
                        {
                            if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_2,r_neighbour_1)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 1 and 2
                                p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);
                                p_neighbour_1->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_0,r_neighbour_1)) < distanceBetweenVerteciesThreshold)
                            {
                                // Merge 0 and 1
                                p_mesh->PerformNodeMerge(p_neighbour_0,p_neighbour_1);
                                p_neighbour_0->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with void");
                            }
                        }
                        // TRACE("line 342");
                    }

                    /* We might have the weird case as bellow:
                    *              \   (v)  /                              \   (v)  /      
                    *               \      /                                \      /
                    *             (3)o    o(0)                               o    o       
                    *                 \  /                                    \  /
                    *                  \/                                      \/
                    *      (Cell)      o(p)  (Cell)   ---------->      (Cell)  o   (Cell)
                    *                 / |                                      |
                    *                /  |                                      |
                    *               /(v)|                                      |
                    *        ------o    o------                        --------o-------
                    *             (1)  (2)                                  (v)
                    * 
                    *   Here, nodes 0,1,2,3 can all be boundary nodes which 
                    *   only belong to a single cell(element) each.
                    */
                    else if(boundary_neighbours.size() == 4)
                    {
                        // Node<DIM>* p_node = p_mesh->GetNode(node_index);


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

                        unsigned boundary_neighbour_3 = boundary_neighbours[3];
                        Node<DIM>* p_neighbour_3= p_mesh->GetNode(boundary_neighbour_3);
                        std::set<unsigned> element_index_set_3 = p_neighbour_3->rGetContainingElementIndices();
                        unsigned element_index_3 = *element_index_set_3.begin();
                        c_vector<double, DIM> r_neighbour_3 = p_neighbour_3->rGetLocation();


                        if(element_index_0 == element_index_1)
                        {
                            if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_0,r_neighbour_2)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 0 and 2
                                p_mesh->PerformNodeMerge(p_neighbour_0,p_neighbour_2);
                                p_neighbour_0->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_0,r_neighbour_3)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 0 and 3
                                p_mesh->PerformNodeMerge(p_neighbour_0,p_neighbour_3);
                                p_neighbour_0->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_2,r_neighbour_1)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 2 and 1
                                p_mesh->PerformNodeMerge(p_neighbour_2,p_neighbour_1);
                                p_neighbour_2->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_3)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 1 and 3
                                p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_3);
                                p_neighbour_1->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                        }
                        else if(element_index_0 == element_index_2)
                        {
                            if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_0,r_neighbour_1)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 0 and 1
                                p_mesh->PerformNodeMerge(p_neighbour_0,p_neighbour_1);
                                p_neighbour_0->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_0,r_neighbour_3)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 0 and 3
                                p_mesh->PerformNodeMerge(p_neighbour_0,p_neighbour_3);
                                p_neighbour_0->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_2,r_neighbour_1)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 2 and 1
                                p_mesh->PerformNodeMerge(p_neighbour_2,p_neighbour_1);
                                p_neighbour_2->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_2,r_neighbour_3)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 2 and 3
                                p_mesh->PerformNodeMerge(p_neighbour_2,p_neighbour_3);
                                p_neighbour_2->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                        }
                        else if(element_index_0 == element_index_3)
                        {
                            if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_0,r_neighbour_1)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 0 and 1
                                p_mesh->PerformNodeMerge(p_neighbour_0,p_neighbour_1);
                                p_neighbour_0->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_0,r_neighbour_2)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 0 and 2
                                p_mesh->PerformNodeMerge(p_neighbour_0,p_neighbour_2);
                                p_neighbour_0->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_3,r_neighbour_1)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 3 and 1
                                p_mesh->PerformNodeMerge(p_neighbour_3,p_neighbour_1);
                                p_neighbour_3->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_2,r_neighbour_3)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 2 and 3
                                p_mesh->PerformNodeMerge(p_neighbour_2,p_neighbour_3);
                                p_neighbour_2->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                        }
                        else if(element_index_1 == element_index_2)
                        {
                            if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_0,r_neighbour_1)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 0 and 1
                                p_mesh->PerformNodeMerge(p_neighbour_0,p_neighbour_1);
                                p_neighbour_0->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_3)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 1 and 3
                                p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_3);
                                p_neighbour_1->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_2,r_neighbour_0)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 2 and 0
                                p_mesh->PerformNodeMerge(p_neighbour_2,p_neighbour_0);
                                p_neighbour_2->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_2,r_neighbour_3)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 2 and 3
                                p_mesh->PerformNodeMerge(p_neighbour_2,p_neighbour_3);
                                p_neighbour_2->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                        }
                        else if(element_index_1 == element_index_3)
                        {
                            if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_0,r_neighbour_1)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 0 and 1
                                p_mesh->PerformNodeMerge(p_neighbour_0,p_neighbour_1);
                                p_neighbour_0->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_2)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 1 and 2
                                p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);
                                p_neighbour_1->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_3,r_neighbour_0)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 3 and 0
                                p_mesh->PerformNodeMerge(p_neighbour_3,p_neighbour_0);
                                p_neighbour_3->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_2,r_neighbour_3)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 2 and 3
                                p_mesh->PerformNodeMerge(p_neighbour_2,p_neighbour_3);
                                p_neighbour_2->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                        }
                        else if(element_index_2 == element_index_3)
                        {
                            if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_0,r_neighbour_2)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 0 and 2
                                p_mesh->PerformNodeMerge(p_neighbour_0,p_neighbour_2);
                                p_neighbour_3->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_2)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 1 and 2
                                p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);
                                p_neighbour_1->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_3,r_neighbour_0)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 3 and 0
                                p_mesh->PerformNodeMerge(p_neighbour_3,p_neighbour_0);
                                p_neighbour_3->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_3)) < distanceBetweenVerteciesThreshold)  
                            {
                                // Merge 1 and 3
                                p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_3);
                                p_neighbour_1->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                TRACE("Performed 2 cell edge stitch with large void");
                            }
                        }
                        // TRACE("line 553");
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

                                                p_neighbour_1->SetAsBoundaryNode(true);
                                                performed_edge_modifier = true;
                                                TRACE("Performed node merge into edge");
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

    if(performed_edge_modifier)
    {
        p_mesh->ReMesh();
        TRACE("1st_mesh");

        unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
        std::stringstream time;
        time << num_timesteps;
        PRINT_VARIABLE(num_timesteps);
        VertexMeshWriter<DIM,DIM> vertexmesh_writer("tmp", "1st_mesh_", false);
        vertexmesh_writer.WriteVtkUsingMesh(*p_mesh, time.str());

        // PRINT_VARIABLE(SimulationTime::Instance()->GetTime());
    }
    performed_edge_modifier = false;

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

                                double distance_p_2_to_edge = DistanceToEdgeFromPoint(r_node,r_neighbour_1,r_neighbour_2);

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
                                        p_node->SetAsBoundaryNode(true);
                                        performed_edge_modifier = true;
                                        TRACE("Performed node merge down edge");
                                    }

                                }
                            }
                            
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_node, r_neighbour_1)) < norm_2(p_mesh->GetVectorFromAtoB(r_node, r_neighbour_2)))
                            {

                                double distance_p_1_to_edge = DistanceToEdgeFromPoint(r_node,r_neighbour_2,r_neighbour_1);

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
                                        p_node->SetAsBoundaryNode(true);
                                        performed_edge_modifier = true;
                                        TRACE("Performed node merge down edge");
                                    }
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
        p_mesh->ReMesh();
        TRACE("2nd_mesh");
        unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
        std::stringstream time;
        time << num_timesteps;
        PRINT_VARIABLE(num_timesteps);
        VertexMeshWriter<DIM,DIM> vertexmesh_writer("tmp", "2nd_mesh", false);
        vertexmesh_writer.WriteVtkUsingMesh(*p_mesh, time.str());
    }
    performed_edge_modifier = false;

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
                            Node<DIM>* p_node = p_mesh->GetNode(node_index);

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
                            TRACE("remove internal node");
                        }
                    }
                }
            }
        }
    }
    if(performed_edge_modifier)
    {
        p_mesh->ReMesh();
        TRACE("3rd_mesh");
        unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
        std::stringstream time;
        time << num_timesteps;
        PRINT_VARIABLE(num_timesteps);
        VertexMeshWriter<DIM,DIM> vertexmesh_writer("tmp", "3rd_mesh", false);
        vertexmesh_writer.WriteVtkUsingMesh(*p_mesh, time.str());
    }
    performed_edge_modifier = false;

    // Void removal
    double mVoidAreaThreshold = 0.01;
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
                    
                    if(node_neighbours.size() == 3)
                    {
                        std::set<unsigned>::const_iterator neigh_it = node_neighbours.begin();
                        unsigned node_neighbour_1 = (*neigh_it);
                        neigh_it++;
                        unsigned node_neighbour_2 = (*neigh_it);
                        neigh_it++;
                        unsigned node_neighbour_3 = (*neigh_it);

                        Node<DIM>* p_neighbour_node_1 = p_mesh->GetNode(node_neighbour_1);
                        Node<DIM>* p_neighbour_node_2 = p_mesh->GetNode(node_neighbour_2);
                        Node<DIM>* p_neighbour_node_3 = p_mesh->GetNode(node_neighbour_3);

                        // Node<DIM>* p_boundary_1;
                        // Node<DIM>* p_boundary_2;



                        if(p_neighbour_node_1->IsBoundaryNode() && p_neighbour_node_2->IsBoundaryNode() && IsNodeABNeighbours(node_neighbour_1,node_neighbour_2, rCellPopulation ))
                        {

                            std::set<unsigned> containing_elements_1 = p_neighbour_node_1->rGetContainingElementIndices();
                            std::set<unsigned> containing_elements_2 = p_neighbour_node_2->rGetContainingElementIndices();

                            std::set<unsigned> shared_elements;
                            std::set_intersection(containing_elements_1.begin(),
                                            containing_elements_1.end(),
                                            containing_elements_2.begin(),
                                            containing_elements_2.end(),
                                            std::inserter(shared_elements, shared_elements.begin()));
                            if(containing_element_indices.size()==2 && containing_elements_2.size()==2 && shared_elements.size()==1)
                            {
                                Node<DIM>* p_boundary_1 = p_neighbour_node_1;
                                Node<DIM>* p_boundary_2 = p_neighbour_node_2;
                                
                                Node<DIM>* p_node = p_mesh->GetNode(node_index);
                                c_vector<double, DIM> r_node = p_node->rGetLocation();
                                c_vector<double, DIM> r_boundary_1 = p_boundary_1->rGetLocation();
                                c_vector<double, DIM> r_boundary_2 = p_boundary_2->rGetLocation();

                                double a = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_1));
                                double b = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_2));
                                double c = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_1, r_boundary_2));
                                double s = 0.5*(a+b+c);

                                double void_area = sqrt(s*(s-a)*(s-b)*(s-c));

                                if(void_area < mVoidAreaThreshold)
                                {
                                    // p_mesh->PerformVoidRemoval(p_node,p_boundary_1,p_boundary_2);
                                    c_vector<double, DIM> nodes_midpoint = p_node->rGetLocation()
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_1->rGetLocation()) / 3.0
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_2->rGetLocation()) / 3.0;


                                    PRINT_VECTOR(nodes_midpoint);

                                    /*          \             /
                                    *            \   (Cell)  /
                                    *             o---------o
                                    *              \       /
                                    *               \ (v) /
                                    *                \   /
                                    *     (Cell)       o     (Cell)
                                    *                /   \
                                    *               /     \
                                    *         -----o  (v)  o-----
                                    */
                                    // Decide whether the last merged node will be a boundary node or not...  CBFed right now...
                                    bool is_boundary = false;
                                    // Check the neighbours of A that are not B or C and check whether any of those are boundary nodes.
                                    // If they are a boundary node, then the merged void will be a boundary node. Otherwise, not!
                                    
                                    unsigned node_A_index = p_node->GetIndex();
                                    unsigned node_B_index = p_boundary_1->GetIndex();
                                    unsigned node_C_index = p_boundary_2->GetIndex();

                                    // Check As neighbours
                                    std::set<unsigned> first_neighbour_indices = p_node->rGetContainingElementIndices();
                                    for (std::set<unsigned>::iterator elem_iter = first_neighbour_indices.begin();
                                            elem_iter != first_neighbour_indices.end();
                                            ++elem_iter)
                                    {
                                        VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(*elem_iter);

                                        for (unsigned local_node_index = 0;
                                                local_node_index < p_element->GetNumNodes();
                                                ++local_node_index)
                                        {
                                            unsigned num_nodes = p_element->GetNumNodes();

                                            unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
                                            unsigned global_node_neigh = p_element->GetNodeGlobalIndex((local_node_index+1)%num_nodes);
                                            if(global_node_index != node_B_index && global_node_index != node_C_index )
                                            {
                                                if(global_node_neigh == node_A_index)
                                                {
                                                    Node<DIM>* p_node = p_mesh->GetNode(global_node_index);
                                                    std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
                                                    if(containing_elements.size() == 1)
                                                    {
                                                        if(p_node->IsBoundaryNode())
                                                        {   
                                                            is_boundary = true;
                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                            
                                        }
                                    }
                                    // Check Bs neighbours
                                    if(!is_boundary)
                                    {
                                        std::set<unsigned> first_neighbour_indices = p_boundary_1->rGetContainingElementIndices();
                                        for (std::set<unsigned>::iterator elem_iter = first_neighbour_indices.begin();
                                                elem_iter != first_neighbour_indices.end();
                                                ++elem_iter)
                                        {
                                            VertexElement<DIM, DIM>*  p_element = p_mesh->GetElement(*elem_iter);
                                            unsigned num_nodes = p_element->GetNumNodes();

                                            for (unsigned local_node_index = 0;
                                                    local_node_index < p_element->GetNumNodes();
                                                    ++local_node_index)
                                            {
                                                // first_neighbour_node_indices.insert(p_element->GetNodeGlobalIndex(local_node_index));
                                                unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
                                                unsigned global_node_neigh = p_element->GetNodeGlobalIndex((local_node_index+1)%num_nodes);
                                                if(global_node_index != node_A_index && global_node_index != node_C_index )
                                                {
                                                    if(global_node_neigh == node_B_index)
                                                    {
                                                        Node<DIM>* p_node = p_mesh->GetNode(global_node_index);
                                                        std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
                                                        if(containing_elements.size() == 1)
                                                        {
                                                            if(p_node->IsBoundaryNode())
                                                            {
                                                                is_boundary = true;
                                                                break;
                                                            }
                                                        }
                                                    }
                                                }
                                                
                                            }
                                        }
                                    }
                                    // Check Cs neighbours
                                    if(!is_boundary)
                                    {
                                        std::set<unsigned> first_neighbour_indices = p_boundary_2->rGetContainingElementIndices();
                                        for (std::set<unsigned>::iterator elem_iter = first_neighbour_indices.begin();
                                                elem_iter != first_neighbour_indices.end();
                                                ++elem_iter)
                                        {
                                            VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(*elem_iter);
                                            unsigned num_nodes = p_element->GetNumNodes();

                                            for (unsigned local_node_index = 0;
                                                    local_node_index < p_element->GetNumNodes();
                                                    ++local_node_index)
                                            {
                                                // first_neighbour_node_indices.insert(p_element->GetNodeGlobalIndex(local_node_index));
                                                unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
                                                unsigned global_node_neigh = p_element->GetNodeGlobalIndex((local_node_index+1)%num_nodes);

                                                if(global_node_index != node_B_index && global_node_index != node_A_index )
                                                {
                                                    if(global_node_neigh == node_C_index)
                                                    {
                                                        Node<DIM>* p_node = p_mesh->GetNode(global_node_index);
                                                        std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
                                                        if(containing_elements.size() == 1)
                                                        {
                                                            if(p_node->IsBoundaryNode())
                                                            {
                                                                is_boundary = true;
                                                                break;
                                                            }
                                                        }
                                                    }
                                                }
                                                
                                            }
                                        }
                                    }

                                    p_mesh->PerformNodeMerge(p_node, p_boundary_1);

                                    Node<DIM>* p_merged_node = p_boundary_1;

                                    if (p_boundary_1->IsDeleted())
                                    {
                                        p_merged_node = p_node;
                                    }

                                    p_mesh->PerformNodeMerge(p_boundary_2, p_merged_node);

                                    if (p_merged_node->IsDeleted())
                                    {
                                        p_merged_node = p_boundary_2;
                                    }

                                    p_merged_node->rGetModifiableLocation() = nodes_midpoint;

                                    // Tag remaining node as non-boundary, or boundary
                                    p_merged_node->SetAsBoundaryNode(is_boundary);

                                    PRINT_VARIABLE(is_boundary);

                                    performed_edge_modifier = true;
                                    TRACE("Performed Void Removel");
                                }
                            }
                        }
                        else if(p_neighbour_node_3->IsBoundaryNode() && p_neighbour_node_2->IsBoundaryNode() && IsNodeABNeighbours(node_neighbour_3,node_neighbour_2, rCellPopulation ))
                        {
                            std::set<unsigned> containing_elements_1 = p_neighbour_node_3->rGetContainingElementIndices();
                            std::set<unsigned> containing_elements_2 = p_neighbour_node_2->rGetContainingElementIndices();

                            std::set<unsigned> shared_elements;
                            std::set_intersection(containing_elements_1.begin(),
                                            containing_elements_1.end(),
                                            containing_elements_2.begin(),
                                            containing_elements_2.end(),
                                            std::inserter(shared_elements, shared_elements.begin()));
                            if(containing_element_indices.size()==2 && containing_elements_2.size()==2 && shared_elements.size()==1)
                            {
                                Node<DIM>* p_boundary_1 = p_neighbour_node_3;
                                Node<DIM>* p_boundary_2 = p_neighbour_node_2;
                                
                                Node<DIM>* p_node = p_mesh->GetNode(node_index);
                                c_vector<double, DIM> r_node = p_node->rGetLocation();
                                c_vector<double, DIM> r_boundary_1 = p_boundary_1->rGetLocation();
                                c_vector<double, DIM> r_boundary_2 = p_boundary_2->rGetLocation();

                                double a = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_1));
                                double b = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_2));
                                double c = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_1, r_boundary_2));
                                double s = 0.5*(a+b+c);

                                double void_area = sqrt(s*(s-a)*(s-b)*(s-c));

                                if(void_area < mVoidAreaThreshold)
                                {
                                    // p_mesh->PerformVoidRemoval(p_node,p_boundary_1,p_boundary_2);
                                    c_vector<double, DIM> nodes_midpoint = p_node->rGetLocation()
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_1->rGetLocation()) / 3.0
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_2->rGetLocation()) / 3.0;


                                    PRINT_VECTOR(nodes_midpoint);
                                    /*          \             /
                                    *            \   (Cell)  /
                                    *             o---------o
                                    *              \       /
                                    *               \ (v) /
                                    *                \   /
                                    *     (Cell)       o     (Cell)
                                    *                /   \
                                    *               /     \
                                    *         -----o  (v)  o-----
                                    */
                                    // Decide whether the last merged node will be a boundary node or not...  CBFed right now...
                                    bool is_boundary = false;
                                    // Check the neighbours of A that are not B or C and check whether any of those are boundary nodes.
                                    // If they are a boundary node, then the merged void will be a boundary node. Otherwise, not!
                                    
                                    unsigned node_A_index = p_node->GetIndex();
                                    unsigned node_B_index = p_boundary_1->GetIndex();
                                    unsigned node_C_index = p_boundary_2->GetIndex();

                                    // Check As neighbours
                                    std::set<unsigned> first_neighbour_indices = p_node->rGetContainingElementIndices();
                                    for (std::set<unsigned>::iterator elem_iter = first_neighbour_indices.begin();
                                            elem_iter != first_neighbour_indices.end();
                                            ++elem_iter)
                                    {
                                        VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(*elem_iter);

                                        for (unsigned local_node_index = 0;
                                                local_node_index < p_element->GetNumNodes();
                                                ++local_node_index)
                                        {
                                            unsigned num_nodes = p_element->GetNumNodes();

                                            unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
                                            unsigned global_node_neigh = p_element->GetNodeGlobalIndex((local_node_index+1)%num_nodes);
                                            if(global_node_index != node_B_index && global_node_index != node_C_index )
                                            {
                                                if(global_node_neigh == node_A_index)
                                                {
                                                    Node<DIM>* p_node = p_mesh->GetNode(global_node_index);
                                                    std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
                                                    if(containing_elements.size() == 1)
                                                    {
                                                        if(p_node->IsBoundaryNode())
                                                        {   
                                                            is_boundary = true;
                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                            
                                        }
                                    }
                                    // Check Bs neighbours
                                    if(!is_boundary)
                                    {
                                        std::set<unsigned> first_neighbour_indices = p_boundary_1->rGetContainingElementIndices();
                                        for (std::set<unsigned>::iterator elem_iter = first_neighbour_indices.begin();
                                                elem_iter != first_neighbour_indices.end();
                                                ++elem_iter)
                                        {
                                            VertexElement<DIM, DIM>*  p_element = p_mesh->GetElement(*elem_iter);
                                            unsigned num_nodes = p_element->GetNumNodes();

                                            for (unsigned local_node_index = 0;
                                                    local_node_index < p_element->GetNumNodes();
                                                    ++local_node_index)
                                            {
                                                // first_neighbour_node_indices.insert(p_element->GetNodeGlobalIndex(local_node_index));
                                                unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
                                                unsigned global_node_neigh = p_element->GetNodeGlobalIndex((local_node_index+1)%num_nodes);
                                                if(global_node_index != node_A_index && global_node_index != node_C_index )
                                                {
                                                    if(global_node_neigh == node_B_index)
                                                    {
                                                        Node<DIM>* p_node = p_mesh->GetNode(global_node_index);
                                                        std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
                                                        if(containing_elements.size() == 1)
                                                        {
                                                            if(p_node->IsBoundaryNode())
                                                            {
                                                                is_boundary = true;
                                                                break;
                                                            }
                                                        }
                                                    }
                                                }
                                                
                                            }
                                        }
                                    }
                                    // Check Cs neighbours
                                    if(!is_boundary)
                                    {
                                        std::set<unsigned> first_neighbour_indices = p_boundary_2->rGetContainingElementIndices();
                                        for (std::set<unsigned>::iterator elem_iter = first_neighbour_indices.begin();
                                                elem_iter != first_neighbour_indices.end();
                                                ++elem_iter)
                                        {
                                            VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(*elem_iter);
                                            unsigned num_nodes = p_element->GetNumNodes();

                                            for (unsigned local_node_index = 0;
                                                    local_node_index < p_element->GetNumNodes();
                                                    ++local_node_index)
                                            {
                                                // first_neighbour_node_indices.insert(p_element->GetNodeGlobalIndex(local_node_index));
                                                unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
                                                unsigned global_node_neigh = p_element->GetNodeGlobalIndex((local_node_index+1)%num_nodes);

                                                if(global_node_index != node_B_index && global_node_index != node_A_index )
                                                {
                                                    if(global_node_neigh == node_C_index)
                                                    {
                                                        Node<DIM>* p_node = p_mesh->GetNode(global_node_index);
                                                        std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
                                                        if(containing_elements.size() == 1)
                                                        {
                                                            if(p_node->IsBoundaryNode())
                                                            {
                                                                is_boundary = true;
                                                                break;
                                                            }
                                                        }
                                                    }
                                                }
                                                
                                            }
                                        }
                                    }

                                    p_mesh->PerformNodeMerge(p_node, p_boundary_1);

                                    Node<DIM>* p_merged_node = p_boundary_1;

                                    if (p_boundary_1->IsDeleted())
                                    {
                                        p_merged_node = p_node;
                                    }

                                    p_mesh->PerformNodeMerge(p_boundary_2, p_merged_node);

                                    if (p_merged_node->IsDeleted())
                                    {
                                        p_merged_node = p_boundary_2;
                                    }

                                    p_merged_node->rGetModifiableLocation() = nodes_midpoint;

                                    // Tag remaining node as non-boundary, or boundary
                                    p_merged_node->SetAsBoundaryNode(is_boundary);

                                    PRINT_VARIABLE(is_boundary);

                                    performed_edge_modifier = true;
                                    TRACE("Performed Void Removel");
                                }
                            }
                        }
                        else if(p_neighbour_node_1->IsBoundaryNode() && p_neighbour_node_3->IsBoundaryNode() && IsNodeABNeighbours(node_neighbour_1,node_neighbour_3, rCellPopulation ))
                        {
                            std::set<unsigned> containing_elements_1 = p_neighbour_node_1->rGetContainingElementIndices();
                            std::set<unsigned> containing_elements_2 = p_neighbour_node_3->rGetContainingElementIndices();

                            std::set<unsigned> shared_elements;
                            std::set_intersection(containing_elements_1.begin(),
                                            containing_elements_1.end(),
                                            containing_elements_2.begin(),
                                            containing_elements_2.end(),
                                            std::inserter(shared_elements, shared_elements.begin()));
                            if(containing_element_indices.size()==2 && containing_elements_2.size()==2 && shared_elements.size()==1)
                            {
                                Node<DIM>* p_boundary_1 = p_neighbour_node_1;
                                Node<DIM>* p_boundary_2 = p_neighbour_node_3;
                                
                                Node<DIM>* p_node = p_mesh->GetNode(node_index);
                                c_vector<double, DIM> r_node = p_node->rGetLocation();
                                c_vector<double, DIM> r_boundary_1 = p_boundary_1->rGetLocation();
                                c_vector<double, DIM> r_boundary_2 = p_boundary_2->rGetLocation();

                                double a = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_1));
                                double b = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_2));
                                double c = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_1, r_boundary_2));
                                double s = 0.5*(a+b+c);

                                double void_area = sqrt(s*(s-a)*(s-b)*(s-c));

                                if(void_area < mVoidAreaThreshold)
                                {
                                    // p_mesh->PerformVoidRemoval(p_node,p_boundary_1,p_boundary_2);
                                    c_vector<double, DIM> nodes_midpoint = p_node->rGetLocation()
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_1->rGetLocation()) / 3.0
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_2->rGetLocation()) / 3.0;


                                    PRINT_VECTOR(nodes_midpoint);

                                    /*          \             /
                                    *            \   (Cell)  /
                                    *             o---------o
                                    *              \       /
                                    *               \ (v) /
                                    *                \   /
                                    *     (Cell)       o     (Cell)
                                    *                /   \
                                    *               /     \
                                    *         -----o  (v)  o-----
                                    */
                                    // Decide whether the last merged node will be a boundary node or not...  CBFed right now...
                                    bool is_boundary = false;
                                    // Check the neighbours of A that are not B or C and check whether any of those are boundary nodes.
                                    // If they are a boundary node, then the merged void will be a boundary node. Otherwise, not!
                                    
                                    unsigned node_A_index = p_node->GetIndex();
                                    unsigned node_B_index = p_boundary_1->GetIndex();
                                    unsigned node_C_index = p_boundary_2->GetIndex();

                                    // Check As neighbours
                                    std::set<unsigned> first_neighbour_indices = p_node->rGetContainingElementIndices();
                                    for (std::set<unsigned>::iterator elem_iter = first_neighbour_indices.begin();
                                            elem_iter != first_neighbour_indices.end();
                                            ++elem_iter)
                                    {
                                        VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(*elem_iter);

                                        for (unsigned local_node_index = 0;
                                                local_node_index < p_element->GetNumNodes();
                                                ++local_node_index)
                                        {
                                            unsigned num_nodes = p_element->GetNumNodes();

                                            unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
                                            unsigned global_node_neigh = p_element->GetNodeGlobalIndex((local_node_index+1)%num_nodes);
                                            if(global_node_index != node_B_index && global_node_index != node_C_index )
                                            {
                                                if(global_node_neigh == node_A_index)
                                                {
                                                    Node<DIM>* p_node = p_mesh->GetNode(global_node_index);
                                                    std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
                                                    if(containing_elements.size() == 1)
                                                    {
                                                        if(p_node->IsBoundaryNode())
                                                        {   
                                                            is_boundary = true;
                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                            
                                        }
                                    }
                                    // Check Bs neighbours
                                    if(!is_boundary)
                                    {
                                        std::set<unsigned> first_neighbour_indices = p_boundary_1->rGetContainingElementIndices();
                                        for (std::set<unsigned>::iterator elem_iter = first_neighbour_indices.begin();
                                                elem_iter != first_neighbour_indices.end();
                                                ++elem_iter)
                                        {
                                            VertexElement<DIM, DIM>*  p_element = p_mesh->GetElement(*elem_iter);
                                            unsigned num_nodes = p_element->GetNumNodes();

                                            for (unsigned local_node_index = 0;
                                                    local_node_index < p_element->GetNumNodes();
                                                    ++local_node_index)
                                            {
                                                // first_neighbour_node_indices.insert(p_element->GetNodeGlobalIndex(local_node_index));
                                                unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
                                                unsigned global_node_neigh = p_element->GetNodeGlobalIndex((local_node_index+1)%num_nodes);
                                                if(global_node_index != node_A_index && global_node_index != node_C_index )
                                                {
                                                    if(global_node_neigh == node_B_index)
                                                    {
                                                        Node<DIM>* p_node = p_mesh->GetNode(global_node_index);
                                                        std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
                                                        if(containing_elements.size() == 1)
                                                        {
                                                            if(p_node->IsBoundaryNode())
                                                            {
                                                                is_boundary = true;
                                                                break;
                                                            }
                                                        }
                                                    }
                                                }
                                                
                                            }
                                        }
                                    }
                                    // Check Cs neighbours
                                    if(!is_boundary)
                                    {
                                        std::set<unsigned> first_neighbour_indices = p_boundary_2->rGetContainingElementIndices();
                                        for (std::set<unsigned>::iterator elem_iter = first_neighbour_indices.begin();
                                                elem_iter != first_neighbour_indices.end();
                                                ++elem_iter)
                                        {
                                            VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(*elem_iter);
                                            unsigned num_nodes = p_element->GetNumNodes();

                                            for (unsigned local_node_index = 0;
                                                    local_node_index < p_element->GetNumNodes();
                                                    ++local_node_index)
                                            {
                                                // first_neighbour_node_indices.insert(p_element->GetNodeGlobalIndex(local_node_index));
                                                unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
                                                unsigned global_node_neigh = p_element->GetNodeGlobalIndex((local_node_index+1)%num_nodes);

                                                if(global_node_index != node_B_index && global_node_index != node_A_index )
                                                {
                                                    if(global_node_neigh == node_C_index)
                                                    {
                                                        Node<DIM>* p_node = p_mesh->GetNode(global_node_index);
                                                        std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
                                                        if(containing_elements.size() == 1)
                                                        {
                                                            if(p_node->IsBoundaryNode())
                                                            {
                                                                is_boundary = true;
                                                                break;
                                                            }
                                                        }
                                                    }
                                                }
                                                
                                            }
                                        }
                                    }

                                    p_mesh->PerformNodeMerge(p_node, p_boundary_1);

                                    Node<DIM>* p_merged_node = p_boundary_1;

                                    if (p_boundary_1->IsDeleted())
                                    {
                                        p_merged_node = p_node;
                                    }

                                    p_mesh->PerformNodeMerge(p_boundary_2, p_merged_node);

                                    if (p_merged_node->IsDeleted())
                                    {
                                        p_merged_node = p_boundary_2;
                                    }

                                    p_merged_node->rGetModifiableLocation() = nodes_midpoint;

                                    // Tag remaining node as non-boundary, or boundary
                                    p_merged_node->SetAsBoundaryNode(is_boundary);

                                    PRINT_VARIABLE(is_boundary);

                                    performed_edge_modifier = true;
                                    TRACE("Performed Void Removel");
                                }
                            }
                        }


                        if(p_neighbour_node_1->IsBoundaryNode() && p_neighbour_node_2->IsBoundaryNode() && !(IsNodeABNeighbours(node_neighbour_1,node_neighbour_2, rCellPopulation)) )
                        {
                            std::set<unsigned> node_neighbours_1 = p_cell_population->GetNeighbouringNodeIndices(node_neighbour_1);
                            std::set<unsigned> node_neighbours_2 = p_cell_population->GetNeighbouringNodeIndices(node_neighbour_2);

                            // Find common elements
                            std::set<unsigned> shared_elements;
                            std::set_intersection(node_neighbours_1.begin(),
                                                node_neighbours_1.end(),
                                                node_neighbours_2.begin(),
                                                node_neighbours_2.end(),
                                                std::inserter(shared_elements, shared_elements.begin()));

                            if(shared_elements.size() == 2)
                            {
                                std::set<unsigned>::const_iterator shared_neigh_it = shared_elements.begin();
                                unsigned shared_neigh_1 = (*shared_neigh_it);
                                shared_neigh_it++;
                                unsigned shared_neigh_2 = (*shared_neigh_it);

                                Node<DIM>* p_tmp_1 = p_mesh->GetNode(shared_neigh_1);
                                Node<DIM>* p_tmp_2 = p_mesh->GetNode(shared_neigh_2);

                                if(node_index == shared_neigh_1 && p_tmp_2->IsBoundaryNode())
                                {
                                    Node<DIM>* p_boundary_3 = p_mesh->GetNode(shared_neigh_2);
                                    Node<DIM>* p_boundary_1 = p_neighbour_node_1;
                                    Node<DIM>* p_boundary_2 = p_neighbour_node_2;
                                    
                                    Node<DIM>* p_node = p_mesh->GetNode(node_index);
                                    c_vector<double, DIM> r_node = p_node->rGetLocation();
                                    c_vector<double, DIM> r_boundary_1 = p_boundary_1->rGetLocation();
                                    c_vector<double, DIM> r_boundary_2 = p_boundary_2->rGetLocation();
                                    c_vector<double, DIM> r_boundary_3 = p_boundary_3->rGetLocation();

                                    double a1 = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_1));
                                    double b1 = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_2));
                                    double c1 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_1, r_boundary_2));
                                    double s1 = 0.5*(a1+b1+c1);

                                    double void_area_1 = sqrt(s1*(s1-a1)*(s1-b1)*(s1-c1));

                                    double a2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_3, r_boundary_1));
                                    double b2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_3, r_boundary_2));
                                    double c2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_1, r_boundary_2));
                                    double s2 = 0.5*(a2+b2+c2);

                                    double void_area_2 = sqrt(s2*(s2-a2)*(s2-b2)*(s2-c2));

                                    if(void_area_1+void_area_2 < mVoidAreaThreshold)
                                    {
                                        c_vector<double, DIM> nodes_midpoint = p_node->rGetLocation()
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_1->rGetLocation()) / 4.0
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_2->rGetLocation()) / 4.0
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_3->rGetLocation()) / 4.0;

                                        p_mesh->PerformNodeMerge(p_node, p_boundary_1);

                                        Node<DIM>* p_merged_node = p_boundary_1;

                                        if (p_boundary_1->IsDeleted())
                                        {
                                            p_merged_node = p_node;
                                        }
                                        p_mesh->PerformNodeMerge(p_boundary_2, p_merged_node);
                                        if (p_merged_node->IsDeleted())
                                        {
                                            p_merged_node = p_boundary_2;
                                        }
                                        p_mesh->PerformNodeMerge(p_boundary_3, p_merged_node);
                                        if (p_merged_node->IsDeleted())
                                        {
                                            p_merged_node = p_boundary_3;
                                        }

                                        p_merged_node->rGetModifiableLocation() = nodes_midpoint;

                                        // Tag remaining node as non-boundary, may need to check this.... hmmm
                                        p_merged_node->SetAsBoundaryNode(false);

                                        performed_edge_modifier = true;
                                        TRACE("Performed Void 2 Removel");

                                    }
                                    
                                }
                                else if(node_index == shared_neigh_2  && p_tmp_1->IsBoundaryNode())
                                {
                                    Node<DIM>* p_boundary_3 = p_mesh->GetNode(shared_neigh_1);
                                    Node<DIM>* p_boundary_1 = p_neighbour_node_1;
                                    Node<DIM>* p_boundary_2 = p_neighbour_node_2;
                                    
                                    Node<DIM>* p_node = p_mesh->GetNode(node_index);
                                    c_vector<double, DIM> r_node = p_node->rGetLocation();
                                    c_vector<double, DIM> r_boundary_1 = p_boundary_1->rGetLocation();
                                    c_vector<double, DIM> r_boundary_2 = p_boundary_2->rGetLocation();
                                    c_vector<double, DIM> r_boundary_3 = p_boundary_3->rGetLocation();

                                    double a1 = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_1));
                                    double b1 = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_2));
                                    double c1 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_1, r_boundary_2));
                                    double s1 = 0.5*(a1+b1+c1);

                                    double void_area_1 = sqrt(s1*(s1-a1)*(s1-b1)*(s1-c1));

                                    double a2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_3, r_boundary_1));
                                    double b2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_3, r_boundary_2));
                                    double c2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_1, r_boundary_2));
                                    double s2 = 0.5*(a2+b2+c2);

                                    double void_area_2 = sqrt(s2*(s2-a2)*(s2-b2)*(s2-c2));

                                    if(void_area_1+void_area_2 < mVoidAreaThreshold)
                                    {
                                        c_vector<double, DIM> nodes_midpoint = p_node->rGetLocation()
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_1->rGetLocation()) / 4.0
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_2->rGetLocation()) / 4.0
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_3->rGetLocation()) / 4.0;

                                        p_mesh->PerformNodeMerge(p_node, p_boundary_1);

                                        Node<DIM>* p_merged_node = p_boundary_1;

                                        if (p_boundary_1->IsDeleted())
                                        {
                                            p_merged_node = p_node;
                                        }
                                        p_mesh->PerformNodeMerge(p_boundary_2, p_merged_node);
                                        if (p_merged_node->IsDeleted())
                                        {
                                            p_merged_node = p_boundary_2;
                                        }
                                        p_mesh->PerformNodeMerge(p_boundary_3, p_merged_node);
                                        if (p_merged_node->IsDeleted())
                                        {
                                            p_merged_node = p_boundary_3;
                                        }

                                        p_merged_node->rGetModifiableLocation() = nodes_midpoint;

                                        // Tag remaining node as non-boundary, may need to check this.... hmmm
                                        p_merged_node->SetAsBoundaryNode(false);

                                        performed_edge_modifier = true;
                                        TRACE("Performed Void 2 Removel");

                                    }

                                }
                            }
                        }
                        else if(p_neighbour_node_1->IsBoundaryNode() && p_neighbour_node_3->IsBoundaryNode() && !(IsNodeABNeighbours(node_neighbour_1,node_neighbour_3, rCellPopulation)) )
                        {
                            std::set<unsigned> node_neighbours_1 = p_cell_population->GetNeighbouringNodeIndices(node_neighbour_1);
                            std::set<unsigned> node_neighbours_2 = p_cell_population->GetNeighbouringNodeIndices(node_neighbour_3);

                            // Find common elements
                            std::set<unsigned> shared_elements;
                            std::set_intersection(node_neighbours_1.begin(),
                                                node_neighbours_1.end(),
                                                node_neighbours_2.begin(),
                                                node_neighbours_2.end(),
                                                std::inserter(shared_elements, shared_elements.begin()));

                            if(shared_elements.size() == 2)
                            {
                                std::set<unsigned>::const_iterator shared_neigh_it = shared_elements.begin();
                                unsigned shared_neigh_1 = (*shared_neigh_it);
                                shared_neigh_it++;
                                unsigned shared_neigh_2 = (*shared_neigh_it);

                                Node<DIM>* p_tmp_1 = p_mesh->GetNode(shared_neigh_1);
                                Node<DIM>* p_tmp_2 = p_mesh->GetNode(shared_neigh_2);

                                if(node_index == shared_neigh_1 && p_tmp_2->IsBoundaryNode())
                                {
                                    Node<DIM>* p_boundary_3 = p_mesh->GetNode(shared_neigh_2);
                                    Node<DIM>* p_boundary_1 = p_neighbour_node_1;
                                    Node<DIM>* p_boundary_2 = p_neighbour_node_2;
                                    
                                    Node<DIM>* p_node = p_mesh->GetNode(node_index);
                                    c_vector<double, DIM> r_node = p_node->rGetLocation();
                                    c_vector<double, DIM> r_boundary_1 = p_boundary_1->rGetLocation();
                                    c_vector<double, DIM> r_boundary_2 = p_boundary_2->rGetLocation();
                                    c_vector<double, DIM> r_boundary_3 = p_boundary_3->rGetLocation();

                                    double a1 = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_1));
                                    double b1 = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_2));
                                    double c1 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_1, r_boundary_2));
                                    double s1 = 0.5*(a1+b1+c1);

                                    double void_area_1 = sqrt(s1*(s1-a1)*(s1-b1)*(s1-c1));

                                    double a2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_3, r_boundary_1));
                                    double b2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_3, r_boundary_2));
                                    double c2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_1, r_boundary_2));
                                    double s2 = 0.5*(a2+b2+c2);

                                    double void_area_2 = sqrt(s2*(s2-a2)*(s2-b2)*(s2-c2));

                                    if(void_area_1+void_area_2 < mVoidAreaThreshold)
                                    {
                                        c_vector<double, DIM> nodes_midpoint = p_node->rGetLocation()
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_1->rGetLocation()) / 4.0
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_2->rGetLocation()) / 4.0
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_3->rGetLocation()) / 4.0;

                                        p_mesh->PerformNodeMerge(p_node, p_boundary_1);

                                        Node<DIM>* p_merged_node = p_boundary_1;

                                        if (p_boundary_1->IsDeleted())
                                        {
                                            p_merged_node = p_node;
                                        }
                                        p_mesh->PerformNodeMerge(p_boundary_2, p_merged_node);
                                        if (p_merged_node->IsDeleted())
                                        {
                                            p_merged_node = p_boundary_2;
                                        }
                                        p_mesh->PerformNodeMerge(p_boundary_3, p_merged_node);
                                        if (p_merged_node->IsDeleted())
                                        {
                                            p_merged_node = p_boundary_3;
                                        }

                                        p_merged_node->rGetModifiableLocation() = nodes_midpoint;

                                        // Tag remaining node as non-boundary, may need to check this.... hmmm
                                        p_merged_node->SetAsBoundaryNode(false);

                                        performed_edge_modifier = true;
                                        TRACE("Performed Void 2 Removel");

                                    }
                                    
                                }
                                else if(node_index == shared_neigh_2  && p_tmp_1->IsBoundaryNode())
                                {
                                    Node<DIM>* p_boundary_3 = p_mesh->GetNode(shared_neigh_1);
                                    Node<DIM>* p_boundary_1 = p_neighbour_node_1;
                                    Node<DIM>* p_boundary_2 = p_neighbour_node_2;
                                    
                                    Node<DIM>* p_node = p_mesh->GetNode(node_index);
                                    c_vector<double, DIM> r_node = p_node->rGetLocation();
                                    c_vector<double, DIM> r_boundary_1 = p_boundary_1->rGetLocation();
                                    c_vector<double, DIM> r_boundary_2 = p_boundary_2->rGetLocation();
                                    c_vector<double, DIM> r_boundary_3 = p_boundary_3->rGetLocation();

                                    double a1 = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_1));
                                    double b1 = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_2));
                                    double c1 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_1, r_boundary_2));
                                    double s1 = 0.5*(a1+b1+c1);

                                    double void_area_1 = sqrt(s1*(s1-a1)*(s1-b1)*(s1-c1));

                                    double a2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_3, r_boundary_1));
                                    double b2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_3, r_boundary_2));
                                    double c2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_1, r_boundary_2));
                                    double s2 = 0.5*(a2+b2+c2);

                                    double void_area_2 = sqrt(s2*(s2-a2)*(s2-b2)*(s2-c2));

                                    if(void_area_1+void_area_2 < mVoidAreaThreshold)
                                    {
                                        c_vector<double, DIM> nodes_midpoint = p_node->rGetLocation()
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_1->rGetLocation()) / 4.0
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_2->rGetLocation()) / 4.0
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_3->rGetLocation()) / 4.0;

                                        p_mesh->PerformNodeMerge(p_node, p_boundary_1);

                                        Node<DIM>* p_merged_node = p_boundary_1;

                                        if (p_boundary_1->IsDeleted())
                                        {
                                            p_merged_node = p_node;
                                        }
                                        p_mesh->PerformNodeMerge(p_boundary_2, p_merged_node);
                                        if (p_merged_node->IsDeleted())
                                        {
                                            p_merged_node = p_boundary_2;
                                        }
                                        p_mesh->PerformNodeMerge(p_boundary_3, p_merged_node);
                                        if (p_merged_node->IsDeleted())
                                        {
                                            p_merged_node = p_boundary_3;
                                        }

                                        p_merged_node->rGetModifiableLocation() = nodes_midpoint;

                                        // Tag remaining node as non-boundary, may need to check this.... hmmm
                                        p_merged_node->SetAsBoundaryNode(false);

                                        performed_edge_modifier = true;
                                        TRACE("Performed Void 2 Removel");

                                    }

                                }
                            }
                        }
                        else if(p_neighbour_node_3->IsBoundaryNode() && p_neighbour_node_2->IsBoundaryNode() && !(IsNodeABNeighbours(node_neighbour_3,node_neighbour_2, rCellPopulation)) )
                        {
                            std::set<unsigned> node_neighbours_1 = p_cell_population->GetNeighbouringNodeIndices(node_neighbour_3);
                            std::set<unsigned> node_neighbours_2 = p_cell_population->GetNeighbouringNodeIndices(node_neighbour_2);

                            // Find common elements
                            std::set<unsigned> shared_elements;
                            std::set_intersection(node_neighbours_1.begin(),
                                                node_neighbours_1.end(),
                                                node_neighbours_2.begin(),
                                                node_neighbours_2.end(),
                                                std::inserter(shared_elements, shared_elements.begin()));

                            if(shared_elements.size() == 2)
                            {
                                std::set<unsigned>::const_iterator shared_neigh_it = shared_elements.begin();
                                unsigned shared_neigh_1 = (*shared_neigh_it);
                                shared_neigh_it++;
                                unsigned shared_neigh_2 = (*shared_neigh_it);

                                Node<DIM>* p_tmp_1 = p_mesh->GetNode(shared_neigh_1);
                                Node<DIM>* p_tmp_2 = p_mesh->GetNode(shared_neigh_2);

                                if(node_index == shared_neigh_1 && p_tmp_2->IsBoundaryNode())
                                {
                                    Node<DIM>* p_boundary_3 = p_mesh->GetNode(shared_neigh_2);
                                    Node<DIM>* p_boundary_1 = p_neighbour_node_1;
                                    Node<DIM>* p_boundary_2 = p_neighbour_node_2;
                                    
                                    Node<DIM>* p_node = p_mesh->GetNode(node_index);
                                    c_vector<double, DIM> r_node = p_node->rGetLocation();
                                    c_vector<double, DIM> r_boundary_1 = p_boundary_1->rGetLocation();
                                    c_vector<double, DIM> r_boundary_2 = p_boundary_2->rGetLocation();
                                    c_vector<double, DIM> r_boundary_3 = p_boundary_3->rGetLocation();

                                    double a1 = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_1));
                                    double b1 = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_2));
                                    double c1 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_1, r_boundary_2));
                                    double s1 = 0.5*(a1+b1+c1);

                                    double void_area_1 = sqrt(s1*(s1-a1)*(s1-b1)*(s1-c1));

                                    double a2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_3, r_boundary_1));
                                    double b2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_3, r_boundary_2));
                                    double c2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_1, r_boundary_2));
                                    double s2 = 0.5*(a2+b2+c2);

                                    double void_area_2 = sqrt(s2*(s2-a2)*(s2-b2)*(s2-c2));

                                    if(void_area_1+void_area_2 < mVoidAreaThreshold)
                                    {
                                        c_vector<double, DIM> nodes_midpoint = p_node->rGetLocation()
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_1->rGetLocation()) / 4.0
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_2->rGetLocation()) / 4.0
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_3->rGetLocation()) / 4.0;

                                        p_mesh->PerformNodeMerge(p_node, p_boundary_1);

                                        Node<DIM>* p_merged_node = p_boundary_1;

                                        if (p_boundary_1->IsDeleted())
                                        {
                                            p_merged_node = p_node;
                                        }
                                        p_mesh->PerformNodeMerge(p_boundary_2, p_merged_node);
                                        if (p_merged_node->IsDeleted())
                                        {
                                            p_merged_node = p_boundary_2;
                                        }
                                        p_mesh->PerformNodeMerge(p_boundary_3, p_merged_node);
                                        if (p_merged_node->IsDeleted())
                                        {
                                            p_merged_node = p_boundary_3;
                                        }

                                        p_merged_node->rGetModifiableLocation() = nodes_midpoint;

                                        // Tag remaining node as non-boundary, may need to check this.... hmmm
                                        p_merged_node->SetAsBoundaryNode(false);

                                        performed_edge_modifier = true;
                                        TRACE("Performed Void 2 Removel");

                                    }
                                    
                                }
                                else if(node_index == shared_neigh_2  && p_tmp_1->IsBoundaryNode())
                                {
                                    Node<DIM>* p_boundary_3 = p_mesh->GetNode(shared_neigh_1);
                                    Node<DIM>* p_boundary_1 = p_neighbour_node_1;
                                    Node<DIM>* p_boundary_2 = p_neighbour_node_2;
                                    
                                    Node<DIM>* p_node = p_mesh->GetNode(node_index);
                                    c_vector<double, DIM> r_node = p_node->rGetLocation();
                                    c_vector<double, DIM> r_boundary_1 = p_boundary_1->rGetLocation();
                                    c_vector<double, DIM> r_boundary_2 = p_boundary_2->rGetLocation();
                                    c_vector<double, DIM> r_boundary_3 = p_boundary_3->rGetLocation();

                                    double a1 = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_1));
                                    double b1 = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_2));
                                    double c1 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_1, r_boundary_2));
                                    double s1 = 0.5*(a1+b1+c1);

                                    double void_area_1 = sqrt(s1*(s1-a1)*(s1-b1)*(s1-c1));

                                    double a2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_3, r_boundary_1));
                                    double b2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_3, r_boundary_2));
                                    double c2 = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_1, r_boundary_2));
                                    double s2 = 0.5*(a2+b2+c2);

                                    double void_area_2 = sqrt(s2*(s2-a2)*(s2-b2)*(s2-c2));

                                    if(void_area_1+void_area_2 < mVoidAreaThreshold)
                                    {
                                        c_vector<double, DIM> nodes_midpoint = p_node->rGetLocation()
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_1->rGetLocation()) / 4.0
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_2->rGetLocation()) / 4.0
                                        + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_3->rGetLocation()) / 4.0;

                                        p_mesh->PerformNodeMerge(p_node, p_boundary_1);

                                        Node<DIM>* p_merged_node = p_boundary_1;

                                        if (p_boundary_1->IsDeleted())
                                        {
                                            p_merged_node = p_node;
                                        }
                                        p_mesh->PerformNodeMerge(p_boundary_2, p_merged_node);
                                        if (p_merged_node->IsDeleted())
                                        {
                                            p_merged_node = p_boundary_2;
                                        }
                                        p_mesh->PerformNodeMerge(p_boundary_3, p_merged_node);
                                        if (p_merged_node->IsDeleted())
                                        {
                                            p_merged_node = p_boundary_3;
                                        }

                                        p_merged_node->rGetModifiableLocation() = nodes_midpoint;

                                        // Tag remaining node as non-boundary, may need to check this.... hmmm
                                        p_merged_node->SetAsBoundaryNode(false);

                                        performed_edge_modifier = true;
                                        TRACE("Performed Void 2 Removel");

                                    }

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
        p_mesh->ReMesh();
        TRACE("4th_mesh");
        unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
        std::stringstream time;
        time << num_timesteps;
        PRINT_VARIABLE(num_timesteps);
        VertexMeshWriter<DIM,DIM> vertexmesh_writer("tmp", "4th_mesh", false);
        vertexmesh_writer.WriteVtkUsingMesh(*p_mesh, time.str());
    }
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

