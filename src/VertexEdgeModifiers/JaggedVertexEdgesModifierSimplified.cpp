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

#include "JaggedVertexEdgesModifierSimplified.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "UblasCustomFunctions.hpp"

#include "VertexMeshWriter.hpp"
#include "MutableMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "Debug.hpp"

template<unsigned DIM>
JaggedVertexEdgesModifierSimplified<DIM>::JaggedVertexEdgesModifierSimplified()
    : AbstractCellBasedSimulationModifier<DIM>(),
    mMaxEdgeLength(0.15), //0.10
    mMinEdgeLength(0.05) //0.025
{
}

template<unsigned DIM>
JaggedVertexEdgesModifierSimplified<DIM>::~JaggedVertexEdgesModifierSimplified()
{
}

template<unsigned DIM>
void JaggedVertexEdgesModifierSimplified<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // TRACE("*****************************************************************************************************************************");
    // PRINT_VARIABLE(SimulationTime::Instance()->GetTime());
    unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    // PRINT_VARIABLE(num_timesteps);
    std::stringstream time;
    time << num_timesteps;

    bool recheck_mesh = true;
    bool recheck_edges = true;
    bool recheck_voids = true;
    bool recheck_node_labels = true;
    bool recheck_int_nodes = true;
    bool recheck_node_edge = true;
    bool recheck_hack = true;
    int how_many_times_here = 0;
    while(recheck_hack || recheck_mesh || recheck_edges || recheck_voids || recheck_node_labels || recheck_int_nodes || recheck_node_edge)
    {   
        // how_many_times_here = 0;
        VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));

        
        std::stringstream t_numb_times_here;
        t_numb_times_here << how_many_times_here;


        // recheck_hack = DirtyHack(rCellPopulation,how_many_times_here);
        recheck_hack = ClearFreeInternalNodes(rCellPopulation,how_many_times_here);

        // VertexBasedCellPopulation<DIM>* p_cell_population1 = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        // MutableVertexMesh<DIM,DIM>* p_mesh1 = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population1->rGetMesh()));
        // VertexMeshWriter<DIM,DIM> vertexmesh_writer1("tmp", "Mesh_1", false);
        // vertexmesh_writer1.WriteVtkUsingMesh(*p_mesh1, time.str() + "_" + t_numb_times_here.str());

        while(false)
        {
            recheck_edges = RefineEdges(rCellPopulation,how_many_times_here);
            how_many_times_here++;
        }
        recheck_edges = false;
        // VertexBasedCellPopulation<DIM>* p_cell_population2 = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        // MutableVertexMesh<DIM,DIM>* p_mesh2 = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population2->rGetMesh()));
        // VertexMeshWriter<DIM,DIM> vertexmesh_writer2("tmp", "Mesh_2", false);
        // vertexmesh_writer2.WriteVtkUsingMesh(*p_mesh2, time.str() + "_" + t_numb_times_here.str());

        // TRACE("Refined Edges");

        while(false)
        {
            recheck_mesh = CheckForVConfig(rCellPopulation,how_many_times_here);
            how_many_times_here++;
        }
        recheck_mesh = false;
        // VertexBasedCellPopulation<DIM>* p_cell_population3 = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        // MutableVertexMesh<DIM,DIM>* p_mesh3 = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population3->rGetMesh()));
        // VertexMeshWriter<DIM,DIM> vertexmesh_writer3("tmp", "Mesh_3", false);
        // vertexmesh_writer3.WriteVtkUsingMesh(*p_mesh3, time.str() + "_" + t_numb_times_here.str());
        // TRACE("Checked for V");

        while(recheck_voids)
        {
            recheck_voids = CheckFor3NodeVoids(rCellPopulation,how_many_times_here);
            how_many_times_here++;
        }
        // VertexBasedCellPopulation<DIM>* p_cell_population4 = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        // MutableVertexMesh<DIM,DIM>* p_mesh4 = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population4->rGetMesh()));
        // VertexMeshWriter<DIM,DIM> vertexmesh_writer4("tmp", "Mesh_4", false);
        // vertexmesh_writer4.WriteVtkUsingMesh(*p_mesh4, time.str() + "_" + t_numb_times_here.str());
    //    TRACE("Checked Voids");

       while(recheck_node_edge)
        {
            recheck_node_edge = CheckForNodeEdgeInt(rCellPopulation,how_many_times_here);
            how_many_times_here++;
        }
        // VertexBasedCellPopulation<DIM>* p_cell_population5 = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        // MutableVertexMesh<DIM,DIM>* p_mesh5 = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population5->rGetMesh()));
        // VertexMeshWriter<DIM,DIM> vertexmesh_writer5("tmp", "Mesh_5", false);
        // vertexmesh_writer5.WriteVtkUsingMesh(*p_mesh5, time.str() + "_" + t_numb_times_here.str());
    //    TRACE("checked Node Edge");

       while(recheck_int_nodes)
        {
            recheck_int_nodes = CheckForInternalNode(rCellPopulation,how_many_times_here);
            how_many_times_here++;
        }
        // VertexBasedCellPopulation<DIM>* p_cell_population6 = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        // MutableVertexMesh<DIM,DIM>* p_mesh6 = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population6->rGetMesh()));
        // VertexMeshWriter<DIM,DIM> vertexmesh_writer6("tmp", "Mesh_6", false);
        // vertexmesh_writer6.WriteVtkUsingMesh(*p_mesh6, time.str() + "_" + t_numb_times_here.str());
    //    TRACE("Internal Nodes");

    //    while(recheck_node_labels)
    //     {
            recheck_node_labels = CleanUpRogueNodes(rCellPopulation,how_many_times_here);
    //     }
    //    TRACE("Cleaned Rogues");
    //    recheck_node_labels = true;
        // PRINT_VARIABLE(how_many_times_here);
        // PRINT_3_VARIABLES(recheck_mesh,recheck_edges,recheck_voids);
        // PRINT_3_VARIABLES(recheck_node_labels,recheck_int_nodes,recheck_node_edge);



        // recheck_edges = RefineEdges(rCellPopulation,how_many_times_here);
        // TRACE("Refined Edges");
        // VertexMeshWriter<DIM,DIM> vertexmesh_writer1("tmp", "Mesh_1_", false);
        // vertexmesh_writer1.WriteVtkUsingMesh(*p_mesh, time.str() + "_" + t_numb_times_here.str());

        // recheck_mesh = CheckForVConfig(rCellPopulation,how_many_times_here);
        // TRACE("Checked for V");
        // VertexMeshWriter<DIM,DIM> vertexmesh_writer2("tmp", "Mesh_2_", false);
        // vertexmesh_writer2.WriteVtkUsingMesh(*p_mesh, time.str() + "_" + t_numb_times_here.str());

        // recheck_voids = CheckFor3NodeVoids(rCellPopulation,how_many_times_here);
        // TRACE("Checked Voids");
        // VertexMeshWriter<DIM,DIM> vertexmesh_writer3("tmp", "Mesh_3_", false);
        // vertexmesh_writer3.WriteVtkUsingMesh(*p_mesh, time.str() + "_" + t_numb_times_here.str());

        // recheck_node_edge = CheckForNodeEdgeInt(rCellPopulation,how_many_times_here);
        // TRACE("checked Node Edge");
        // VertexMeshWriter<DIM,DIM> vertexmesh_writer5("tmp", "Mesh_4_", false);
        // vertexmesh_writer5.WriteVtkUsingMesh(*p_mesh, time.str() + "_" + t_numb_times_here.str());

        // recheck_int_nodes = CheckForInternalNode(rCellPopulation,how_many_times_here);
        // TRACE("Internal Nodes");
        // VertexMeshWriter<DIM,DIM> vertexmesh_writer4("tmp", "Mesh_5_", false);
        // vertexmesh_writer4.WriteVtkUsingMesh(*p_mesh, time.str() + "_" + t_numb_times_here.str());

        // recheck_node_labels = CleanUpRogueNodes(rCellPopulation,how_many_times_here);
        // // recheck_node_labels = false;
        // TRACE("Cleaned Rogues");
        // VertexMeshWriter<DIM,DIM> vertexmesh_writer6("tmp", "Mesh_6_", false);
        // vertexmesh_writer6.WriteVtkUsingMesh(*p_mesh, time.str() + "_" + t_numb_times_here.str());

        // PRINT_3_VARIABLES(how_many_times_here, recheck_edges, recheck_mesh);
        // PRINT_VARIABLE(how_many_times_here);
        // PRINT_3_VARIABLES(recheck_edges,recheck_mesh,recheck_voids);
        // PRINT_3_VARIABLES(recheck_node_edge,recheck_int_nodes,recheck_node_labels);
        
        // Ensure we dont get stuck here for ever... (Happens in growing monolayer...)
        // how_many_times_here++;
        
        if(how_many_times_here > 15)
        {
            // PRINT_3_VARIABLES(recheck_mesh,recheck_edges,recheck_voids);
            // PRINT_3_VARIABLES(recheck_node_labels,recheck_int_nodes,recheck_node_edge);
            recheck_mesh= false;
            recheck_edges = false;
            recheck_voids = false;
            recheck_node_labels = false;
            recheck_int_nodes = false;
            recheck_node_edge = false;
            recheck_hack  = false;
            // TRACE("Exited the loop");
            // PRINT_VARIABLE(SimulationTime::Instance()->GetTime());
        }
        // TRACE("*******************");

    }
    
}

template<unsigned DIM>
void JaggedVertexEdgesModifierSimplified<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    
    bool recheck_mesh = true;
    bool recheck_edges = true;

}

template<unsigned DIM>
double JaggedVertexEdgesModifierSimplified<DIM>::DistanceToEdgeFromPoint( c_vector<double, DIM> P1,  c_vector<double, DIM> P2,  c_vector<double, DIM> P0)
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
bool JaggedVertexEdgesModifierSimplified<DIM>::IsNodeABNeighbours(unsigned node_a,unsigned node_b, AbstractCellPopulation<DIM,DIM>& rCellPopulation )
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
bool JaggedVertexEdgesModifierSimplified<DIM>::Remove3NodeVoid(bool performed_edge_mod,unsigned node_index, Node<DIM>* p_neighbour_node_1, Node<DIM>* p_neighbour_node_2, VertexBasedCellPopulation<DIM>* p_cell_population, MutableVertexMesh<DIM,DIM>* p_mesh)
{
    bool performed_edge_modifier = performed_edge_mod;

    double mVoidAreaThreshold = 0.1;

    // VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    // MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));

    Node<DIM>* p_node = p_mesh->GetNode(node_index);

    std::set<unsigned> containing_element_indices = p_node->rGetContainingElementIndices();

    std::set<unsigned> containing_elements_1 = p_neighbour_node_1->rGetContainingElementIndices();
    std::set<unsigned> containing_elements_2 = p_neighbour_node_2->rGetContainingElementIndices();

    // std::set<unsigned> shared_elements;
    // std::set_intersection(containing_elements_1.begin(),
    //     containing_elements_1.end(),
    //     containing_elements_2.begin(),
    //     containing_elements_2.end(),
    //     std::inserter(shared_elements, shared_elements.begin()));
    // if(containing_elements_1.size()==2 && containing_elements_2.size()==2 && shared_elements.size()==1)
    // {
        Node<DIM>* p_boundary_1 = p_neighbour_node_1;
        Node<DIM>* p_boundary_2 = p_neighbour_node_2;
                    
        c_vector<double, DIM> r_node = p_node->rGetLocation();
        c_vector<double, DIM> r_boundary_1 = p_boundary_1->rGetLocation();
        c_vector<double, DIM> r_boundary_2 = p_boundary_2->rGetLocation();

        PRINT_VECTOR(r_node);
        PRINT_VECTOR(r_boundary_1);
        PRINT_VECTOR(r_boundary_2);

        double a = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_1));
        double b = norm_2(p_mesh->GetVectorFromAtoB(r_node, r_boundary_2));
        double c = norm_2(p_mesh->GetVectorFromAtoB(r_boundary_1, r_boundary_2));
        double s = 0.5*(a+b+c);

        double void_area = sqrt(s*(s-a)*(s-b)*(s-c));

        if(void_area < mVoidAreaThreshold)
        {
            c_vector<double, DIM> nodes_midpoint = p_node->rGetLocation()
            + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_1->rGetLocation()) / 3.0
            + p_mesh->GetVectorFromAtoB(p_node->rGetLocation(), p_boundary_2->rGetLocation()) / 3.0;

            // PRINT_VECTOR(nodes_midpoint);

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
            // std::set<unsigned> first_neighbour_indices = p_node->rGetContainingElementIndices();
            // for (std::set<unsigned>::iterator elem_iter = first_neighbour_indices.begin();
            //         elem_iter != first_neighbour_indices.end();
            //         ++elem_iter)
            // {
            //     VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(*elem_iter);

            //     for (unsigned local_node_index = 0;
            //             local_node_index < p_element->GetNumNodes();
            //             ++local_node_index)
            //         {
            //             unsigned num_nodes = p_element->GetNumNodes();

            //             unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
            //             unsigned global_node_neigh = p_element->GetNodeGlobalIndex((local_node_index+1)%num_nodes);
            //             if(global_node_index != node_B_index && global_node_index != node_C_index )
            //             {
            //                 if(global_node_neigh == node_A_index)
            //                 {
            //                     Node<DIM>* p_node = p_mesh->GetNode(global_node_index);
            //                     std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
            //                     if(containing_elements.size() == 1)
            //                     {
            //                         if(p_node->IsBoundaryNode())
            //                         {   
            //                             is_boundary = true;
            //                             break;
            //                         }
            //                     }
            //                 }
            //             }
                                            
            //         }
            // }
            // // Check Bs neighbours
            // if(!is_boundary)
            // {
            //     std::set<unsigned> first_neighbour_indices = p_boundary_1->rGetContainingElementIndices();
            //     for (std::set<unsigned>::iterator elem_iter = first_neighbour_indices.begin();
            //         elem_iter != first_neighbour_indices.end();
            //         ++elem_iter)
            //     {
            //         VertexElement<DIM, DIM>*  p_element = p_mesh->GetElement(*elem_iter);
            //         unsigned num_nodes = p_element->GetNumNodes();

            //         for (unsigned local_node_index = 0;
            //             local_node_index < p_element->GetNumNodes();
            //             ++local_node_index)
            //         {
            //             // first_neighbour_node_indices.insert(p_element->GetNodeGlobalIndex(local_node_index));
            //             unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
            //             unsigned global_node_neigh = p_element->GetNodeGlobalIndex((local_node_index+1)%num_nodes);
            //             if(global_node_index != node_A_index && global_node_index != node_C_index )
            //             {
            //                 if(global_node_neigh == node_B_index)
            //                 {
            //                     Node<DIM>* p_node = p_mesh->GetNode(global_node_index);
            //                     std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
            //                     if(containing_elements.size() == 1)
            //                     {
            //                         if(p_node->IsBoundaryNode())
            //                         {
            //                             is_boundary = true;
            //                             break;
            //                         }
            //                     }
            //                 }
            //             }
                                                
            //         }
            //     }
            // }
            // // Check Cs neighbours
            // if(!is_boundary)
            // {
            //     std::set<unsigned> first_neighbour_indices = p_boundary_2->rGetContainingElementIndices();
            //     for (std::set<unsigned>::iterator elem_iter = first_neighbour_indices.begin();
            //         elem_iter != first_neighbour_indices.end();
            //         ++elem_iter)
            //     {
            //         VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(*elem_iter);
            //         unsigned num_nodes = p_element->GetNumNodes();

            //         for (unsigned local_node_index = 0;
            //             local_node_index < p_element->GetNumNodes();
            //             ++local_node_index)
            //         {
            //             // first_neighbour_node_indices.insert(p_element->GetNodeGlobalIndex(local_node_index));
            //             unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
            //             unsigned global_node_neigh = p_element->GetNodeGlobalIndex((local_node_index+1)%num_nodes);

            //             if(global_node_index != node_B_index && global_node_index != node_A_index )
            //             {
            //                 if(global_node_neigh == node_C_index)
            //                 {
            //                     Node<DIM>* p_node = p_mesh->GetNode(global_node_index);
            //                     std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
            //                     if(containing_elements.size() == 1)
            //                     {
            //                         if(p_node->IsBoundaryNode())
            //                         {
            //                             is_boundary = true;
            //                             break;
            //                         }
            //                     }
            //                 }
            //             }
                                                
            //         }
            //     }
            // }

            p_mesh->PerformNodeMerge(p_node, p_boundary_1);
            p_node->SetAsBoundaryNode(true);
            Node<DIM>* p_merged_node = p_boundary_1;

            if (p_boundary_1->IsDeleted())
            {
                p_merged_node = p_node;
            }
            p_mesh->PerformNodeMerge(p_boundary_2, p_merged_node);
            p_boundary_2->SetAsBoundaryNode(true);
            if (p_merged_node->IsDeleted())
            {
                p_merged_node = p_boundary_2;
            }
            p_merged_node->rGetModifiableLocation() = nodes_midpoint;
            p_merged_node->SetAsBoundaryNode(true);
            
            // p_mesh->ReMesh();

            // unsigned numb_boundary_neighs = 0;
            // if(is_boundary==false)
            // {
            //     unsigned merged_node_index = p_merged_node->GetIndex();
            //     std::set<unsigned> merged_node_neighbours = p_cell_population->GetNeighbouringNodeIndices(merged_node_index);
                
            //     for (std::set<unsigned>::iterator neighbour_iter = merged_node_neighbours.begin();
            //         neighbour_iter != merged_node_neighbours.end();
            //         ++neighbour_iter)
            //     {
            //         unsigned neighbour_index = *neighbour_iter;
            //         Node<DIM>* p_neighbour_node = p_mesh->GetNode(neighbour_index);
            //         std::set<unsigned> element_index_set_n = p_neighbour_node->rGetContainingElementIndices();

            //         if(p_neighbour_node->IsBoundaryNode() )
            //         {
            //             numb_boundary_neighs++;
            //         }

            //     }
                

            // }
            // if(numb_boundary_neighs >= 2)
            // {
            //     is_boundary = true;
            // }


            // // Tag remaining node as non-boundary, or boundary
            // // p_merged_node->SetAsBoundaryNode(is_boundary);
            // p_merged_node->SetAsBoundaryNode(true);

            performed_edge_modifier = true;
            
            
            return performed_edge_modifier;
        }
    // }

    return performed_edge_modifier;

}

template<unsigned DIM>
bool JaggedVertexEdgesModifierSimplified<DIM>::Remove4NodeVoid(bool performed_edge_mod,unsigned node_index, unsigned node_neighbour_1, unsigned node_neighbour_2, VertexBasedCellPopulation<DIM>* p_cell_population, MutableVertexMesh<DIM,DIM>* p_mesh)
{
    bool performed_edge_modifier = performed_edge_mod;

    double mVoidAreaThreshold = 0.05;

    // VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    // MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));

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
            Node<DIM>* p_boundary_1 = p_mesh->GetNode(node_neighbour_1);
            Node<DIM>* p_boundary_2 = p_mesh->GetNode(node_neighbour_2);
                                    
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


                // Decide whether the last merged node will be a boundary node or not...  CBFed right now...
                bool is_boundary = false;
                // Check the neighbours of A that are not B or C and check whether any of those are boundary nodes.
                // If they are a boundary node, then the merged void will be a boundary node. Otherwise, not!
                                        
                unsigned node_A_index = p_node->GetIndex();
                unsigned node_B_index = p_boundary_1->GetIndex();
                unsigned node_C_index = p_boundary_2->GetIndex();
                unsigned node_D_index = shared_neigh_2;

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
                            if(global_node_index != node_B_index && global_node_index != node_C_index && global_node_index != node_D_index)
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
                            if(global_node_index != node_A_index && global_node_index != node_C_index  && global_node_index != node_D_index)
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

                            if(global_node_index != node_B_index && global_node_index != node_A_index  && global_node_index != node_D_index)
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
                p_mesh->PerformNodeMerge(p_boundary_3, p_merged_node);
                if (p_merged_node->IsDeleted())
                {
                    p_merged_node = p_boundary_3;
                }

                p_merged_node->rGetModifiableLocation() = nodes_midpoint;

                // Tag remaining node as non-boundary, may need to check this.... hmmm
                p_merged_node->SetAsBoundaryNode(is_boundary);
                performed_edge_modifier = true;
                
                p_mesh->ReMesh();
                return performed_edge_modifier;
            }
                                    
        }
        else if(node_index == shared_neigh_2  && p_tmp_1->IsBoundaryNode())
        {
            Node<DIM>* p_boundary_3 = p_mesh->GetNode(shared_neigh_1);
            Node<DIM>* p_boundary_1 = p_mesh->GetNode(node_neighbour_1);
            Node<DIM>* p_boundary_2 = p_mesh->GetNode(node_neighbour_2);
                                    
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

                // Decide whether the last merged node will be a boundary node or not...  CBFed right now...
                bool is_boundary = false;
                // Check the neighbours of A that are not B or C and check whether any of those are boundary nodes.
                // If they are a boundary node, then the merged void will be a boundary node. Otherwise, not!
                                        
                unsigned node_A_index = p_node->GetIndex();
                unsigned node_B_index = p_boundary_1->GetIndex();
                unsigned node_C_index = p_boundary_2->GetIndex();
                unsigned node_D_index = shared_neigh_2;

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
                            if(global_node_index != node_B_index && global_node_index != node_C_index && global_node_index != node_D_index)
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
                            if(global_node_index != node_A_index && global_node_index != node_C_index  && global_node_index != node_D_index)
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

                            if(global_node_index != node_B_index && global_node_index != node_A_index  && global_node_index != node_D_index)
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
                p_mesh->PerformNodeMerge(p_boundary_3, p_merged_node);
                if (p_merged_node->IsDeleted())
                {
                    p_merged_node = p_boundary_3;
                }

                p_merged_node->rGetModifiableLocation() = nodes_midpoint;

                // Tag remaining node as non-boundary, may need to check this.... hmmm
                p_merged_node->SetAsBoundaryNode(is_boundary);
                performed_edge_modifier = true;

                p_mesh->ReMesh();
                return performed_edge_modifier;
            }

        }
    }
    return performed_edge_modifier;
}

template<unsigned DIM>
bool JaggedVertexEdgesModifierSimplified<DIM>::RefineEdges(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here)
{

    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("VertexBoundaryRefinementModifier is to be used with a VertexBasedCellPopulation only");
    }

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));
    
    // TRACE("Initial_mesh");
    //     unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    //     std::stringstream time;
    //     time << num_timesteps;
    //     std::stringstream t_numb_times_here;
    //     t_numb_times_here << numb_times_here;
    //     PRINT_VARIABLE(num_timesteps);
    //     VertexMeshWriter<DIM,DIM> vertexmesh_writer("tmp", "Initial_mesh", false);
    //     // vertexmesh_writer.WriteVtkUsingMesh(*p_mesh, time.str());
    //     vertexmesh_writer.WriteVtkUsingMesh(*p_mesh, time.str() + "_" + t_numb_times_here.str());



    bool recheck_edges = true;

    // while (recheck_edges)
    // {
        recheck_edges = false;

        for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_mesh->GetElementIteratorBegin();
            elem_iter != p_mesh->GetElementIteratorEnd();
            ++elem_iter)
        {

            for (unsigned node_local_index = 0; node_local_index < elem_iter->GetNumNodes(); node_local_index++)
            {
                unsigned next_node_local_index = (node_local_index+1) % (elem_iter->GetNumNodes());

                unsigned node_global_index = elem_iter->GetNodeGlobalIndex(node_local_index);
                unsigned next_node_global_index = elem_iter->GetNodeGlobalIndex(next_node_local_index);

                Node<DIM>* p_node_a = p_mesh->GetNode(node_global_index);
                Node<DIM>* p_node_b = p_mesh->GetNode(next_node_global_index);
                                
                if (p_node_a->IsBoundaryNode() && p_node_b->IsBoundaryNode())
                {
                    // Find the sets of elements containing nodes A and B
                    std::set<unsigned> node_a_elem_indices = p_node_a->rGetContainingElementIndices();
                    std::set<unsigned> node_b_elem_indices = p_node_b->rGetContainingElementIndices();

                    // Find common elements
                    std::set<unsigned> shared_elements;
                    std::set_intersection(node_a_elem_indices.begin(),
                            node_a_elem_indices.end(),
                            node_b_elem_indices.begin(),
                            node_b_elem_indices.end(),
                            std::inserter(shared_elements, shared_elements.begin()));

                    // assert(shared_elements.size()>0); //otherwise not in the same element at all 

                    if(shared_elements.size() == 1)
                    {
                        // Here we have a boundary edge so add new node if needed.
                        c_vector<double,DIM> edge = p_mesh->GetVectorFromAtoB(p_node_a->rGetLocation(),p_node_b->rGetLocation());
                        
                        if (norm_2(edge) > mMaxEdgeLength)
                        {
                            // TRACE("Divided an edge");
                            // PRINT_VECTOR(p_node_a->rGetLocation());
                            // PRINT_VECTOR(p_node_b->rGetLocation());
                            p_mesh->DivideEdge(p_node_a, p_node_b);
                            recheck_edges = true;
                            break;
                        }
                        else if (norm_2(edge) < mMinEdgeLength)
                        {
                            // TRACE("help");
                            // Check to make sure we don't delete a node that is a vertex, only a free boundary node
                            if(node_a_elem_indices.size() >= 2 && node_b_elem_indices.size() == 1)
                            {
                                // Delete neighbour b
                                // std::set<unsigned> element_index_set_b = p_node_b->rGetContainingElementIndices();
                                // unsigned elem_index_1 = (*element_index_set_b.begin());
                                // VertexElement<DIM,DIM>* p_element_1 = p_mesh->GetElement(elem_index_1);
                                // p_element_1->DeleteNode(p_element_1->GetNodeLocalIndex(next_node_global_index));
                                // p_mesh->DeleteNodePriorToReMesh(next_node_global_index);
                                
                                p_mesh->PerformNodeMerge(p_node_a,p_node_b);
                                p_node_a->SetAsBoundaryNode(true);

                                recheck_edges = true;
                                break;
                            }
                            else if(node_b_elem_indices.size() >= 2 && node_a_elem_indices.size() == 1)
                            {
                                // Delete neighbour a
                                // std::set<unsigned> element_index_set_a = p_node_a->rGetContainingElementIndices();
                                // unsigned elem_index_1 = (*element_index_set_a.begin());
                                // VertexElement<DIM,DIM>* p_element_1 = p_mesh->GetElement(elem_index_1);
                                // p_element_1->DeleteNode(p_element_1->GetNodeLocalIndex(node_global_index));
                                // p_mesh->DeleteNodePriorToReMesh(node_global_index);

                                p_mesh->PerformNodeMerge(p_node_b,p_node_a);
                                p_node_b->SetAsBoundaryNode(true);

                                recheck_edges = true;
                                break;

                            }
                        }   
                    }
                }
            }
            if(recheck_edges)
            {
                break;
            }
        }
        p_mesh->ReMesh();
    // }
    // unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    // std::stringstream time;
    // time << num_timesteps;
    // std::stringstream t_numb_times_here;
    // t_numb_times_here << numb_times_here;
    // VertexMeshWriter<DIM,DIM> vertexmesh_writer1("tmp", "Mesh_2", false);
    // vertexmesh_writer1.WriteVtkUsingMesh(*p_mesh, time.str() + "_" + t_numb_times_here.str());
    return recheck_edges;
    
}



template<unsigned DIM>
bool JaggedVertexEdgesModifierSimplified<DIM>::ClearFreeInternalNodes(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here)
{
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));

    bool ReCheck_Mesh = false;
    bool performed_edge_modifier = false;
    for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter_1 = p_mesh->GetNodeIteratorBegin();
        node_iter_1 != p_mesh->GetNodeIteratorEnd();
        ++node_iter_1)
    {
        std::set<unsigned> containing_element_indices = node_iter_1->rGetContainingElementIndices();

        if(!(node_iter_1->IsBoundaryNode()) && containing_element_indices.size()==1)
        {
            double min_dist = 10.0;
            unsigned node_to_merge;
            unsigned node_index = node_iter_1->GetIndex();
            std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);
            
            Node<DIM>* p_node = p_mesh->GetNode(node_index);

            for (std::set<unsigned>::iterator neighbour_iter = node_neighbours.begin();
                neighbour_iter != node_neighbours.end();
                ++neighbour_iter)
            {
                unsigned neighbour_index = *neighbour_iter;
                Node<DIM>* p_neighbour_node = p_mesh->GetNode(neighbour_index);

                c_vector<double, DIM> r_node = p_node->rGetLocation();
                c_vector<double, DIM> r_neighbour = p_neighbour_node->rGetLocation();
                
                if(norm_2(p_mesh->GetVectorFromAtoB(r_node,r_neighbour)) < min_dist) 
                {
                    min_dist = norm_2(p_mesh->GetVectorFromAtoB(r_node,r_neighbour));
                    node_to_merge = neighbour_index;
                }

            }
            if(min_dist<10)
            {
                Node<DIM>* p_merge_node = p_mesh->GetNode(node_to_merge);
                p_mesh->PerformNodeMerge(p_merge_node,p_node);
                performed_edge_modifier = true;
                break;
            }
            
        }

    }
    p_mesh->ReMesh();
    // unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    // std::stringstream time;
    // time << num_timesteps;
    // std::stringstream t_numb_times_here;
    // t_numb_times_here << numb_times_here;
    // VertexMeshWriter<DIM,DIM> vertexmesh_writer1("tmp", "Mesh_1", false);
    // vertexmesh_writer1.WriteVtkUsingMesh(*p_mesh, time.str() + "_" + t_numb_times_here.str());

    return ReCheck_Mesh;
}

template<unsigned DIM>
bool JaggedVertexEdgesModifierSimplified<DIM>::DirtyHack(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here)
{
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));

    bool ReCheck_Mesh = false;
    bool performed_edge_modifier = false;
    double distanceBetweenVerteciesThreshold = 1.25*p_mesh->GetCellRearrangementThreshold();
    for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter_1 = p_mesh->GetNodeIteratorBegin();
        node_iter_1 != p_mesh->GetNodeIteratorEnd();
        ++node_iter_1)
    {
        if(!(node_iter_1->IsBoundaryNode()))
        {
            for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter_2 = p_mesh->GetNodeIteratorBegin();
                node_iter_2 != p_mesh->GetNodeIteratorEnd();
                ++node_iter_2)
            {
                std::set<unsigned> containing_element_indices = node_iter_2->rGetContainingElementIndices();

                if(node_iter_2->IsBoundaryNode() && containing_element_indices.size()==1)
                {
                    Node<DIM>* p_node_1 = p_mesh->GetNode(node_iter_1->GetIndex());
                    Node<DIM>* p_node_2 = p_mesh->GetNode(node_iter_2->GetIndex());
                    if(norm_2(p_mesh->GetVectorFromAtoB(p_node_1->rGetLocation(),p_node_2->rGetLocation() )) < distanceBetweenVerteciesThreshold)
                    {
                        
                        p_mesh->PerformNodeMerge(p_node_1,p_node_2);
                        performed_edge_modifier = true;
                        break;
                    }
                }
            }
        }
        

    }
    p_mesh->ReMesh();
    return ReCheck_Mesh;
}


template<unsigned DIM>
bool JaggedVertexEdgesModifierSimplified<DIM>::CheckForVConfig(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here)
{
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));

    bool ReCheck_Mesh = false;
    bool performed_edge_modifier = false;
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
        double distanceBetweenVerteciesThreshold = 5.0*p_mesh->GetCellRearrangementThreshold();
        // double distanceBetweenVerteciesThreshold = 0.02; //0.075
        double distanceToCommonVertexThreshold = 0.5; // must be greater than mMaxEdgeLength in VertexBoundaryRefinementModifier.cpp

        
            for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
                node_iter != p_mesh->GetNodeIteratorEnd();
                ++node_iter)
            {
                if (node_iter->IsBoundaryNode())
                {
                    std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
                    
                    c_vector<double, DIM> r_node = node_iter->rGetLocation();

                    if(containing_element_indices.size() == 2)
                    {
                        
                        unsigned node_index = node_iter->GetIndex();
                        std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);

                        // Get the boundary nodes neighbouring boundary nodes
                        std::vector<unsigned> free_boundary_neighbours;
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
                                    free_boundary_neighbours.push_back(neighbour_index);
                                }

                            }
                            if(p_neighbour_node->IsBoundaryNode())
                            {
                                boundary_neighbours.push_back(neighbour_index);
                            }

                        }

                        // If we have two nodes which are closer to the common vertex than thethreshold distanceToCommonVertexThreshold,
                        // we might be able to delete them if they are close enough to eachother.
                        // if(free_boundary_neighbours.size() == 2 && (boundary_neighbours.size() < node_neighbours.size()) )
                        // {
                        //     // TRACE("Getting elements!");
                        //     Node<DIM>* p_neighbour_1 = p_mesh->GetNode(free_boundary_neighbours[0]);
                        //     Node<DIM>* p_neighbour_2 = p_mesh->GetNode(free_boundary_neighbours[1]);
                        //     Node<DIM>* p_node = p_mesh->GetNode(node_index);

                        //     c_vector<double, DIM> r_neighbour_1 = p_neighbour_1->rGetLocation();
                        //     c_vector<double, DIM> r_neighbour_2 = p_neighbour_2->rGetLocation();

                        
                        //     if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_2)) < distanceBetweenVerteciesThreshold )  
                        //     {


                                
                        //         // TRACE("Performed 3.1 cell edge stitch");
                        //         // PRINT_VECTOR(r_neighbour_1);
                        //         // PRINT_VECTOR(r_neighbour_2);
                        //         // Delete neighbour 1
                        //         std::set<unsigned> element_index_set_1 = p_neighbour_1->rGetContainingElementIndices();
                        //         unsigned elem_index_1 = (*element_index_set_1.begin());
                        //         VertexElement<DIM,DIM>* p_element_1 = p_mesh->GetElement(elem_index_1);
                        //         p_element_1->DeleteNode(p_element_1->GetNodeLocalIndex(free_boundary_neighbours[0]));
                        //         p_mesh->DeleteNodePriorToReMesh(free_boundary_neighbours[0]);

                        //         // Delete neighbour 2
                        //         std::set<unsigned> element_index_set_2 = p_neighbour_2->rGetContainingElementIndices();
                        //         unsigned elem_index_2 = (*element_index_set_2.begin());
                        //         VertexElement<DIM,DIM>* p_element_2 = p_mesh->GetElement(elem_index_2);
                        //         p_element_2->DeleteNode(p_element_2->GetNodeLocalIndex(free_boundary_neighbours[1]));
                        //         p_mesh->DeleteNodePriorToReMesh(free_boundary_neighbours[1]);

                        //         // Move node to where the other 2 used to be:
                        //         p_node->rGetModifiableLocation() = r_neighbour_1 + 0.5 * p_mesh->GetVectorFromAtoB(r_neighbour_1, r_neighbour_2);
                        //         p_node->SetAsBoundaryNode(true);

                        //         performed_edge_modifier = true;
                        //         // TRACE("Performed 3.2 cell edge stitch");
                        //         break;
                                    

                        //     }
                            

                        // }
                        // else if(free_boundary_neighbours.size() == 2 && (boundary_neighbours.size() == node_neighbours.size()) )
                        if(free_boundary_neighbours.size() == 2 )
                        {
                            // TRACE("Performed 3.2 cell edge stitch");
                            // PRINT_VARIABLE(numb_times_here);

                            // unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
                            // std::stringstream time;
                            // time << num_timesteps;
                            // VertexMeshWriter<DIM,DIM> vertexmesh_writer("tmp", "32_mesh", false);
                            // std::stringstream out;
                            // out << numb_times_here;
                            // std::string write = time.str() + "_" + out.str();
                            // vertexmesh_writer.WriteVtkUsingMesh(*p_mesh, write);

                            // PRINT_3_VARIABLES(free_boundary_neighbours.size(),boundary_neighbours.size(),node_neighbours.size());
                            
                            Node<DIM>* p_neighbour_1 = p_mesh->GetNode(free_boundary_neighbours[0]);
                            Node<DIM>* p_neighbour_2 = p_mesh->GetNode(free_boundary_neighbours[1]);
                            Node<DIM>* p_node = p_mesh->GetNode(node_index);

                            c_vector<double, DIM> r_neighbour_1 = p_neighbour_1->rGetLocation();
                            c_vector<double, DIM> r_neighbour_2 = p_neighbour_2->rGetLocation();

                            // PRINT_VECTOR(r_neighbour_1);
                            // PRINT_VECTOR(r_neighbour_2);
                        
                            if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_2)) < distanceBetweenVerteciesThreshold )  
                            {
                                // TRACE("doing the thing 3.2");
                                p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);
                                p_neighbour_1->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                break;

                            }


                        }
                        else if(free_boundary_neighbours.size() == 4 && (boundary_neighbours.size() == node_neighbours.size()) )
                        {
                            Node<DIM>* p_neighbour_1 = p_mesh->GetNode(free_boundary_neighbours[0]);
                            Node<DIM>* p_neighbour_2 = p_mesh->GetNode(free_boundary_neighbours[1]);
                            Node<DIM>* p_neighbour_3 = p_mesh->GetNode(free_boundary_neighbours[2]);
                            Node<DIM>* p_neighbour_4 = p_mesh->GetNode(free_boundary_neighbours[3]);
                            c_vector<double, DIM> r_neighbour_1 = p_neighbour_1->rGetLocation();
                            c_vector<double, DIM> r_neighbour_2 = p_neighbour_2->rGetLocation();
                            c_vector<double, DIM> r_neighbour_3 = p_neighbour_3->rGetLocation();
                            c_vector<double, DIM> r_neighbour_4 = p_neighbour_4->rGetLocation();

                            if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_2)) < distanceBetweenVerteciesThreshold )  
                            {
                                // TRACE("doing the thing 3.2");
                                p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);
                                p_neighbour_1->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                break;
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_3)) < distanceBetweenVerteciesThreshold )  
                            {
                                // TRACE("doing the thing 3.2");
                                p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_3);
                                p_neighbour_1->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                break;
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_1,r_neighbour_4)) < distanceBetweenVerteciesThreshold )  
                            {
                                // TRACE("doing the thing 3.2");
                                p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_4);
                                p_neighbour_1->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                break;
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_2,r_neighbour_3)) < distanceBetweenVerteciesThreshold )  
                            {
                                // TRACE("doing the thing 3.2");
                                p_mesh->PerformNodeMerge(p_neighbour_2,p_neighbour_3);
                                p_neighbour_1->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                break;
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_2,r_neighbour_4)) < distanceBetweenVerteciesThreshold )  
                            {
                                // TRACE("doing the thing 3.2");
                                p_mesh->PerformNodeMerge(p_neighbour_2,p_neighbour_4);
                                p_neighbour_1->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                break;
                            }
                            else if(norm_2(p_mesh->GetVectorFromAtoB(r_neighbour_3,r_neighbour_4)) < distanceBetweenVerteciesThreshold )  
                            {
                                // TRACE("doing the thing 3.2");
                                p_mesh->PerformNodeMerge(p_neighbour_3,p_neighbour_4);
                                p_neighbour_1->SetAsBoundaryNode(true);
                                performed_edge_modifier = true;
                                break;
                            }

                        }

                    }
                    
                }
                if(performed_edge_modifier)
                {
                    break;
                }
            }
        

        if(performed_edge_modifier)
        {
            ReCheck_Mesh = true;
            p_mesh->ReMesh();

            // unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
            // std::stringstream time2;
            // time2 << num_timesteps;
            // VertexMeshWriter<DIM,DIM> vertexmesh_writer2("tmp", "post_32_mesh", false);
            // std::stringstream out2;
            // out2 << numb_times_here;
            // std::string write2 = time2.str() + "_" + out2.str();
            // vertexmesh_writer2.WriteVtkUsingMesh(*p_mesh, write2);
        }
    // unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    // std::stringstream time;
    // time << num_timesteps;
    // std::stringstream t_numb_times_here;
    // t_numb_times_here << numb_times_here;
    // VertexMeshWriter<DIM,DIM> vertexmesh_writer1("tmp", "Mesh_3", false);
    // vertexmesh_writer1.WriteVtkUsingMesh(*p_mesh, time.str() + "_" + t_numb_times_here.str());
    return ReCheck_Mesh;

}

template<unsigned DIM>
bool JaggedVertexEdgesModifierSimplified<DIM>::CheckFor3NodeVoids(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here)
{
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));

    bool ReCheck_Mesh = false;

    bool performed_edge_modifier = false;
    
    unsigned node_to_print = 0;

    double mVoidAreaThreshold = 0.05;
    for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
        node_iter != p_mesh->GetNodeIteratorEnd();
        ++node_iter)
    {
        if (node_iter->IsBoundaryNode())
        {
            std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
            if(containing_element_indices.size() >= 2)
            {
                unsigned node_index = node_iter->GetIndex();
                std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);
                
                // PRINT_VARIABLE(node_neighbours.size());

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
                    if(p_neighbour_node_1->IsBoundaryNode() && p_neighbour_node_2->IsBoundaryNode())
                    {
                        std::set<unsigned> containing_element_indices_1 = p_neighbour_node_1->rGetContainingElementIndices();
                        std::set<unsigned> containing_element_indices_2 = p_neighbour_node_2->rGetContainingElementIndices();

                        // Find common elements
                        std::set<unsigned> shared_elements;
                        std::set_intersection(containing_element_indices_1.begin(),
                                            containing_element_indices_1.end(),
                                            containing_element_indices_2.begin(),
                                            containing_element_indices_2.end(),
                                            std::inserter(shared_elements, shared_elements.begin()));
                        // i.e. they are neighbours
                        if(shared_elements.size()==1 && containing_element_indices_1.size()==2 && containing_element_indices_2.size()==2)
                        {
                            performed_edge_modifier = Remove3NodeVoid(performed_edge_modifier, node_index, p_neighbour_node_1, p_neighbour_node_2, p_cell_population, p_mesh);
                            
                            if(performed_edge_modifier)
                            {
                                if (p_neighbour_node_1->IsDeleted() && p_neighbour_node_2->IsDeleted())
                                {
                                    node_iter->SetAsBoundaryNode(true);
                                }
                                else if (node_iter->IsDeleted() && p_neighbour_node_2->IsDeleted())
                                {
                                    p_neighbour_node_1->SetAsBoundaryNode(true);
                                }
                                else if (node_iter->IsDeleted() && p_neighbour_node_1->IsDeleted())
                                {
                                    p_neighbour_node_2->SetAsBoundaryNode(true);
                                }
                                break;
                            }
                        }
                                                    
                              
                    }
                    else if(p_neighbour_node_3->IsBoundaryNode() && p_neighbour_node_2->IsBoundaryNode() )
                    {
                        std::set<unsigned> containing_element_indices_1 = p_neighbour_node_3->rGetContainingElementIndices();
                        std::set<unsigned> containing_element_indices_2 = p_neighbour_node_2->rGetContainingElementIndices();
                        // Find common elements
                        std::set<unsigned> shared_elements;
                        std::set_intersection(containing_element_indices_1.begin(),
                                            containing_element_indices_1.end(),
                                            containing_element_indices_2.begin(),
                                            containing_element_indices_2.end(),
                                            std::inserter(shared_elements, shared_elements.begin()));
                        // i.e. they are neighbours
                        if(shared_elements.size()==1 && containing_element_indices_1.size()==2 && containing_element_indices_2.size()==2)
                        {
                            performed_edge_modifier = Remove3NodeVoid(performed_edge_modifier, node_index, p_neighbour_node_3, p_neighbour_node_2, p_cell_population, p_mesh);

                            if(performed_edge_modifier)
                            {
                                // TRACE("Removed Node!!    2");
                                if (p_neighbour_node_3->IsDeleted() && p_neighbour_node_2->IsDeleted())
                                {
                                    node_iter->SetAsBoundaryNode(true);
                                }
                                else if (node_iter->IsDeleted() && p_neighbour_node_2->IsDeleted())
                                {
                                    p_neighbour_node_3->SetAsBoundaryNode(true);
                                }
                                else if (node_iter->IsDeleted() && p_neighbour_node_3->IsDeleted())
                                {
                                    // TRACE("3");
                                    PRINT_VARIABLE(p_neighbour_node_2->IsDeleted());
                                    PRINT_VARIABLE(p_neighbour_node_2->IsBoundaryNode());
                                    PRINT_VECTOR(p_neighbour_node_2->rGetLocation());
                                    PRINT_VARIABLE(node_neighbour_2);
                                    node_to_print = node_neighbour_2;
                                    p_neighbour_node_2->SetAsBoundaryNode(true);
                                }
                                break;
                            }
                        }
                    }
                    else if(p_neighbour_node_1->IsBoundaryNode() && p_neighbour_node_3->IsBoundaryNode())
                    {
                        std::set<unsigned> containing_element_indices_1 = p_neighbour_node_1->rGetContainingElementIndices();
                        std::set<unsigned> containing_element_indices_2 = p_neighbour_node_3->rGetContainingElementIndices();
                        // Find common elements
                        std::set<unsigned> shared_elements;
                        std::set_intersection(containing_element_indices_1.begin(),
                                            containing_element_indices_1.end(),
                                            containing_element_indices_2.begin(),
                                            containing_element_indices_2.end(),
                                            std::inserter(shared_elements, shared_elements.begin()));
                        // i.e. they are neighbours
                        if(shared_elements.size()==1 && containing_element_indices_1.size()==2 && containing_element_indices_2.size()==2)
                        {
                            performed_edge_modifier = Remove3NodeVoid(performed_edge_modifier, node_index, p_neighbour_node_1, p_neighbour_node_3, p_cell_population, p_mesh);

                            if(performed_edge_modifier)
                            {
                                if (p_neighbour_node_1->IsDeleted() && p_neighbour_node_3->IsDeleted())
                                {
                                    node_iter->SetAsBoundaryNode(true);
                                }
                                else if (node_iter->IsDeleted() && p_neighbour_node_3->IsDeleted())
                                {
                                    p_neighbour_node_1->SetAsBoundaryNode(true);
                                }
                                else if (node_iter->IsDeleted() && p_neighbour_node_1->IsDeleted())
                                {
                                    p_neighbour_node_3->SetAsBoundaryNode(true);
                                }
                                break;
                            }
                        }
                    }

                }
                else if(node_neighbours.size() == 4)
                {
                    std::set<unsigned>::const_iterator neigh_it = node_neighbours.begin();
                    unsigned node_neighbour_1 = (*neigh_it);
                    neigh_it++;
                    unsigned node_neighbour_2 = (*neigh_it);
                    neigh_it++;
                    unsigned node_neighbour_3 = (*neigh_it);
                    neigh_it++;
                    unsigned node_neighbour_4 = (*neigh_it);

                    Node<DIM>* p_neighbour_node_1 = p_mesh->GetNode(node_neighbour_1);
                    Node<DIM>* p_neighbour_node_2 = p_mesh->GetNode(node_neighbour_2);
                    Node<DIM>* p_neighbour_node_3 = p_mesh->GetNode(node_neighbour_3);
                    Node<DIM>* p_neighbour_node_4 = p_mesh->GetNode(node_neighbour_4);

                    std::set<unsigned> containing_element_indices_1 = p_neighbour_node_1->rGetContainingElementIndices();
                    std::set<unsigned> containing_element_indices_2 = p_neighbour_node_2->rGetContainingElementIndices();
                    std::set<unsigned> containing_element_indices_3 = p_neighbour_node_3->rGetContainingElementIndices();
                    std::set<unsigned> containing_element_indices_4 = p_neighbour_node_4->rGetContainingElementIndices();

                    // Find common elements
                    std::set<unsigned> shared_elements_12;
                    std::set_intersection(containing_element_indices_1.begin(),containing_element_indices_1.end(),
                        containing_element_indices_2.begin(),containing_element_indices_2.end(),std::inserter(shared_elements_12, shared_elements_12.begin()));

                    std::set<unsigned> shared_elements_13;
                    std::set_intersection(containing_element_indices_1.begin(),containing_element_indices_1.end(),
                        containing_element_indices_3.begin(),containing_element_indices_3.end(),std::inserter(shared_elements_13, shared_elements_13.begin()));

                    std::set<unsigned> shared_elements_14;
                    std::set_intersection(containing_element_indices_1.begin(),containing_element_indices_1.end(),
                        containing_element_indices_4.begin(),containing_element_indices_4.end(),std::inserter(shared_elements_14, shared_elements_14.begin()));
                    
                    std::set<unsigned> shared_elements_23;
                    std::set_intersection(containing_element_indices_2.begin(),containing_element_indices_2.end(),
                        containing_element_indices_3.begin(),containing_element_indices_3.end(),std::inserter(shared_elements_23, shared_elements_23.begin()));
                    
                    std::set<unsigned> shared_elements_24;
                    std::set_intersection(containing_element_indices_2.begin(),containing_element_indices_2.end(),
                        containing_element_indices_4.begin(),containing_element_indices_4.end(),std::inserter(shared_elements_24, shared_elements_24.begin()));

                    std::set<unsigned> shared_elements_34;
                    std::set_intersection(containing_element_indices_3.begin(),containing_element_indices_3.end(),
                        containing_element_indices_4.begin(),containing_element_indices_4.end(),std::inserter(shared_elements_34, shared_elements_34.begin()));

                    if(p_neighbour_node_1->IsBoundaryNode() && p_neighbour_node_2->IsBoundaryNode() && shared_elements_12.size()==1)
                    {
                        // i.e. they are neighbours
                        if(containing_element_indices_1.size()==2 && containing_element_indices_2.size()==2)
                        {
                            performed_edge_modifier = Remove3NodeVoid(performed_edge_modifier, node_index, p_neighbour_node_1, p_neighbour_node_2, p_cell_population, p_mesh);
                            break;
                        }
                    }
                    else if(p_neighbour_node_1->IsBoundaryNode() && p_neighbour_node_3->IsBoundaryNode() && shared_elements_13.size()==1)
                    {
                        // i.e. they are neighbours
                        if(containing_element_indices_1.size()==2 && containing_element_indices_3.size()==2)
                        {
                            performed_edge_modifier = Remove3NodeVoid(performed_edge_modifier, node_index, p_neighbour_node_1, p_neighbour_node_3, p_cell_population, p_mesh);
                            break;
                        }
                    }
                    else if(p_neighbour_node_1->IsBoundaryNode() && p_neighbour_node_4->IsBoundaryNode() && shared_elements_14.size()==1)
                    {
                        // i.e. they are neighbours
                        if(containing_element_indices_1.size()==2 && containing_element_indices_4.size()==2)
                        {
                            performed_edge_modifier = Remove3NodeVoid(performed_edge_modifier, node_index, p_neighbour_node_1, p_neighbour_node_4, p_cell_population, p_mesh);
                            break;
                        }
                    }
                    else if(p_neighbour_node_2->IsBoundaryNode() && p_neighbour_node_3->IsBoundaryNode() && shared_elements_23.size()==1)
                    {
                        // i.e. they are neighbours
                        if(containing_element_indices_2.size()==2 && containing_element_indices_3.size()==2)
                        {
                            performed_edge_modifier = Remove3NodeVoid(performed_edge_modifier, node_index, p_neighbour_node_2, p_neighbour_node_3, p_cell_population, p_mesh);
                            break;
                        }
                    }
                    else if(p_neighbour_node_2->IsBoundaryNode() && p_neighbour_node_4->IsBoundaryNode() && shared_elements_24.size()==1)
                    {
                        // i.e. they are neighbours
                        if(containing_element_indices_2.size()==2 && containing_element_indices_4.size()==2)
                        {
                            performed_edge_modifier = Remove3NodeVoid(performed_edge_modifier, node_index, p_neighbour_node_2, p_neighbour_node_4, p_cell_population, p_mesh);
                            break;
                        }
                    }
                    else if(p_neighbour_node_3->IsBoundaryNode() && p_neighbour_node_4->IsBoundaryNode() && shared_elements_34.size()==1)
                    {
                        // i.e. they are neighbours
                        if(containing_element_indices_3.size()==2 && containing_element_indices_4.size()==2)
                        {
                            performed_edge_modifier = Remove3NodeVoid(performed_edge_modifier, node_index, p_neighbour_node_3, p_neighbour_node_4, p_cell_population, p_mesh);
                            break;
                        }
                    }

                }
            }
                    /*
                    *
                    *     |           
                    *   --o------o
                    *       \    |   
                    *        \   |  
                    *         \  |  
                    *          \ |
                    *            o
                    *            |
                    */          
                    // else if(containing_element_indices.size() == 1)
                    // {
                    //     unsigned node_index = node_iter->GetIndex();
                    //     std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);
                        
                    //     if(node_neighbours.size() == 2)
                    //     {
                    //         std::set<unsigned>::const_iterator neigh_it = node_neighbours.begin();
                    //         unsigned node_neighbour_1 = (*neigh_it);
                    //         neigh_it++;
                    //         unsigned node_neighbour_2 = (*neigh_it);

                    //         Node<DIM>* p_neighbour_node_1 = p_mesh->GetNode(node_neighbour_1);
                    //         Node<DIM>* p_neighbour_node_2 = p_mesh->GetNode(node_neighbour_2);

                    //         if(p_neighbour_node_1->IsBoundaryNode() && p_neighbour_node_2->IsBoundaryNode() && IsNodeABNeighbours(node_neighbour_1,node_neighbour_2, rCellPopulation ))
                    //         {
                    //             // Delete neighbour 1
                    //             std::set<unsigned> element_index_set_1 = node_iter->rGetContainingElementIndices();
                    //             unsigned elem_index_1 = (*element_index_set_1.begin());
                    //             VertexElement<DIM,DIM>* p_element_1 = p_mesh->GetElement(elem_index_1);
                    //             p_element_1->DeleteNode(p_element_1->GetNodeLocalIndex(node_index));
                    //             p_mesh->DeleteNodePriorToReMesh(node_index);

                    //             // Check if the boundary nodes need relabeling
                    //             std::set<unsigned> node_neighbours_1 = p_cell_population->GetNeighbouringNodeIndices(node_neighbour_1);
                    //             bool is_neighbour_node_1_boundary = false;
                    //             unsigned numb_boundary_neighs_1 = 0;
                    //             for (std::set<unsigned>::iterator neighbour_iter = node_neighbours_1.begin();
                    //                                     neighbour_iter != node_neighbours_1.end();
                    //                                     ++neighbour_iter)
                    //             {
                    //                 unsigned node_neigh_index = *neighbour_iter;
                    //                 Node<DIM>* p_node_iter = p_mesh->GetNode(node_neigh_index);
                    //                 std::set<unsigned> containing_element_indices_neigh = p_node_iter->rGetContainingElementIndices();

                    //                 if(p_node_iter->IsBoundaryNode() && node_neigh_index != node_index && containing_element_indices_neigh.size()==1)
                    //                 {
                    //                     // PRINT_VECTOR(p_node_iter->rGetLocation());
                    //                     is_neighbour_node_1_boundary = true;
                    //                 }
                    //                 if(p_node_iter->IsBoundaryNode())
                    //                 {
                    //                     numb_boundary_neighs_1++;
                    //                 }
                    //             }
                    //             if(is_neighbour_node_1_boundary==false && numb_boundary_neighs_1==(node_neighbours_1.size()))
                    //             {
                    //                 is_neighbour_node_1_boundary = true;
                    //             }
                    //             // p_neighbour_node_1->SetAsBoundaryNode(is_neighbour_node_1_boundary);
                    //             // PRINT_VECTOR(p_neighbour_node_1->rGetLocation());
                    //             // PRINT_VARIABLE(is_neighbour_node_1_boundary);

                    //             std::set<unsigned> node_neighbours_2 = p_cell_population->GetNeighbouringNodeIndices(node_neighbour_2);
                    //             bool is_neighbour_node_2_boundary = false;
                    //             unsigned numb_boundary_neighs_2 = 0;
                    //             for (std::set<unsigned>::iterator neighbour_iter = node_neighbours_2.begin();
                    //                                     neighbour_iter != node_neighbours_2.end();
                    //                                     ++neighbour_iter)
                    //             {
                    //                 unsigned node_neigh_index = *neighbour_iter;
                    //                 Node<DIM>* p_node_iter = p_mesh->GetNode(node_neigh_index);
                    //                 std::set<unsigned> containing_element_indices_neigh = p_node_iter->rGetContainingElementIndices();

                    //                 if(p_node_iter->IsBoundaryNode() && node_neigh_index != node_index && containing_element_indices_neigh.size()==1)
                    //                 {
                    //                     // PRINT_VECTOR(p_node_iter->rGetLocation());
                    //                     is_neighbour_node_2_boundary = true;
                    //                 }
                    //                 if(p_node_iter->IsBoundaryNode())
                    //                 {
                    //                     numb_boundary_neighs_2++;
                    //                 }
                    //             }
                    //             if(is_neighbour_node_2_boundary==false && numb_boundary_neighs_2==(node_neighbours_1.size()))
                    //             {
                    //                 is_neighbour_node_2_boundary = true;
                    //             }

                    //             p_neighbour_node_1->SetAsBoundaryNode(is_neighbour_node_1_boundary);
                    //             // PRINT_VECTOR(p_neighbour_node_1->rGetLocation());
                    //             // PRINT_VARIABLE(is_neighbour_node_1_boundary);

                    //             p_neighbour_node_2->SetAsBoundaryNode(is_neighbour_node_2_boundary);
                    //             // PRINT_VECTOR(p_neighbour_node_2->rGetLocation());
                    //             // PRINT_VARIABLE(is_neighbour_node_2_boundary);

                                
                    //             performed_edge_modifier = true;
                    //             // TRACE("Performed 1 node Void Removel");
                                
                    //         }
                    //     }
                    // }
        }
        if(performed_edge_modifier)
        {
            break;
        }
    }

        

    if(performed_edge_modifier)
    {
        // Node<DIM>* p_neighbour_node_1 = p_mesh->GetNode(node_to_print);
        // PRINT_VARIABLE(p_neighbour_node_1->IsBoundaryNode());
        ReCheck_Mesh = true;
        p_mesh->ReMesh();
        // PRINT_VARIABLE(p_neighbour_node_1->IsBoundaryNode());
        
    }
    // unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    // std::stringstream time;
    // time << num_timesteps;
    // std::stringstream t_numb_times_here;
    // t_numb_times_here << numb_times_here;
    // VertexMeshWriter<DIM,DIM> vertexmesh_writer1("tmp", "Mesh_4", false);
    // vertexmesh_writer1.WriteVtkUsingMesh(*p_mesh, time.str() + "_" + t_numb_times_here.str());

    return ReCheck_Mesh;

}

template<unsigned DIM>
bool JaggedVertexEdgesModifierSimplified<DIM>::CleanUpRogueNodes(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here)
{
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));
    
    bool ReCheck_Mesh = false;


            for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
                node_iter != p_mesh->GetNodeIteratorEnd();
                ++node_iter)
            {
                bool is_boundary = false;
                bool is_def_bound = false;
                // bool is_boundary_tmp = false;
                // if (node_iter->IsBoundaryNode())
                // {
                    std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
                    if(containing_element_indices.size() == 1)
                    {
                        is_boundary = true;
                    }
                    else if(containing_element_indices.size() >= 2)
                    {
                        unsigned node_index = node_iter->GetIndex();
                        std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);
                        if(node_neighbours.size() > containing_element_indices.size())
                        {
                            is_boundary = true;
                        }
                    }

                    // else if(containing_element_indices.size() == 2)
                    // {
                    //     unsigned node_index = node_iter->GetIndex();
                    //     std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);
                    //     unsigned number_of_boundary_neighbours = 0;
                    //     std::vector<unsigned> boundary_neighbours;
                    //     for (std::set<unsigned>::iterator neighbour_iter = node_neighbours.begin();
                    //         neighbour_iter != node_neighbours.end();
                    //         ++neighbour_iter)
                    //     {
                    //         unsigned neighbour_index = *neighbour_iter;
                    //         Node<DIM>* p_neighbour_node = p_mesh->GetNode(neighbour_index);

                    //         if(p_neighbour_node->IsBoundaryNode())
                    //         {
                    //             boundary_neighbours.push_back(neighbour_index);
                    //             number_of_boundary_neighbours++;
                    //             std::set<unsigned> containing_element_indices_neighbour = p_neighbour_node->rGetContainingElementIndices();
                    //             if(containing_element_indices_neighbour.size() == 1)
                    //             {
                    //                 is_boundary = true;
                    //             }

                    //         }
                    //     }
                    //     if(is_boundary == false)
                    //     {
                    //         if(boundary_neighbours.size()==2)
                    //         {
                    //             Node<DIM>* p_neighbour_1 = p_mesh->GetNode(boundary_neighbours[0]);
                    //             Node<DIM>* p_neighbour_2 = p_mesh->GetNode(boundary_neighbours[1]);
                    //             std::set<unsigned> containing_element_indices_0 = node_iter->rGetContainingElementIndices();
                    //             std::set<unsigned> containing_element_indices_1 = p_neighbour_1->rGetContainingElementIndices();
                    //             std::set<unsigned> containing_element_indices_2 = p_neighbour_2->rGetContainingElementIndices();

                    //             // Find common elements
                    //             std::set<unsigned> shared_elements_0_1;
                    //             std::set_intersection(containing_element_indices_0.begin(),
                    //                             containing_element_indices_0.end(),
                    //                             containing_element_indices_1.begin(),
                    //                             containing_element_indices_1.end(),
                    //                             std::inserter(shared_elements_0_1, shared_elements_0_1.begin()));
                    //             std::set<unsigned> shared_elements_0_2;
                    //             std::set_intersection(containing_element_indices_0.begin(),
                    //                             containing_element_indices_0.end(),
                    //                             containing_element_indices_2.begin(),
                    //                             containing_element_indices_2.end(),
                    //                             std::inserter(shared_elements_0_2, shared_elements_0_2.begin()));
                    //             std::set<unsigned> shared_elements_1_2;
                    //             std::set_intersection(containing_element_indices_1.begin(),
                    //                             containing_element_indices_1.end(),
                    //                             containing_element_indices_2.begin(),
                    //                             containing_element_indices_2.end(),
                    //                             std::inserter(shared_elements_1_2, shared_elements_1_2.begin()));
                    //             if(containing_element_indices_0.size()==2 && containing_element_indices_1.size()==2 && containing_element_indices_2.size()==2)
                    //             {
                    //                 if(shared_elements_0_1.size()==2 && shared_elements_0_2.size()==2 && shared_elements_1_2.size()==1)
                    //                 {
                    //                     // case 1
                    //                     is_boundary = false;

                    //                 }
                    //                 else if(shared_elements_0_1.size()==1 && shared_elements_0_2.size()==1 && shared_elements_1_2.size()==1)
                    //                 {
                    //                     // case 2
                    //                     is_boundary = true;
                    //                 }
                    //                 else if(shared_elements_0_1.size()==1 && shared_elements_0_2.size()==1 && shared_elements_1_2.size()==0)
                    //                 {
                    //                     // case 3
                    //                     is_boundary = true;
                    //                 }
                    //             }
                    //             else if(containing_element_indices_0.size()>2 && containing_element_indices_1.size()>2 && containing_element_indices_2.size()>2)
                    //             {
                    //                 if(shared_elements_0_1.size()==2 && shared_elements_0_2.size()==2 && shared_elements_1_2.size()==1)
                    //                 {
                    //                     is_boundary = false;

                    //                 }
                    //                 else if(shared_elements_0_1.size()==1 && shared_elements_0_2.size()==1 && shared_elements_1_2.size()==1)
                    //                 {
                    //                     is_boundary = true;
                    //                 }
                    //             }
                                

                    //         }
                            

                    //     }

                    // }
                    // else if(containing_element_indices.size() == 3)
                    // {
                    //     unsigned node_index = node_iter->GetIndex();
                    //     std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);
                    //     unsigned number_of_boundary_neighbours = 0;
                    //     std::vector<unsigned> boundary_neighbours;
                    //     for (std::set<unsigned>::iterator neighbour_iter = node_neighbours.begin();
                    //         neighbour_iter != node_neighbours.end();
                    //         ++neighbour_iter)
                    //     {
                    //         unsigned neighbour_index = *neighbour_iter;
                    //         Node<DIM>* p_neighbour_node = p_mesh->GetNode(neighbour_index);

                    //         if(p_neighbour_node->IsBoundaryNode())
                    //         {
                    //             boundary_neighbours.push_back(neighbour_index);
                    //             number_of_boundary_neighbours++;
                    //             std::set<unsigned> containing_element_indices_neighbour = p_neighbour_node->rGetContainingElementIndices();
                    //             if(containing_element_indices_neighbour.size() == 1)
                    //             {
                    //                 is_boundary = true;
                    //             }

                    //         }
                    //     }
                    //     if(is_boundary == false)
                    //     {
                    //         if(boundary_neighbours.size()==2)
                    //         {
                    //             Node<DIM>* p_neighbour_1 = p_mesh->GetNode(boundary_neighbours[0]);
                    //             Node<DIM>* p_neighbour_2 = p_mesh->GetNode(boundary_neighbours[1]);
                    //             std::set<unsigned> containing_element_indices_0 = node_iter->rGetContainingElementIndices();
                    //             std::set<unsigned> containing_element_indices_1 = p_neighbour_1->rGetContainingElementIndices();
                    //             std::set<unsigned> containing_element_indices_2 = p_neighbour_2->rGetContainingElementIndices();

                    //             // Find common elements
                    //             std::set<unsigned> shared_elements_0_1;
                    //             std::set_intersection(containing_element_indices_0.begin(),
                    //                             containing_element_indices_0.end(),
                    //                             containing_element_indices_1.begin(),
                    //                             containing_element_indices_1.end(),
                    //                             std::inserter(shared_elements_0_1, shared_elements_0_1.begin()));
                    //             std::set<unsigned> shared_elements_0_2;
                    //             std::set_intersection(containing_element_indices_0.begin(),
                    //                             containing_element_indices_0.end(),
                    //                             containing_element_indices_2.begin(),
                    //                             containing_element_indices_2.end(),
                    //                             std::inserter(shared_elements_0_2, shared_elements_0_2.begin()));
                    //             std::set<unsigned> shared_elements_1_2;
                    //             std::set_intersection(containing_element_indices_1.begin(),
                    //                             containing_element_indices_1.end(),
                    //                             containing_element_indices_2.begin(),
                    //                             containing_element_indices_2.end(),
                    //                             std::inserter(shared_elements_1_2, shared_elements_1_2.begin()));
                    //             if(containing_element_indices_0.size()>=2 && containing_element_indices_1.size()>=2 && containing_element_indices_2.size()>=2)
                    //             {
                    //                 if(shared_elements_0_1.size()==2 && shared_elements_0_2.size()==2 && shared_elements_1_2.size()==1)
                    //                 {
                    //                     is_boundary = false;
                    //                 }
                    //                 else if(shared_elements_0_1.size()==1 && shared_elements_0_2.size()==1 && shared_elements_1_2.size()==1)
                    //                 {
                    //                     is_boundary = true;
                    //                 }
                    //             }


                    //         }
                    //     }
                        
                    // }                    
                // }

                node_iter->SetAsBoundaryNode(is_boundary);
            }

            // for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
            //     node_iter != p_mesh->GetNodeIteratorEnd();
            //     ++node_iter)
            // {
            //     bool is_boundary = false;
            //     bool is_def_bound = false;
            //     // bool is_boundary_tmp = false;
            //     if (node_iter->IsBoundaryNode())
            //     {
            //         std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
            //         if(containing_element_indices.size() == 1)
            //         {
            //             is_boundary = true;
            //         }
            //         else if(containing_element_indices.size() == 2)
            //         {
            //             unsigned node_index = node_iter->GetIndex();
            //             std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);
            //             unsigned number_of_boundary_neighbours = 0;
            //             std::vector<unsigned> boundary_neighbours;
            //             for (std::set<unsigned>::iterator neighbour_iter = node_neighbours.begin();
            //                 neighbour_iter != node_neighbours.end();
            //                 ++neighbour_iter)
            //             {
            //                 unsigned neighbour_index = *neighbour_iter;
            //                 Node<DIM>* p_neighbour_node = p_mesh->GetNode(neighbour_index);

            //                 if(p_neighbour_node->IsBoundaryNode())
            //                 {
            //                     boundary_neighbours.push_back(neighbour_index);
            //                     number_of_boundary_neighbours++;
            //                     std::set<unsigned> containing_element_indices_neighbour = p_neighbour_node->rGetContainingElementIndices();
            //                     if(containing_element_indices_neighbour.size() == 1)
            //                     {
            //                         is_def_bound = true;
            //                     }

            //                 }
            //             }

            //             if(is_def_bound==true)
            //             {
            //                 is_boundary = true;
            //             }
            //             else
            //             {
            //                 if(number_of_boundary_neighbours == 2)
            //                 {
            //                     Node<DIM>* p_neighbour_1 = p_mesh->GetNode(boundary_neighbours[0]);
            //                     Node<DIM>* p_neighbour_2 = p_mesh->GetNode(boundary_neighbours[1]);

            //                     std::set<unsigned> containing_element_indices_1 = p_neighbour_1->rGetContainingElementIndices();
            //                     std::set<unsigned> containing_element_indices_2 = p_neighbour_2->rGetContainingElementIndices();

            //                     // Find common elements
            //                     std::set<unsigned> shared_elements;
            //                     std::set_intersection(containing_element_indices_1.begin(),
            //                                     containing_element_indices_1.end(),
            //                                     containing_element_indices_2.begin(),
            //                                     containing_element_indices_2.end(),
            //                                     std::inserter(shared_elements, shared_elements.begin()));

            //                     if(shared_elements.size()==1)
            //                     {
            //                         is_boundary = true;
            //                         unsigned neigh_1 = boundary_neighbours[0];
            //                         unsigned neigh_2 = boundary_neighbours[1];

            //                         std::set<unsigned> node_neighbours_1 = p_cell_population->GetNeighbouringNodeIndices(neigh_1);
            //                         std::set<unsigned> node_neighbours_2 = p_cell_population->GetNeighbouringNodeIndices(neigh_2);
                                    
            //                         // if neighbours are also eachothers neighbours, we have a void!
            //                         if((std::count(node_neighbours_1.begin(), node_neighbours_1.end(), neigh_2)) && (std::count(node_neighbours_2.begin(), node_neighbours_2.end(), neigh_1)) )
            //                         {
            //                             is_def_bound = true;
            //                         }
            //                         else
            //                         {
            //                             is_boundary = true;
            //                         }

            //                     }
            //                     else if(containing_element_indices_1.size()>=2 || containing_element_indices_2.size()>=2)
            //                     {
            //                         is_boundary = true;
            //                     }
            //                 }
            //                 if(is_def_bound==true)
            //                 {
            //                     is_boundary = true;
            //                 }
                            
            //             }
            //         }
            //         else if(containing_element_indices.size() == 3)
            //         {
            //             unsigned node_index = node_iter->GetIndex();
            //             std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);
            //             unsigned number_of_boundary_neighbours = 0;
            //             std::vector<unsigned> boundary_neighbours;
            //             for (std::set<unsigned>::iterator neighbour_iter = node_neighbours.begin();
            //                 neighbour_iter != node_neighbours.end();
            //                 ++neighbour_iter)
            //             {
            //                 unsigned neighbour_index = *neighbour_iter;
            //                 Node<DIM>* p_neighbour_node = p_mesh->GetNode(neighbour_index);

            //                 if(p_neighbour_node->IsBoundaryNode())
            //                 {
            //                     boundary_neighbours.push_back(neighbour_index);
            //                     number_of_boundary_neighbours++;
            //                     std::set<unsigned> containing_element_indices_neighbour = p_neighbour_node->rGetContainingElementIndices();
            //                     if(containing_element_indices_neighbour.size() == 1)
            //                     {
            //                         is_def_bound = true;
            //                     }

            //                 }
            //             }
            //             if(is_def_bound==true)
            //             {
            //                 is_boundary = true;
            //             }
            //             else
            //             {
            //                 if(number_of_boundary_neighbours == 2)
            //                 {
            //                     Node<DIM>* p_neighbour_1 = p_mesh->GetNode(boundary_neighbours[0]);
            //                     Node<DIM>* p_neighbour_2 = p_mesh->GetNode(boundary_neighbours[1]);

            //                     std::set<unsigned> containing_element_indices_1 = p_neighbour_1->rGetContainingElementIndices();
            //                     std::set<unsigned> containing_element_indices_2 = p_neighbour_2->rGetContainingElementIndices();

            //                     // Find common elements
            //                     std::set<unsigned> shared_elements;
            //                     std::set_intersection(containing_element_indices_1.begin(),
            //                                     containing_element_indices_1.end(),
            //                                     containing_element_indices_2.begin(),
            //                                     containing_element_indices_2.end(),
            //                                     std::inserter(shared_elements, shared_elements.begin()));

            //                     if(shared_elements.size()==1)
            //                     {
            //                         is_boundary = true;
            //                         unsigned neigh_1 = boundary_neighbours[0];
            //                         unsigned neigh_2 = boundary_neighbours[1];

            //                         std::set<unsigned> node_neighbours_1 = p_cell_population->GetNeighbouringNodeIndices(neigh_1);
            //                         std::set<unsigned> node_neighbours_2 = p_cell_population->GetNeighbouringNodeIndices(neigh_2);
                                    
            //                         // if neighbours are also eachothers neighbours, we have a void!
            //                         if((std::count(node_neighbours_1.begin(), node_neighbours_1.end(), neigh_2)) && (std::count(node_neighbours_2.begin(), node_neighbours_2.end(), neigh_1)) )
            //                         {
            //                             is_def_bound = true;
            //                         }

            //                     }
            //                     else if(containing_element_indices_1.size()>=2 || containing_element_indices_2.size()>=2)
            //                     {
            //                         is_boundary = true;
            //                     }
            //                 }
            //                 if(is_def_bound==true)
            //                 {
            //                     is_boundary = true;
            //                 }
                            
            //             }
            //         }
            //         // else
            //         // {
            //         //     unsigned node_index = node_iter->GetIndex();
            //         //     std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);
            //         //     unsigned number_of_boundary_neighbours = 0;
            //         //     std::vector<unsigned> boundary_neighbours;

            //         //     for (std::set<unsigned>::iterator neighbour_iter = node_neighbours.begin();
            //         //         neighbour_iter != node_neighbours.end();
            //         //         ++neighbour_iter)
            //         //     {
            //         //         unsigned neighbour_index = *neighbour_iter;
            //         //         Node<DIM>* p_neighbour_node = p_mesh->GetNode(neighbour_index);

            //         //         if(p_neighbour_node->IsBoundaryNode())
            //         //         {
            //         //             boundary_neighbours.push_back(neighbour_index);
            //         //             number_of_boundary_neighbours++;
            //         //             std::set<unsigned> containing_element_indices_neighbour = p_neighbour_node->rGetContainingElementIndices();
            //         //             if(containing_element_indices_neighbour.size() == 1)
            //         //             {
            //         //                 is_def_bound = true;
            //         //             }

            //         //         }
            //         //     }
            //         //     // if(number_of_boundary_neighbours == node_neighbours.size())
            //         //     if(number_of_boundary_neighbours == 2)
            //         //     {
            //         //         is_boundary = true;
            //         //         if(boundary_neighbours.size()==2)
            //         //         {
            //         //             Node<DIM>* p_neighbour_1 = p_mesh->GetNode(boundary_neighbours[0]);
            //         //             Node<DIM>* p_neighbour_2 = p_mesh->GetNode(boundary_neighbours[1]);

            //         //             std::set<unsigned> containing_element_indices_1 = p_neighbour_1->rGetContainingElementIndices();
            //         //             std::set<unsigned> containing_element_indices_2 = p_neighbour_2->rGetContainingElementIndices();

            //         //             // Find common elements
            //         //             std::set<unsigned> shared_elements;
            //         //             std::set_intersection(containing_element_indices_1.begin(),
            //         //                             containing_element_indices_1.end(),
            //         //                             containing_element_indices_2.begin(),
            //         //                             containing_element_indices_2.end(),
            //         //                             std::inserter(shared_elements, shared_elements.begin()));

            //         //             if(shared_elements.size()==1)
            //         //             {
            //         //                 unsigned neigh_1 = boundary_neighbours[0];
            //         //                 unsigned neigh_2 = boundary_neighbours[1];

            //         //                 std::set<unsigned> node_neighbours_1 = p_cell_population->GetNeighbouringNodeIndices(neigh_1);
            //         //                 std::set<unsigned> node_neighbours_2 = p_cell_population->GetNeighbouringNodeIndices(neigh_2);

            //         //                 // if neighbours are also eachothers neighbours, we have a void!
            //         //                 if((std::count(node_neighbours_1.begin(), node_neighbours_1.end(), neigh_2)) && (std::count(node_neighbours_2.begin(), node_neighbours_2.end(), neigh_1)) )
            //         //                 {
            //         //                     is_def_bound = true;
            //         //                 }

            //         //             }
            //         //         }
            //         //     }

            //         //     if(is_boundary==true && containing_element_indices.size()==3 && number_of_boundary_neighbours==(node_neighbours.size()-1))
            //         //     {
            //         //         if(is_def_bound == false)
            //         //         {
            //         //             is_boundary = false;
            //         //         }
                            
            //         //     }
            //         // }

                    
            //     }
            //     else
            //     {
            //         std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
            //         if(containing_element_indices.size() == 1)
            //         {
            //             is_boundary = true;
            //         }
            //         else if(containing_element_indices.size() > 1)
            //         {
            //             unsigned node_index = node_iter->GetIndex();
            //             std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);
        
            //             for (std::set<unsigned>::iterator neighbour_iter = node_neighbours.begin();
            //                 neighbour_iter != node_neighbours.end();
            //                 ++neighbour_iter)
            //             {
            //                 unsigned neighbour_index = *neighbour_iter;
            //                 Node<DIM>* p_neighbour_node = p_mesh->GetNode(neighbour_index);
            //                 std::set<unsigned> element_index_set_n = p_neighbour_node->rGetContainingElementIndices();

            //                 if(p_neighbour_node->IsBoundaryNode())
            //                 {
            //                     std::vector<unsigned> boundary_neighbours;
            //                     if( (element_index_set_n.size()==1) )
            //                     {
            //                         is_boundary = true;
            //                     }
                                
            //                 }
            //             }
            //         }
            //     }
            //     if(is_def_bound)
            //     {
            //         is_boundary = true;
            //     }

            //     node_iter->SetAsBoundaryNode(is_boundary);
            // }
        


    return ReCheck_Mesh;
}


template<unsigned DIM>
bool JaggedVertexEdgesModifierSimplified<DIM>::CheckForInternalNode(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here)
{
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));
    
    bool ReCheck_Mesh = false;
    bool performed_edge_modifier = false;
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
        if(true)
        {
            for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
                node_iter != p_mesh->GetNodeIteratorEnd();
                ++node_iter)
            {
                // if (node_iter->IsBoundaryNode())
                // {
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
                                // PRINT_VECTOR(p_node->rGetLocation())
                                // PRINT_VARIABLE(p_node->IsBoundaryNode());

                                std::set<unsigned>::const_iterator elem_it = containing_element_indices.begin();

                                unsigned elem_index_1 = (*elem_it);
                                VertexElement<DIM,DIM>* p_element_1 = p_mesh->GetElement(elem_index_1);
                                p_element_1->DeleteNode(p_element_1->GetNodeLocalIndex(node_index));

                                elem_it++;
                                unsigned elem_index_2 = (*elem_it);
                                VertexElement<DIM,DIM>* p_element_2 = p_mesh->GetElement(elem_index_2);
                                p_element_2->DeleteNode(p_element_2->GetNodeLocalIndex(node_index));

                                p_mesh->DeleteNodePriorToReMesh(node_index);

                                // Check if the boundary nodes need relabeling
                                if(p_neighbour_node_1->IsBoundaryNode())
                                {
                                    std::set<unsigned> node_neighbours_1 = p_cell_population->GetNeighbouringNodeIndices(node_neighbour_1);
                                    bool is_neighbour_node_1_boundary = false;
                                    unsigned numb_boundary_neighs_1 = 0;
                                    for (std::set<unsigned>::iterator neighbour_iter = node_neighbours_1.begin();
                                                            neighbour_iter != node_neighbours_1.end();
                                                            ++neighbour_iter)
                                    {
                                        unsigned node_neigh_index = *neighbour_iter;
                                        Node<DIM>* p_node_iter = p_mesh->GetNode(node_neigh_index);
                                        std::set<unsigned> containing_element_indices_neigh = p_node_iter->rGetContainingElementIndices();

                                        if(p_node_iter->IsBoundaryNode() && node_neigh_index != node_index && containing_element_indices_neigh.size()==1)
                                        {
                                            // PRINT_VECTOR(p_node_iter->rGetLocation());
                                            is_neighbour_node_1_boundary = true;
                                        }
                                        if(p_node_iter->IsBoundaryNode())
                                        {
                                            numb_boundary_neighs_1++;
                                        }
                                    }
                                    if(is_neighbour_node_1_boundary==false && numb_boundary_neighs_1==(node_neighbours_1.size()-1))
                                    {
                                        is_neighbour_node_1_boundary = true;
                                    }
                                    else if(is_neighbour_node_1_boundary==false && numb_boundary_neighs_1==3)
                                    {
                                        is_neighbour_node_1_boundary = true;
                                    }
                                    p_neighbour_node_1->SetAsBoundaryNode(is_neighbour_node_1_boundary);
                                // PRINT_VARIABLE(is_neighbour_node_1_boundary);
                                }
                                if(p_neighbour_node_2->IsBoundaryNode())
                                {
                                    std::set<unsigned> node_neighbours_2 = p_cell_population->GetNeighbouringNodeIndices(node_neighbour_2);
                                    bool is_neighbour_node_2_boundary = false;
                                    unsigned numb_boundary_neighs_2 = 0;
                                    for (std::set<unsigned>::iterator neighbour_iter = node_neighbours_2.begin();
                                                            neighbour_iter != node_neighbours_2.end();
                                                            ++neighbour_iter)
                                    {
                                        unsigned node_neigh_index = *neighbour_iter;
                                        Node<DIM>* p_node_iter = p_mesh->GetNode(node_neigh_index);
                                        std::set<unsigned> containing_element_indices_neigh = p_node_iter->rGetContainingElementIndices();

                                        if(p_node_iter->IsBoundaryNode() && node_neigh_index != node_index && containing_element_indices_neigh.size()==1)
                                        {
                                            // PRINT_VECTOR(p_node_iter->rGetLocation());
                                            is_neighbour_node_2_boundary = true;
                                        }
                                        if(p_node_iter->IsBoundaryNode())
                                        {
                                            numb_boundary_neighs_2++;
                                        }
                                    }
                                    if(is_neighbour_node_2_boundary==false && numb_boundary_neighs_2==(node_neighbours_2.size()-1))
                                    {
                                        is_neighbour_node_2_boundary = true;
                                    }
                                    else if(is_neighbour_node_2_boundary==false && numb_boundary_neighs_2==3)
                                    {
                                        is_neighbour_node_2_boundary = true;
                                    }
                                    p_neighbour_node_2->SetAsBoundaryNode(is_neighbour_node_2_boundary);
                                // PRINT_VARIABLE(is_neighbour_node_2_boundary);
                                }
                                performed_edge_modifier = true;
                                // TRACE("remove internal node");
                            }
                        }
                    }
                // }
            }
        }
        if(performed_edge_modifier)
        {
            ReCheck_Mesh = true;
            p_mesh->ReMesh();

        }
        performed_edge_modifier = false;

    // unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    // std::stringstream time;
    // time << num_timesteps;
    // std::stringstream t_numb_times_here;
    // t_numb_times_here << numb_times_here;
    // VertexMeshWriter<DIM,DIM> vertexmesh_writer1("tmp", "Mesh_6", false);
    // vertexmesh_writer1.WriteVtkUsingMesh(*p_mesh, time.str() + "_" + t_numb_times_here.str());
    return ReCheck_Mesh;
}

template<unsigned DIM>
bool SmoothVertexEdgesModifierSimplified<DIM>::CheckForNodeEdgeInt(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here)
{
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));
    
    bool ReCheck_Mesh = false;
    bool performed_edge_modifier = false;
 
        
    double mDistanceFromNodeToEdge = p_mesh->GetCellRearrangementThreshold();
    // double mDistanceFromNodeToEdge = 0.02;
    double mDistanceFromNodeToNodeCheck = 1.0;

    /*  | (v) |                       |(v)|
    *   oD    oB                      oD  oB
    *    \    |  Cell                   \ |  Cell
    *     \   |                   Cell   \|
    * Cell \  |                    _______oA
    * _____Ao |         ---->             |
    *         |                           |
    *    (v)  |                       (v) |
    *         oC                          oC
    */
    // push a node into an edge of a neighbouring cell...

    for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
        node_iter != p_mesh->GetNodeIteratorEnd();
        ++node_iter)
    {
        if (node_iter->IsBoundaryNode())
        {
            std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();

            if(containing_element_indices.size() >= 1)
            {
                unsigned node_1_index = node_iter->GetIndex();

                for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter_2 = p_mesh->GetNodeIteratorBegin();
                    node_iter_2 != p_mesh->GetNodeIteratorEnd();
                     ++node_iter_2)
                {
                    unsigned node_2_index = node_iter_2->GetIndex();
                    // Make sure node_2 is a boundary node

                    if (node_iter_2->IsBoundaryNode() && node_1_index!=node_2_index)
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
                                    if(norm_2(p_mesh->GetVectorFromAtoB(r_node_1,r_node_2)) < mDistanceFromNodeToNodeCheck)
                                    {
                                        // unsigned node_2_index = node_iter_2->GetIndex();
                                        std::set<unsigned> node_2_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_2_index);

                                        for (std::set<unsigned>::iterator neighbour_iter = node_2_neighbours.begin();
                                            neighbour_iter != node_2_neighbours.end();
                                            ++neighbour_iter)
                                        {
                                            unsigned neighbour_index = *neighbour_iter;
                                            Node<DIM>* p_neighbour_node = p_mesh->GetNode(neighbour_index);
                                            // std::set<unsigned> element_index_set_n = p_neighbour_node->rGetContainingElementIndices();

                                            if(p_neighbour_node->IsBoundaryNode() && neighbour_index!= node_1_index)
                                            {
                                                
                                                c_vector<double, DIM> r_neighbour = p_neighbour_node->rGetLocation();
                                                                                            
                                                double distance_to_edge = DistanceToEdgeFromPoint(r_neighbour,r_node_2,r_node_1);
                                                
                                                /*             \
                                                *               \
                                                *                 o-----
                                                *                      
                                                *      o----------o
                                                *
                                                * Need a way to implement this threshold:
                                                */
                                                // if(distance_to_edge < mDistanceFromNodeToEdge  && norm_2(p_mesh->GetVectorFromAtoB(r_node_1,r_neighbour)) < norm_2(p_mesh->GetVectorFromAtoB(r_node_2,r_neighbour)))
                                                // {
                                                //     PRINT_VARIABLE(SimulationTime::Instance()->GetTime());
                                                //     TRACE("Doinging the thing");
                                                //     PRINT_VECTOR(r_node_1);

                                                //     Node<DIM>* p_neighbour_1= p_mesh->GetNode(node_1_index);
                                                //     Node<DIM>* p_neighbour_2= p_mesh->GetNode(node_2_index);


                                                    
                                                //     p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);

                                                //     // Check if p_neighbour_node is still a boundary node
                                                    
                                                //     unsigned p_neighbour_index = p_neighbour_node->GetIndex();
                                                //     std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(p_neighbour_index);
                                                //     bool is_boundary_neigh = false;
                                                //     unsigned numb_boundary_neigh = 0;
                                                //     for (std::set<unsigned>::iterator neighbour_iter = node_neighbours.begin();
                                                //         neighbour_iter != node_neighbours.end();
                                                //         ++neighbour_iter)
                                                //     {
                                                //         unsigned node_neigh = *neighbour_iter;
                                                //         if(node_neigh != node_2_index)
                                                //         {
                                                //             Node<DIM>* p_neighbour= p_mesh->GetNode(node_neigh);
                                                //             std::set<unsigned> containing_elements_neigh = p_neighbour->rGetContainingElementIndices();

                                                //             if(p_neighbour->IsBoundaryNode())
                                                //             {
                                                //                 numb_boundary_neigh++ ;
                                                //             }

                                                //             if(p_neighbour->IsBoundaryNode() && containing_elements_neigh.size() == 1 && node_neigh != node_1_index)
                                                //             {
                                                //                 is_boundary_neigh = true;
                                                //             }
                                                //             else if(p_neighbour->IsBoundaryNode() && containing_elements_neigh.size() == 2)
                                                //             {
                                                //                 std::set<unsigned> containing_elements_iter = p_neighbour_node->rGetContainingElementIndices();
                                                //                 std::set<unsigned> shared_elements_2;
                                                //                 std::set_intersection(containing_elements_iter.begin(),
                                                //                                     containing_elements_iter.end(),
                                                //                                     containing_elements_neigh.begin(),
                                                //                                     containing_elements_neigh.end(),
                                                //                                     std::inserter(shared_elements_2, shared_elements_2.begin()));
                                                //                 if(shared_elements_2.size() == 2 && containing_elements_iter.size()==2)
                                                //                 {
                                                //                     is_boundary_neigh = true;
                                                //                 }
                                                //             }

                                                //         }
                                                //     }
                                                //     // This is a stiched edge, which can somethimes happen.
                                                //     if(is_boundary_neigh==false && numb_boundary_neigh == node_2_neighbours.size())
                                                //     {
                                                //         is_boundary_neigh = true;
                                                //     }
                                                //     p_neighbour_node->SetAsBoundaryNode(is_boundary_neigh);
                                                //     // PRINT_VARIABLE(is_boundary_neigh);
                                                //     // PRINT_VECTOR(p_neighbour_node->rGetLocation());

                                                    
                                                    
                                                    
                                                //     // p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);
                                                    
                                                //     // unsigned node_index_1 = p_neighbour_1->GetIndex(); 
                                                //     // std::set<unsigned> node_1_neighbs = p_cell_population->GetNeighbouringNodeIndices(node_1_index);
                                                //     // bool is_boundary = false;
                                                //     // unsigned numb_boundary_nodes = 0;
                                                //     // for (std::set<unsigned>::iterator neighbour_iter = node_1_neighbs.begin();
                                                //     //     neighbour_iter != node_1_neighbs.end();
                                                //     //     ++neighbour_iter)
                                                //     // {
                                                //     //     unsigned node_neigh = *neighbour_iter;
                                                //     //     if(node_neigh != node_2_index)
                                                //     //     {
                                                //     //         Node<DIM>* p_neighbour_1= p_mesh->GetNode(node_neigh);
                                                //     //         std::set<unsigned> containing_elements_neigh = p_neighbour_1->rGetContainingElementIndices();

                                                //     //         if(p_neighbour_1->IsBoundaryNode())
                                                //     //         {
                                                //     //             numb_boundary_nodes++ ;
                                                //     //         }

                                                //     //         if(p_neighbour_1->IsBoundaryNode() && containing_elements_neigh.size() == 1)
                                                //     //         {
                                                //     //             is_boundary = true;
                                                //     //         }
                                                //     //         else if(p_neighbour_1->IsBoundaryNode() && containing_elements_neigh.size() == 2)
                                                //     //         {
                                                //     //             std::set<unsigned> containing_elements_iter = node_iter->rGetContainingElementIndices();
                                                //     //             std::set<unsigned> shared_elements_2;
                                                //     //             std::set_intersection(containing_elements_iter.begin(),
                                                //     //                                 containing_elements_iter.end(),
                                                //     //                                 containing_elements_neigh.begin(),
                                                //     //                                 containing_elements_neigh.end(),
                                                //     //                                 std::inserter(shared_elements_2, shared_elements_2.begin()));
                                                //     //             if(shared_elements_2.size() == 2)
                                                //     //             {
                                                //     //                 is_boundary = true;
                                                //     //             }
                                                //     //         }
                                                //     //     }
                                                //     // }
                                                //     // // This is a stiched edge, which can somethimes happen.
                                                //     // if(is_boundary==false && numb_boundary_nodes == node_2_neighbours.size())
                                                //     // {
                                                //     //     is_boundary = true;
                                                //     // }

                                                //     // Node<DIM>* p_neighbour_1= p_mesh->GetNode(node_1_index);
                                                //     // Node<DIM>* p_neighbour_2= p_mesh->GetNode(node_2_index);
                                                    
                                                //     // p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);

                                                //     p_neighbour_1->SetAsBoundaryNode(true);

                                                //     TRACE("Done");
                                                    
                                                //     // PRINT_VARIABLE(is_boundary);
                                                //     // PRINT_VECTOR(p_neighbour_1->rGetLocation());
                                                //     performed_edge_modifier = true;
                                                //     // TRACE("Performed node merge into edge");
                                                //     break;
                                                // }
                                                double length_of_edge = norm_2(p_mesh->GetVectorFromAtoB(r_node_2,r_neighbour));
                                                
                                                if(distance_to_edge < mDistanceFromNodeToEdge  && norm_2(p_mesh->GetVectorFromAtoB(r_node_1,r_neighbour)) < length_of_edge && norm_2(p_mesh->GetVectorFromAtoB(r_node_1,r_node_2)) < length_of_edge)
                                                {

                                                    // PRINT_VARIABLE(SimulationTime::Instance()->GetTime());
                                                    TRACE("Got a node too close to an edge!!!");
                                                    PRINT_VECTOR(r_node_1);
                                                    PRINT_VECTOR(r_neighbour);
                                                    PRINT_VECTOR(r_node_2);

                                                    Node<DIM>* p_neighbour_1= p_mesh->GetNode(node_1_index);
                                                    Node<DIM>* p_neighbour_2= p_mesh->GetNode(node_2_index);

                                                    // p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);

                                                    // Find the indices of the elements owned by each node
                                                    std::set<unsigned> elements_containing_nodeA = p_neighbour_2->rGetContainingElementIndices();
                                                    std::set<unsigned> elements_containing_nodeB = p_neighbour_node->rGetContainingElementIndices();

                                                    // Find common elements
                                                    std::set<unsigned> shared_elements;
                                                    std::set_intersection(elements_containing_nodeA.begin(),
                                                                        elements_containing_nodeA.end(),
                                                                        elements_containing_nodeB.begin(),
                                                                        elements_containing_nodeB.end(),
                                                                        std::inserter(shared_elements, shared_elements.begin()));
                                                    
                                                    // Iterate over common elements
                                                    unsigned node_A_index = p_neighbour_2->GetIndex();
                                                    unsigned node_B_index = p_neighbour_node->GetIndex();
                                                    for (std::set<unsigned>::iterator iter = shared_elements.begin();
                                                        iter != shared_elements.end();
                                                        ++iter)
                                                    {
                                                        VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(*iter);

                                                        // Find which node has the lower local index in this element
                                                        unsigned local_indexA = p_element->GetNodeLocalIndex(node_A_index);
                                                        unsigned local_indexB = p_element->GetNodeLocalIndex(node_B_index);

                                                        unsigned index = local_indexB;

                                                        // If node B has a higher index then use node A's index...
                                                        if (local_indexB > local_indexA)
                                                        {
                                                            index = local_indexA;

                                                            // ...unless nodes A and B share the element's last edge
                                                            if ((local_indexA == 0) && (local_indexB == p_element->GetNumNodes()-1))
                                                            {
                                                                index = local_indexB;
                                                            }
                                                        }
                                                        else if ((local_indexB == 0) && (local_indexA == p_element->GetNumNodes()-1))
                                                        {
                                                            // ...otherwise use node B's index, unless nodes A and B share the element's last edge
                                                            index = local_indexA;
                                                        }

                                                        // Add new node to this element
                                                        p_mesh->GetElement(*iter)->AddNode(p_neighbour_1, index);
                                                    }


                                                    // Check if p_neighbour_node is still a boundary node
                                                    
                                                    unsigned p_neighbour_index = p_neighbour_node->GetIndex();
                                                    std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(p_neighbour_index);
                                                    bool is_boundary_neigh = false;
                                                    unsigned numb_boundary_neigh = 0;
                                                    for (std::set<unsigned>::iterator neighbour_iter = node_neighbours.begin();
                                                        neighbour_iter != node_neighbours.end();
                                                        ++neighbour_iter)
                                                    {
                                                        unsigned node_neigh = *neighbour_iter;
                                                        if(node_neigh != node_2_index)
                                                        {
                                                            Node<DIM>* p_neighbour= p_mesh->GetNode(node_neigh);
                                                            std::set<unsigned> containing_elements_neigh = p_neighbour->rGetContainingElementIndices();

                                                            if(p_neighbour->IsBoundaryNode())
                                                            {
                                                                numb_boundary_neigh++ ;
                                                            }

                                                            if(p_neighbour->IsBoundaryNode() && containing_elements_neigh.size() == 1 && node_neigh != node_1_index)
                                                            {
                                                                is_boundary_neigh = true;
                                                            }
                                                            else if(p_neighbour->IsBoundaryNode() && containing_elements_neigh.size() == 2)
                                                            {
                                                                std::set<unsigned> containing_elements_iter = p_neighbour_node->rGetContainingElementIndices();
                                                                std::set<unsigned> shared_elements_2;
                                                                std::set_intersection(containing_elements_iter.begin(),
                                                                                    containing_elements_iter.end(),
                                                                                    containing_elements_neigh.begin(),
                                                                                    containing_elements_neigh.end(),
                                                                                    std::inserter(shared_elements_2, shared_elements_2.begin()));
                                                                if(shared_elements_2.size() == 2 && containing_elements_iter.size()==2)
                                                                {
                                                                    is_boundary_neigh = true;
                                                                }
                                                            }

                                                        }
                                                    }
                                                    // This is a stiched edge, which can somethimes happen.
                                                    if(is_boundary_neigh==false && numb_boundary_neigh == node_2_neighbours.size())
                                                    {
                                                        is_boundary_neigh = true;
                                                    }
                                                    p_neighbour_node->SetAsBoundaryNode(is_boundary_neigh);
                                                    

                                                    p_neighbour_1->SetAsBoundaryNode(true);

                                                    // TRACE("Done");
                                                    
                                                    performed_edge_modifier = true;
                                                    TRACE("Performed node merge into edge");
                                                    break;
                                                }

                                            }
                                        }

                                    }

                        }
                        else if(node_2_and_node_in_same_element)
                        {
                            c_vector<double, DIM> r_node_1 = node_iter->rGetLocation();
                            c_vector<double, DIM> r_node_2 = node_iter_2->rGetLocation();

                            
                            

                                    // 1 and 2 are close, check if 1 is close to the edge of node 2 and its neighbours
                                    if(norm_2(p_mesh->GetVectorFromAtoB(r_node_1,r_node_2)) < mDistanceFromNodeToNodeCheck)
                                    {
                                        // unsigned node_2_index = node_iter_2->GetIndex();
                                        std::set<unsigned> node_2_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_2_index);

                                        for (std::set<unsigned>::iterator neighbour_iter = node_2_neighbours.begin();
                                            neighbour_iter != node_2_neighbours.end();
                                            ++neighbour_iter)
                                        {
                                            unsigned neighbour_index = *neighbour_iter;
                                            Node<DIM>* p_neighbour_node = p_mesh->GetNode(neighbour_index);
                                            // std::set<unsigned> element_index_set_n = p_neighbour_node->rGetContainingElementIndices();

                                            if(p_neighbour_node->IsBoundaryNode() && neighbour_index!= node_1_index)
                                            {
                                                
                                                c_vector<double, DIM> r_neighbour = p_neighbour_node->rGetLocation();
                                                                                            
                                                double distance_to_edge = DistanceToEdgeFromPoint(r_neighbour,r_node_2,r_node_1);
                                                
                                                /*             \
                                                *               \
                                                *                 o-----
                                                *                      
                                                *      o----------o
                                                *
                                                * Need a way to implement this threshold:
                                                */
                                                double length_of_edge = norm_2(p_mesh->GetVectorFromAtoB(r_node_2,r_neighbour));
                                                
                                                if(distance_to_edge < mDistanceFromNodeToEdge  && norm_2(p_mesh->GetVectorFromAtoB(r_node_1,r_neighbour)) < length_of_edge && norm_2(p_mesh->GetVectorFromAtoB(r_node_1,r_node_2)) < length_of_edge)
                                                {
                                                    PRINT_VECTOR(r_node_1);
                                                    PRINT_VECTOR(r_node_2);
                                                    PRINT_VECTOR(r_neighbour);

                                                    Node<DIM>* p_neighbour_1= p_mesh->GetNode(node_1_index);
                                                    Node<DIM>* p_neighbour_2= p_mesh->GetNode(node_2_index);

                                                    // p_mesh->PerformNodeMerge(p_neighbour_1,p_neighbour_2);

                                                    // Find the indices of the elements owned by each node
                                                    std::set<unsigned> elements_containing_nodeA = p_neighbour_2->rGetContainingElementIndices();
                                                    std::set<unsigned> elements_containing_nodeB = p_neighbour_node->rGetContainingElementIndices();

                                                    // Find common elements
                                                    std::set<unsigned> shared_elements;
                                                    std::set_intersection(elements_containing_nodeA.begin(),
                                                                        elements_containing_nodeA.end(),
                                                                        elements_containing_nodeB.begin(),
                                                                        elements_containing_nodeB.end(),
                                                                        std::inserter(shared_elements, shared_elements.begin()));
                                                    
                                                    // Iterate over common elements
                                                    unsigned node_A_index = p_neighbour_2->GetIndex();
                                                    unsigned node_B_index = p_neighbour_node->GetIndex();
                                                    for (std::set<unsigned>::iterator iter = shared_elements.begin();
                                                        iter != shared_elements.end();
                                                        ++iter)
                                                    {
                                                        VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(*iter);

                                                        // Find which node has the lower local index in this element
                                                        unsigned local_indexA = p_element->GetNodeLocalIndex(node_A_index);
                                                        unsigned local_indexB = p_element->GetNodeLocalIndex(node_B_index);

                                                        unsigned index = local_indexB;

                                                        // If node B has a higher index then use node A's index...
                                                        if (local_indexB > local_indexA)
                                                        {
                                                            index = local_indexA;

                                                            // ...unless nodes A and B share the element's last edge
                                                            if ((local_indexA == 0) && (local_indexB == p_element->GetNumNodes()-1))
                                                            {
                                                                index = local_indexB;
                                                            }
                                                        }
                                                        else if ((local_indexB == 0) && (local_indexA == p_element->GetNumNodes()-1))
                                                        {
                                                            // ...otherwise use node B's index, unless nodes A and B share the element's last edge
                                                            index = local_indexA;
                                                        }

                                                        // Add new node to this element
                                                        p_mesh->GetElement(*iter)->AddNode(p_neighbour_1, index);
                                                    }


                                                    // Check if p_neighbour_node is still a boundary node
                                                    
                                                    unsigned p_neighbour_index = p_neighbour_node->GetIndex();
                                                    std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(p_neighbour_index);
                                                    bool is_boundary_neigh = false;
                                                    unsigned numb_boundary_neigh = 0;
                                                    for (std::set<unsigned>::iterator neighbour_iter = node_neighbours.begin();
                                                        neighbour_iter != node_neighbours.end();
                                                        ++neighbour_iter)
                                                    {
                                                        unsigned node_neigh = *neighbour_iter;
                                                        if(node_neigh != node_2_index)
                                                        {
                                                            Node<DIM>* p_neighbour= p_mesh->GetNode(node_neigh);
                                                            std::set<unsigned> containing_elements_neigh = p_neighbour->rGetContainingElementIndices();

                                                            if(p_neighbour->IsBoundaryNode())
                                                            {
                                                                numb_boundary_neigh++ ;
                                                            }

                                                            if(p_neighbour->IsBoundaryNode() && containing_elements_neigh.size() == 1 && node_neigh != node_1_index)
                                                            {
                                                                is_boundary_neigh = true;
                                                            }
                                                            else if(p_neighbour->IsBoundaryNode() && containing_elements_neigh.size() == 2)
                                                            {
                                                                std::set<unsigned> containing_elements_iter = p_neighbour_node->rGetContainingElementIndices();
                                                                std::set<unsigned> shared_elements_2;
                                                                std::set_intersection(containing_elements_iter.begin(),
                                                                                    containing_elements_iter.end(),
                                                                                    containing_elements_neigh.begin(),
                                                                                    containing_elements_neigh.end(),
                                                                                    std::inserter(shared_elements_2, shared_elements_2.begin()));
                                                                if(shared_elements_2.size() == 2 && containing_elements_iter.size()==2)
                                                                {
                                                                    is_boundary_neigh = true;
                                                                }
                                                            }

                                                        }
                                                    }
                                                    // This is a stiched edge, which can somethimes happen.
                                                    if(is_boundary_neigh==false && numb_boundary_neigh == node_2_neighbours.size())
                                                    {
                                                        is_boundary_neigh = true;
                                                    }
                                                    p_neighbour_node->SetAsBoundaryNode(is_boundary_neigh);
                                                    

                                                    p_neighbour_1->SetAsBoundaryNode(true);

                                                    // TRACE("Done");
                                                    
                                                    performed_edge_modifier = true;
                                                    TRACE("Performed node merge into edge - second imp");
                                                    break;
                                                }

                                            }
                                        }

                                    }
                            
                        }


                    }
                    
                    if(performed_edge_modifier)
                    {
                        break;
                    }

                }
            }
        }
        if(performed_edge_modifier)
        {
            break;
        }
    }
        

    if(performed_edge_modifier)
    {
        TRACE("Let's Mesh");
        ReCheck_Mesh = true;
        p_mesh->ReMesh();
        TRACE("ReMesh Completed");

    }

    performed_edge_modifier = false;

    // unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    // std::stringstream time;
    // time << num_timesteps;
    // std::stringstream t_numb_times_here;
    // t_numb_times_here << numb_times_here;
    // VertexMeshWriter<DIM,DIM> vertexmesh_writer1("tmp", "Mesh_5", false);
    // vertexmesh_writer1.WriteVtkUsingMesh(*p_mesh, time.str() + "_" + t_numb_times_here.str());

    return ReCheck_Mesh;

}


template<unsigned DIM>
bool JaggedVertexEdgesModifierSimplified<DIM>::CheckT3Distrupt(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here)
{
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));
    
    bool ReCheck_Mesh = false;
    bool performed_edge_modifier = false;
    
    double distanceBetweenVerteciesThreshold_2 = 2.0*p_mesh->GetCellRearrangementThreshold();
    double distrubThreshold = 0.1;
    // double distanceBetweenVerteciesThreshold_2 = 0.05;

    for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter_1 = p_mesh->GetNodeIteratorBegin();
        node_iter_1 != p_mesh->GetNodeIteratorEnd();
        ++node_iter_1)
    {
        std::set<unsigned> containing_element_indices_1 = node_iter_1->rGetContainingElementIndices();

        if (node_iter_1->IsBoundaryNode() && containing_element_indices_1.size() == 2)
        {

            for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter_2 = p_mesh->GetNodeIteratorBegin();
            node_iter_2 != p_mesh->GetNodeIteratorEnd();
            ++node_iter_2)
            {
                unsigned node_index_1 = node_iter_1->GetIndex();
                unsigned node_index_2 = node_iter_2->GetIndex();
                std::set<unsigned> containing_element_indices_2 = node_iter_2->rGetContainingElementIndices();

                if (node_iter_2->IsBoundaryNode() && containing_element_indices_2.size() == 2 && node_index_1 != node_index_2 )
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
                        r_node_12 = r_node_12/norm_2(r_node_12);
                        c_vector<double, DIM> r_orth_r_node_12 = zero_vector<double>(DIM);
                        r_orth_r_node_12[0] = r_node_12[1];
                        r_orth_r_node_12[1] = r_node_12[0];

                        p_node_1->rGetModifiableLocation()[0] = r_node_1[0] + distrubThreshold*r_orth_r_node_12[0];
                        p_node_1->rGetModifiableLocation()[1] = r_node_1[1] + distrubThreshold*r_orth_r_node_12[1];

                        p_node_2->rGetModifiableLocation()[0] = r_node_2[0] - distrubThreshold*r_orth_r_node_12[0];
                        p_node_2->rGetModifiableLocation()[1] = r_node_2[1] - distrubThreshold*r_orth_r_node_12[1];

                        // p_mesh->PerformNodeMerge(p_node_1,p_node_2);

                        // p_node_1->SetAsBoundaryNode(true);
                        performed_edge_modifier = true;
                        break;
                    }
                }
            }
        }
    }
        
        
    if(performed_edge_modifier)
    {
        ReCheck_Mesh = true;
        p_mesh->ReMesh();
    }

    performed_edge_modifier = false;

    return ReCheck_Mesh;

}

template<unsigned DIM>
void JaggedVertexEdgesModifierSimplified<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{

    // Call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class JaggedVertexEdgesModifierSimplified<1>;
template class JaggedVertexEdgesModifierSimplified<2>;
template class JaggedVertexEdgesModifierSimplified<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(JaggedVertexEdgesModifierSimplified)

