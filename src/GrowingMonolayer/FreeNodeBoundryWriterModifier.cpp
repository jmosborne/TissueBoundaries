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

#include "FreeNodeBoundryWriterModifier.hpp"
#include "VtkMeshWriter.hpp"

#include "VertexBasedCellPopulation.hpp"
#include "MutableVertexMesh.hpp"

#include "Cylindrical2dVertexMesh.hpp"

#include "Debug.hpp"



template<unsigned DIM>
FreeNodeBoundryWriterModifier<DIM>::FreeNodeBoundryWriterModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
    mOutputDirectory("")
{
}

template<unsigned DIM>
void FreeNodeBoundryWriterModifier<DIM>::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;
}

// std::string MeshModifier::GetOutputDirectory()
// {
// 	return mOutputDirectory;
// }

template<unsigned DIM>
FreeNodeBoundryWriterModifier<DIM>::~FreeNodeBoundryWriterModifier()
{
}

template<unsigned DIM>
void FreeNodeBoundryWriterModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
}

template<unsigned DIM>
void FreeNodeBoundryWriterModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));



    std::stringstream time;
    time << SimulationTime::Instance()->GetTimeStepsElapsed();

    VertexMeshWriter<DIM, DIM> mesh_writer(mOutputDirectory, "MY_results_"+time.str(), false);

    
    // mesh_writer.SetParallelFiles(*p_mesh);

    auto num_real_nodes = p_mesh->GetNumNodes();


    // For each cell that this process owns, find the corresponding node index, which we only want to calculate once
    std::vector<unsigned> node_indices_in_cell_order;
    std::vector<double> vtk_cell_data(num_real_nodes);
    

    // bool have_node = false;
    for (typename VertexMesh<DIM,DIM>::NodeIterator node_iter_original = p_mesh->GetNodeIteratorBegin();
        node_iter_original != p_mesh->GetNodeIteratorEnd();
        ++node_iter_original)
    {
        unsigned node_index = node_iter_original->GetIndex();

        node_indices_in_cell_order.emplace_back(node_index);
        vtk_cell_data[node_index] = node_iter_original->IsBoundaryNode();
    }
        
            
    mesh_writer.AddPointData("Is_Boundary_Node", vtk_cell_data);

    // mesh_writer.WriteFilesUsingMesh(*p_cylindrical_mesh);
    // mesh_writer.WriteVtkUsingMesh(*p_mesh, time.str());
    mesh_writer.WriteVtkUsingMesh(*p_mesh);


}

template<unsigned DIM>
void FreeNodeBoundryWriterModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateAtEndOfOutputTimeStep(rCellPopulation);

}


template<unsigned DIM>
void FreeNodeBoundryWriterModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
// template class FreeNodeBoundryWriterModifier<1>;
template class FreeNodeBoundryWriterModifier<2>;
// template class FreeNodeBoundryWriterModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FreeNodeBoundryWriterModifier)

