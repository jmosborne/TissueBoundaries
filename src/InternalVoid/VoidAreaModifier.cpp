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

#include "VoidAreaModifier.hpp"

#include "NodeBasedCellPopulation.hpp"

#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"

#include "VertexBasedCellPopulation.hpp"


#include "Debug.hpp"



template<unsigned DIM>
VoidAreaModifier<DIM>::VoidAreaModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
    mOutputDirectory(""),
    mCutoff(1.5),
    mPixelSeparation(0.05)
{
}

template<unsigned DIM>
void VoidAreaModifier<DIM>::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;

    OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
    out_stream locationFile = output_file_handler.OpenOutputFile("voidArea.dat");
    *locationFile << "time \t";
    *locationFile << "VoidArea" << "\t";
    *locationFile << "\n";
    locationFile->close();
}

template<unsigned DIM>
void VoidAreaModifier<DIM>::SetCutoff(double cutoff)
{
	mCutoff = cutoff;
}

template<unsigned DIM>
void VoidAreaModifier<DIM>::SetPixelSeparation(double pixelSeparation)
{
	mPixelSeparation = pixelSeparation;
}

// std::string MeshModifier::GetOutputDirectory()
// {
// 	return mOutputDirectory;
// }

template<unsigned DIM>
VoidAreaModifier<DIM>::~VoidAreaModifier()
{
}

template<unsigned DIM>
void VoidAreaModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
}

template<unsigned DIM>
void VoidAreaModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
    if(DIM == 2)
    {
        double void_area = 0.0;

        // both NodeBased and MeshBased (with no ghosts) populations are treated the same here, as both just have a cutoff radius
        if( (bool(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation)) 
        || bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)))
        && !(bool(dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))) )
        {
            double pixel_radial_reach = 0.5*mCutoff;
            double separation_between_pixels = mPixelSeparation;
            double area_of_pixel = pow(separation_between_pixels,2);

            // minus 2 is to ensure we are in the tissue and not detecting gaps due to periodicity from each side
            int pixel_tissue_width = (rCellPopulation.rGetMesh().GetWidth(0) - 2)/separation_between_pixels;
            int pixel_tissue_depth = (rCellPopulation.rGetMesh().GetWidth(1) - 2)/separation_between_pixels;

            int number_of_pixels_with_cell = 0;

            for(unsigned pixel_i = 0; pixel_i<pixel_tissue_width; pixel_i++)
            {
                for(unsigned pixel_j = 0; pixel_j<pixel_tissue_depth; pixel_j++)
                {
                    double pixel_x_coordinate = 1.0 + pixel_i*separation_between_pixels;
                    double pixel_y_coordinate = 1.0 + pixel_j*separation_between_pixels;
                    bool does_pixel_contain_cell = false;
                    // Iterate over cell population
                    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                        cell_iter != rCellPopulation.End();
                        ++cell_iter)
                    {
                        c_vector<double, DIM>  cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                        if( pow(cell_location[0] - pixel_x_coordinate,2) + pow(cell_location[1] - pixel_y_coordinate,2) < pow(pixel_radial_reach,2) )
                        {
                            does_pixel_contain_cell = true;
                            break;
                        }
                    }

                    if(!does_pixel_contain_cell)
                    {
                        void_area = void_area + area_of_pixel;
                    }
                }

            }

        }

        // if(bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)))
        // {

        // }
        
        // We handle Vertex model a little differently as there is a very natural cell area defined
        if(bool(dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation)) 
        || bool(dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation)))
        {
            double tissue_area = 0.0;

            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
            {
                // Get the volume of this cell
                double cell_volume = rCellPopulation.GetVolumeOfCell(*cell_iter);

                tissue_area = tissue_area + cell_volume;
            }

            void_area = (rCellPopulation.rGetMesh().GetWidth(0))*(rCellPopulation.rGetMesh().GetWidth(1)) - tissue_area;
            
        }

        OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
        out_stream locationFile = output_file_handler.OpenOutputFile("voidArea.dat", std::ios::app);
        SimulationTime* p_time = SimulationTime::Instance();
        *locationFile << p_time->GetTime() << "\t";
        *locationFile << void_area << "\t";
        *locationFile << "\n";
        locationFile->close();

    }

}



template<unsigned DIM>
void VoidAreaModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateAtEndOfOutputTimeStep(rCellPopulation);


}


template<unsigned DIM>
void VoidAreaModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class VoidAreaModifier<1>;
template class VoidAreaModifier<2>;
template class VoidAreaModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VoidAreaModifier)

