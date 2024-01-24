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

#include "CicularityCalcModifier.hpp"

#include "NodeBasedCellPopulation.hpp"

#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"

#include "VertexBasedCellPopulation.hpp"

#include <cmath>

#include "Debug.hpp"



template<unsigned DIM>
CicularityCalcModifier<DIM>::CicularityCalcModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
    mOutputDirectory(""),
    mCutoff(2.0),
    mPixelSeparation(0.25)
{
}

template<unsigned DIM>
void CicularityCalcModifier<DIM>::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;

    OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
    out_stream locationFile = output_file_handler.OpenOutputFile("CicularityCalc.dat");
    *locationFile << "time \t";
    *locationFile << "Cell_Area" << "\t";
    *locationFile << "Cell_Perimiter" << "\t";
    *locationFile << "Cell_Circularity" << "\t";
    *locationFile << "Number_Cells" << "\t";
    *locationFile << "\n";
    locationFile->close();


    OutputFileHandler output_file_handler_2(mOutputDirectory+"/", false);
    out_stream locationFile_2 = output_file_handler_2.OpenOutputFile("CircularityContour.dat");
    locationFile_2->close();
}

template<unsigned DIM>
void CicularityCalcModifier<DIM>::SetCutoff(double cutoff)
{
	mCutoff = cutoff;
}

template<unsigned DIM>
void CicularityCalcModifier<DIM>::SetPixelSeparation(double pixelSeparation)
{
	mPixelSeparation = pixelSeparation;
}

// std::string MeshModifier::GetOutputDirectory()
// {
// 	return mOutputDirectory;
// }

template<unsigned DIM>
CicularityCalcModifier<DIM>::~CicularityCalcModifier()
{
}

template<unsigned DIM>
void CicularityCalcModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
}

template<unsigned DIM>
void CicularityCalcModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
    if(DIM == 2)
    {
        // double pixel_area = 0.0;
        // double pixel_perimiter = 0.0;

        double tissue_area = 0.0;
        double tissue_perimiter = 0.0;
        unsigned number_cells = 0.0;

        OutputFileHandler output_file_handler_2(mOutputDirectory+"/", false);
        out_stream locationFile_2 = output_file_handler_2.OpenOutputFile("CircularityContour.dat", std::ios::app);
        SimulationTime* p_time = SimulationTime::Instance();
        *locationFile_2 << p_time->GetTime() << "\t";

        // both NodeBased and MeshBased (with no ghosts) populations are treated the same here, as both just have a cutoff radius
        // if( (bool(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation)) 
        // || bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)))
        // && !(bool(dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))) )

        {
            std::vector<unsigned> boundary_cells;
            std::vector<double> x_boundary;
            std::vector<double> y_boundary;

            if(bool(dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation)))
            {
                for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                    cell_iter != rCellPopulation.End();
                    ++cell_iter)
                {
                    cell_iter->GetCellData()->SetItem("is_boundary", 0.0);
                    number_cells++;
                }

                MeshBasedCellPopulationWithGhostNodes<DIM>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);
                MutableMesh<DIM,DIM>* p_mesh = static_cast<MutableMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));

                for (typename MutableMesh<DIM,DIM>::EdgeIterator edge_iterator = p_mesh->EdgesBegin();
                    edge_iterator != p_mesh->EdgesEnd();
                    ++edge_iterator)
                {
                    unsigned nodeA_global_index = edge_iterator.GetNodeA()->GetIndex();
                    unsigned nodeB_global_index = edge_iterator.GetNodeB()->GetIndex();


                    if(p_cell_population->IsGhostNode(nodeA_global_index))
                    {
                        if(!(p_cell_population->IsGhostNode(nodeB_global_index)))
                        {
                            if(!(std::count(boundary_cells.begin(), boundary_cells.end(), nodeB_global_index)))
                            {
                                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(nodeB_global_index);
                                c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(p_cell);

                                boundary_cells.push_back(nodeB_global_index);

                                x_boundary.push_back(real_node_location[0]);
                                y_boundary.push_back(real_node_location[1]);

                                p_cell->GetCellData()->SetItem("is_boundary", 1.0);

                            }

                        }


                    }
                    else if(p_cell_population->IsGhostNode(nodeB_global_index))
                    {
                        if(!(p_cell_population->IsGhostNode(nodeA_global_index)))
                        {
                            if(!(std::count(boundary_cells.begin(), boundary_cells.end(), nodeA_global_index)))
                            {
                                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(nodeA_global_index);
                                c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(p_cell);

                                boundary_cells.push_back(nodeA_global_index);

                                x_boundary.push_back(real_node_location[0]);
                                y_boundary.push_back(real_node_location[1]);

                                p_cell->GetCellData()->SetItem("is_boundary", 1.0);

                            }

                        }

                    }


                }


            }

            else if(bool(dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation)))
            {
                VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
                MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));

                for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                    cell_iter != rCellPopulation.End();
                    ++cell_iter)
                {
                    bool check_node = true;
                    cell_iter->GetCellData()->SetItem("is_boundary", 0.0);
                    number_cells++;

                    VertexElement<DIM, DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(*cell_iter);
                    for (unsigned i=0; i<p_element->GetNumNodes(); i++)
                    {

                        unsigned node_global_index = p_element->GetNodeGlobalIndex(i);

                        Node<DIM>* p_node_a = p_mesh->GetNode(node_global_index);

                        if (p_node_a->IsBoundaryNode() && check_node)
                        {
                            unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                            c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                            
                            boundary_cells.push_back(real_node_index);

                            x_boundary.push_back(real_node_location[0]);
                            y_boundary.push_back(real_node_location[1]);

                            cell_iter->GetCellData()->SetItem("is_boundary", 1.0);
                            check_node = false;
                        }
                    }




                }

            }
            
            else //NodeBasedCellPopulation, MeshBasedCellPopulation
            {

                double pixel_radial_reach = 0.5*mCutoff;
                double separation_between_pixels = mPixelSeparation;
                // double area_of_pixel = separation_between_pixels*separation_between_pixels;

                // minus 2 is to ensure we are in the tissue and not detecting gaps due to periodicity from each side
                double x_max = rCellPopulation.rGetMesh().GetWidth(0);
                double y_max = rCellPopulation.rGetMesh().GetWidth(1);
                
                int pixel_tissue_width = (2*(x_max) + 2)/separation_between_pixels;
                int pixel_tissue_depth = (2*(y_max) + 2 )/separation_between_pixels;

                // PRINT_2_VARIABLES(pixel_tissue_width,pixel_tissue_depth);

                unsigned **pixel_grid = new unsigned*[pixel_tissue_depth];
                for(int i = 0; i < pixel_tissue_depth; ++i) 
                {
                    pixel_grid[i] = new unsigned[pixel_tissue_width];
                }

                // int number_of_pixels_with_cell = 0;
                // Calculate the area of the tissue
                for(int pixel_i = 0; pixel_i<pixel_tissue_width; pixel_i++)
                {
                    for(int pixel_j = 0; pixel_j<pixel_tissue_depth; pixel_j++)
                    {
                        pixel_grid[pixel_i][pixel_j] = 0;

                        double pixel_x_coordinate = -1.0 - x_max + pixel_i*separation_between_pixels;
                        double pixel_y_coordinate = -1.0 - y_max + pixel_j*separation_between_pixels;
                        // Iterate over cell population
                        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                            cell_iter != rCellPopulation.End();
                            ++cell_iter)
                        {
                            c_vector<double, DIM>  cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                            if( pow(cell_location[0] - pixel_x_coordinate,2) + pow(cell_location[1] - pixel_y_coordinate,2) < pow(pixel_radial_reach,2) )
                            {
                                pixel_grid[pixel_i][pixel_j] = 1;
                                // pixel_area = pixel_area + area_of_pixel;
                                break;
                            }
                        }

                    }

                }


                // Iterate over cell population
                for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                    cell_iter != rCellPopulation.End();
                    ++cell_iter)
                {
                    cell_iter->GetCellData()->SetItem("is_boundary", 0.0);
                    number_cells++;

                    for(int pixel_i = 0; pixel_i<pixel_tissue_width; pixel_i++)
                    {
                        bool break_out_loop = false;

                        for(int pixel_j = 0; pixel_j<pixel_tissue_depth; pixel_j++)
                        {
                            double pixel_x_coordinate = -1.0 - x_max + pixel_i*separation_between_pixels;
                            double pixel_y_coordinate = -1.0 - y_max + pixel_j*separation_between_pixels;
                                    
                            c_vector<double, DIM>  cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                            if(pixel_grid[pixel_i][pixel_j] == 0)
                            {
                                if( pow(cell_location[0] - pixel_x_coordinate,2) + pow(cell_location[1] - pixel_y_coordinate,2) < pow(pixel_radial_reach+2*separation_between_pixels,2) )
                                {
                                    unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                                    boundary_cells.push_back(real_node_index);
                                    x_boundary.push_back(cell_location[0]);
                                    y_boundary.push_back(cell_location[1]);

                                    cell_iter->GetCellData()->SetItem("is_boundary", 1.0);


                                    break_out_loop = true;
                                    break;
                                }
                            }

                            if(break_out_loop)
                            {
                                break;
                            }

                        }

                        if(break_out_loop)
                        {
                            break;
                        }

                    }
                }

                for(int i = 0; i < pixel_tissue_depth; ++i) 
                {
                    delete [] pixel_grid[i];
                }
                delete [] pixel_grid;
            }




            // Now we have the set of boundary nodes, find the circularity
            double x_centre = 0.0;
            double y_centre = 0.0;
            unsigned number_of_boundary_cells = boundary_cells.size();
            
            for(unsigned i=0; i<number_of_boundary_cells; i++)
            {
                x_centre = x_centre + x_boundary[i];
                y_centre = y_centre + y_boundary[i];
            }
            x_centre = x_centre/number_of_boundary_cells;
            y_centre = y_centre/number_of_boundary_cells;

			double *x_boundary_c = new double[number_of_boundary_cells];
			double *y_boundary_c = new double[number_of_boundary_cells];
            double *theta_coordinates = new double[number_of_boundary_cells];
            unsigned *ordered_set = new unsigned[number_of_boundary_cells];


            for(unsigned i=0; i<number_of_boundary_cells; i++)
            {
                x_boundary_c[i] = x_boundary[i] - x_centre;
                y_boundary_c[i] = y_boundary[i] - y_centre;
            }


            for(unsigned i=0; i<number_of_boundary_cells; i++)
            {
                double length_to_i = sqrt(x_boundary_c[i]*x_boundary_c[i] + y_boundary_c[i]*y_boundary_c[i]);

                if(x_boundary_c[i] >= 0 && y_boundary_c[i] > 0)
                {
                    theta_coordinates[i] = acos(x_boundary_c[i]/length_to_i);
                }
                else if(x_boundary_c[i] < 0 && y_boundary_c[i] > 0)
                {
                    theta_coordinates[i] = 3.141592653589793 - asin(y_boundary_c[i]/length_to_i);
                }
                else if(x_boundary_c[i] < 0 && y_boundary_c[i] < 0)
                {
                    theta_coordinates[i] = atan(y_boundary_c[i]/x_boundary_c[i]) + 3.141592653589793;
                }
                else if(x_boundary_c[i] >= 0 && y_boundary_c[i] < 0)
                {
                    theta_coordinates[i] = 2*3.141592653589793 - acos(x_boundary_c[i]/length_to_i);
                }

            }

            bool *flag_is_in = new bool[number_of_boundary_cells];
            for(unsigned i=0; i<number_of_boundary_cells; i++)
            {
                flag_is_in[i] = false;
            }

            unsigned count = 0;
            
            for(double theta_j=0.0; theta_j<6.29; theta_j+=0.01)
            {
                for(unsigned i=0; i<number_of_boundary_cells; i++)
                {
                    if(flag_is_in[i]==false)
                    {
                        if(theta_j > theta_coordinates[i])
                        {
                            ordered_set[count] = i;
                            // PRINT_4_VARIABLES(x_boundary_c[i],y_boundary_c[i],i,theta_coordinates[i]);
                            flag_is_in[i] = true;
                            count++;
                        }
                    }
                }
            }

            // TRACE("pop");
            for(unsigned i=0; i<number_of_boundary_cells; i++)
            {
                unsigned ii, jj;
                double xii, xjj, yii, yjj;
                if(i == number_of_boundary_cells-1)
                {
                    ii=ordered_set[i];
                    jj=ordered_set[0];
                    
                }
                else
                {
                    ii=ordered_set[i];
                    jj=ordered_set[i+1];
                }
                // PRINT_2_VARIABLES(ii,jj);

                xii = x_boundary_c[ii];
                yii = y_boundary_c[ii];

                xjj = x_boundary_c[jj];
                yjj = y_boundary_c[jj];

                *locationFile_2 << xii << "\t";
                *locationFile_2 << yii << "\t";
                *locationFile_2 << xjj << "\t";
                *locationFile_2 << yjj << "\t";

                tissue_area = tissue_area + 0.5*(xii*yjj-xjj*yii);
                tissue_perimiter = tissue_perimiter + sqrt(pow(xii-xjj,2) + pow(yii-yjj,2));

            }
            delete[] x_boundary_c;
			delete[] y_boundary_c;
            delete[] theta_coordinates;
            delete[] ordered_set;
            delete[] flag_is_in;

        }

        OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
        out_stream locationFile = output_file_handler.OpenOutputFile("CicularityCalc.dat", std::ios::app);
        // SimulationTime* p_time = SimulationTime::Instance();
        *locationFile << p_time->GetTime() << "\t";
        *locationFile << tissue_area << "\t";
        *locationFile << tissue_perimiter << "\t";
        *locationFile << 4*3.141592653589793*tissue_area/(tissue_perimiter*tissue_perimiter) << "\t";
        *locationFile << number_cells << "\t";
        *locationFile << "\n";
        locationFile->close();


        *locationFile_2 << "\n";
        locationFile_2->close();

    }

}



template<unsigned DIM>
void CicularityCalcModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateAtEndOfOutputTimeStep(rCellPopulation);


}


template<unsigned DIM>
void CicularityCalcModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class CicularityCalcModifier<1>;
template class CicularityCalcModifier<2>;
template class CicularityCalcModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CicularityCalcModifier)

