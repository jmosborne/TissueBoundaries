#ifndef TESTGROWINGMONOLAYER_HPP_
#define TESTGROWINGMONOLAYER_HPP_


#include <cxxtest/TestSuite.h> 

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp" 

#include "SmartPointers.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"

#include "HoneycombMeshGenerator.hpp"
#include "MutableMesh.hpp"
#include "NodesOnlyMesh.hpp"
#include "HoneycombVertexMeshGenerator.hpp"

#include "CellsGenerator.hpp"

#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "Cylindrical2dNodesOnlyMesh.hpp"

#include "ContactInhibitionCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "CellProliferativeTypesCountWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellAncestorWriter.hpp"

#include "OffLatticeSimulation.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"

#include "VolumeTrackingModifier.hpp"
#include "SimpleTargetAreaModifier.hpp"
// #include "VertexBoundaryRefinementModifier.hpp"
#include "CurvedVertexEdgesModifierSimplified.hpp"
#include "SmoothVertexEdgesModifierSimplified.hpp"
#include "JaggedVertexEdgesModifierSimplified.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Warnings.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"

// #include "CicularityCalcModifier.hpp"
// #include "BoundaryCellWriter.hpp"

#include "GhostNodeCollisionRemovalModifier.hpp"

#include "Debug.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"
#include "CellLabel.hpp"
#include "CellLabelWriter.hpp"

#include "CellMutationStatesWriter.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "BernoulliTrialWithContactInhibitionCellCycleModel.hpp"
#include "ModifiedBernoulliTrialWithContactInhibitionCellCycleModel.hpp"
#include "TwoPopulationLineModifyer.hpp"

#include "NodeBoundryWriterModifier.hpp"


/*
 *  This is where you can set parameters to be used in all the simulations.
 */

static const double M_END_TIME = 1; //100
static const double M_DT_TIME = 0.001;
static const double M_SAMPLE_TIME = 20;

static const double M_CONTACT_INHIBITION_LEVEL = 0.8;
static const double M_STEM_CELL_DIVISION_PROBABILITY = 0.1;
static const double M_STEM_CELL_MINIMUM_DIVISION_AGE = 1.0;

static const double M_INITIAL_WIDTH = 12;
static const double M_INITIAL_LENGTH = 14;

static const double M_DOMAIN_WIDTH = M_INITIAL_WIDTH;
static const double M_DOMAIN_HEIGHT = M_INITIAL_LENGTH*sqrt(3.0/4.0);


static const std::string M_HEAD_FOLDER = "Competing_Populations_Cylindrical";

class TestCompetingPopulations_Loop : public AbstractCellBasedWithTimingsTestSuite
{
private:


    /**
    * Helper method. Smooth out edges of a vertex mesh - move this to .hpp, .cpp
    * 
    * @param rCellPopulation a cell population
    */
    void SmoothVertexMeshEdges(AbstractCellPopulation<2>& rCellPopulation)
    {
        MutableVertexMesh<2, 2>& r_mesh = static_cast<VertexBasedCellPopulation<2>* >(&rCellPopulation)->rGetMesh();

        for (VertexMesh<2,2>::NodeIterator node_iter = r_mesh.GetNodeIteratorBegin();
            node_iter != r_mesh.GetNodeIteratorEnd();
            ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();
            std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
            if (containing_element_indices.size() == 1)
            {
                // Get this element
                unsigned elem_index = (*containing_element_indices.begin());

                VertexElement<2,2>* p_element = r_mesh.GetElement(elem_index);

                // Remove node from this element and delete the node
                p_element->DeleteNode(p_element->GetNodeLocalIndex(node_index));
                r_mesh.DeleteNodePriorToReMesh(node_index);
            }
        }
            r_mesh.ReMesh();
    }

    /**
    * Helper method. Iterate over all cells and 'round out' the domain by
    * killing those cells whose centres are located outside a given region.
    * 
    * @param rCellPopulation a cell population
    */
    void RoundOutCellPopulation(AbstractCellPopulation<2>& rCellPopulation, double x1, double y1, double x2, double y2)
    {
        

        if (bool(dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>*>(&rCellPopulation)))
        {
            std::set<unsigned> location_indices;
            std::set<unsigned> ghost_node_indices = dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>*>(&rCellPopulation)->GetGhostNodeIndices();

            for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
                cell_iter != rCellPopulation.rGetCells().end();)
            {
                // Get the coordinates of this cell centre
                c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                double x = centre_of_cell[0];
                double y = centre_of_cell[1];
                unsigned location_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
 
                // if ((fabs(y-x)>M_DOMAINWIDTH) && (x<M_DOMAIN_X_MIN) && (x>M_DOMAIN_X_MAX) && (y<M_DOMAIN_Y_MIN) && (y>M_DOMAIN_Y_MAX))
                // if ( pow(x-x1,2) + pow(y-y1,2) < pow(2.0,2) || pow(x-x2,2) + pow(y-y2,2) < pow(2.0,2))
                // if (  x <= x1 || x >= x2 )
                if (  y <= y1 || y >= y2 )
                {   
                    ++cell_iter;
                }
                else
                {
                    // Delete cell and store it as a ghost node
                    rCellPopulation.RemoveCellUsingLocationIndex(location_index, (*cell_iter));
                    // Update vector of cells
                    cell_iter = rCellPopulation.rGetCells().erase(cell_iter);
 
                    // Change to chost node            
                    ghost_node_indices.insert(location_index);
                    
                }
            }
            dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>* >(&rCellPopulation)->SetGhostNodes(ghost_node_indices);
            dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>* >(&rCellPopulation)->RemoveDeadCells();
            dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>* >(&rCellPopulation)->Update();
        }
        else
        {

            for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
                    cell_iter != rCellPopulation.End();
                    ++cell_iter)
            {
                // Get the coordinates of this cell centre
                c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                double x = centre_of_cell[0];
                double y = centre_of_cell[1];


                // Use circles
                // if(pow(x-x1,2) + pow(y-y1,2) < pow(2.0,2))
                // {
                //     // label as 1
                // }
                // else if(pow(x-x2,2) + pow(y-y2,2) < pow(2.0,2))
                // {
                //     // label as 2
                // }

                // Use Plane Boundary
                // if( x <= x1 )
                // {
                //     // label as 1
                // }
                // else if( x >= x2 )
                // {
                //     // label as 2
                // }
                if( y <= y1 )
                {
                    // label as 1
                }
                else if( y >= y2 )
                {
                    // label as 2
                }
                else
                {
                    cell_iter->Kill();
                }
            }
            
            /* Need to remove the deleted cells and call update note this is usually
            * performed in the Solve() method of the simulation class.
            */
            if (bool(dynamic_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation)))
            {
                dynamic_cast<NodeBasedCellPopulation<2>* >(&rCellPopulation)->RemoveDeadCells();
                dynamic_cast<NodeBasedCellPopulation<2>* >(&rCellPopulation)->Update();
            }
            else if (bool(dynamic_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation)))
            {
                dynamic_cast<MeshBasedCellPopulation<2>* >(&rCellPopulation)->RemoveDeadCells();
                dynamic_cast<MeshBasedCellPopulation<2>* >(&rCellPopulation)->Update();
            }
            else if (bool(dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation)))
            {
                dynamic_cast<VertexBasedCellPopulation<2>* >(&rCellPopulation)->RemoveDeadCells();
                dynamic_cast<VertexBasedCellPopulation<2>* >(&rCellPopulation)->Update();
            }

        }
    }
    
    void LabelCells(AbstractCellPopulation<2>& rCellPopulation , boost::shared_ptr<AbstractCellProperty> pLabel, double x1, double y1, double x2, double y2)
    {       

        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
        {
            // Get the coordinates of this cell centre
            c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            double x = centre_of_cell[0];
            double y = centre_of_cell[1];

            // Use circles
            // if(pow(x-x1,2) + pow(y-y1,2) < pow(2.0,2))
            // {
            //     (*cell_iter)->AddCellProperty(pLabel);
            //     // label as 1
            // }
            // else if(pow(x-x2,2) + pow(y-y2,2) < pow(2.0,2))
            // {
            //     // label as 2
            // }

            // Use Plane Boundary
            // if( x <= x1 )
            // {
            //     (*cell_iter)->AddCellProperty(pLabel);
            //     // label as 1
            // }
            // else if( x >= x2 )
            // {
            //     // label as 2
            // }
            if( y <= y1 )
            {
                (*cell_iter)->AddCellProperty(pLabel);
                // label as 1
            }
            else if( y >= y2 )
            {
                // label as 2
            }
        }

    }

    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells, double equilibriumVolume, double quiescentVolumeFraction, 
            double stemCellDivisionProbability, double stemCellMinimumDivisionAge)
    {

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_transit_cell_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_stem_cell_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        for (unsigned i=0; i<num_cells; i++)
        {
            // BernoulliTrialWithContactInhibitionCellCycleModel* p_model = new BernoulliTrialWithContactInhibitionCellCycleModel();
            ModifiedBernoulliTrialWithContactInhibitionCellCycleModel* p_model = new ModifiedBernoulliTrialWithContactInhibitionCellCycleModel();
            
            p_model->SetDimension(2);
            p_model->SetEquilibriumVolume(equilibriumVolume);
            p_model->SetQuiescentVolumeFraction(quiescentVolumeFraction);
            p_model->SetStemCellDivisionProbability(stemCellDivisionProbability);
            p_model->SetStemCellMinimumDivisionAge(stemCellMinimumDivisionAge);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_cell_type);
            double ave_stem_cell_cycle_duration = p_model->GetAverageStemCellCycleTime();
            double birth_time = - RandomNumberGenerator::Instance()->ranf() * ave_stem_cell_cycle_duration;
            p_cell->SetBirthTime(birth_time);

            // Set Target Area so dont need to use a growth model in vertex simulations
            p_cell->GetCellData()->SetItem("target area", sqrt(3.0)/2.0);
            // p_cell->GetCellData()->SetItem("volume", M_PI*0.25);

            rCells.push_back(p_cell);
        }
    }

public:

    void TestCompetingCellPopulations()

    {

        //Command line argument stuff: Get the seed parameters  
        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-run_index"));
        unsigned start_index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-run_index");
        
        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-num_runs"));
        unsigned num_runs = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-num_runs");
        
        // unsigned start_index = 0;
        // unsigned num_runs = 10;

        // Loop over the random seed.
        for(unsigned sim_index=start_index; sim_index < start_index + num_runs; sim_index++)
        {

            std::cout << " Run number " << sim_index << "... \n" << std::flush;

            // Reseed the random number generator
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            
            std::stringstream out;
            out << sim_index;

            //////


            /* 
            * == Mesh Ghosts ==
            */
            {
                std::cout << "Mesh Ghosts \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Mesh/Ghosts/Run_" +  out.str();

                // Create a simple 2D MutableVertexMesh
                unsigned thickness_of_ghost_layer = 4;
                CylindricalHoneycombMeshGenerator generator(M_INITIAL_WIDTH, M_INITIAL_LENGTH, thickness_of_ghost_layer);
                Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

                std::vector<unsigned> location_indices = generator.GetCellLocationIndices();


                std::vector<CellPtr> cells;
                GenerateCells(location_indices.size(), cells, sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: M_PI*0.25 as r=0.5


                MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

                // Set population to output all data to results files
                cell_population.AddPopulationWriter<VoronoiDataWriter>();
                cell_population.SetWriteVtkAsPoints(true);

                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellMutationStatesWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();


                double x1 = 0; double y1 = (0.5/6.0)*M_DOMAIN_HEIGHT;
                double x2 = 0; double y2 = (5.0/6.0)*M_DOMAIN_HEIGHT;

                RoundOutCellPopulation(cell_population, x1,  y1,  x2,  y2);

                // Set up cell-based simulation and output directory
                OffLatticeSimulation<2> simulator(cell_population);

                // Set up force law and pass it to the simulation
                MAKE_PTR(DifferentialAdhesionGeneralisedLinearSpringForce<2>, p_differential_adhesion_force);
                p_differential_adhesion_force->SetMeinekeSpringStiffness(50.0);
                p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
                p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.1);
                p_differential_adhesion_force->SetCutOffLength(1.5);
                simulator.AddForce(p_differential_adhesion_force);


                simulator.SetOutputDirectory(output_directory);

                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Create some boundary conditions and pass them to the simulation
                c_vector<double,2> point = zero_vector<double>(2);
                c_vector<double,2> normal = zero_vector<double>(2);
                normal(1) = -1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal)); // y>0
                simulator.AddCellPopulationBoundaryCondition(p_bc1);
                point(1) = M_DOMAIN_HEIGHT;
                normal(1) = 1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                simulator.AddCellPopulationBoundaryCondition(p_bc2);
                // normal(0) = 1.0; normal(1) = 0.0;
                // point(0) = M_DOMAIN_WIDTH; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal)); // y>0
                // simulator.AddCellPopulationBoundaryCondition(p_bc3);
                // normal(0) = -1.0; normal(1) = 0.0;
                // point(0) = 0.0; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                // simulator.AddCellPopulationBoundaryCondition(p_bc4);

                // Now label some cells

                boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
                LabelCells(cell_population, p_state,  x1,  y1,  x2,  y2);

                MAKE_PTR(GhostNodeCollisionRemovalModifier<2>, GNR_modifier);
                simulator.AddSimulationModifier(GNR_modifier);   

                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME);

                MAKE_PTR(TwoPopulationLineModifyer<2>, p_line_modifier);
                p_line_modifier->SetOutputDirectory(output_directory);
                p_line_modifier->SetCutoff(1.5);
                simulator.AddSimulationModifier(p_line_modifier);
                
                simulator.Solve();
                
                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

            }

            
            /* 
            * == Mesh Finite ==
            */
            {
                std::cout << "Mesh Finite \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Mesh/Finite/Run_" +  out.str();

                // Create a simple 2D MutableVertexMesh
                CylindricalHoneycombMeshGenerator generator(M_INITIAL_WIDTH, M_INITIAL_LENGTH, 0);
                Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumNodes(), cells, sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: M_PI*0.25 as r=0.5

                MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

                // Set population to output all data to results files
                cell_population.AddPopulationWriter<VoronoiDataWriter>();
                cell_population.SetWriteVtkAsPoints(true);
                cell_population.SetBoundVoronoiTessellation(true);


                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellMutationStatesWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();


                double x1 = 0; double y1 = (0.5/6.0)*M_DOMAIN_HEIGHT;
                double x2 = 0; double y2 = (5.0/6.0)*M_DOMAIN_HEIGHT;

                RoundOutCellPopulation(cell_population, x1,  y1,  x2,  y2);

                // Set up cell-based simulation and output directory
                OffLatticeSimulation<2> simulator(cell_population);

                // Set up force law and pass it to the simulation
                MAKE_PTR(DifferentialAdhesionGeneralisedLinearSpringForce<2>, p_differential_adhesion_force);
                p_differential_adhesion_force->SetMeinekeSpringStiffness(50.0);
                p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
                p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.1);
                p_differential_adhesion_force->SetCutOffLength(1.5);
                simulator.AddForce(p_differential_adhesion_force);

                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                simulator.SetOutputDirectory(output_directory);


                // Create some boundary conditions and pass them to the simulation
                c_vector<double,2> point = zero_vector<double>(2);
                c_vector<double,2> normal = zero_vector<double>(2);
                normal(1) = -1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal)); // y>0
                simulator.AddCellPopulationBoundaryCondition(p_bc1);
                point(1) = M_DOMAIN_HEIGHT;
                normal(1) = 1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                simulator.AddCellPopulationBoundaryCondition(p_bc2);
                // normal(0) = 1.0; normal(1) = 0.0;
                // point(0) = M_DOMAIN_WIDTH; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal)); // y>0
                // simulator.AddCellPopulationBoundaryCondition(p_bc3);
                // normal(0) = -1.0; normal(1) = 0.0;
                // point(0) = 0.0; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                // simulator.AddCellPopulationBoundaryCondition(p_bc4);

                // Now label some cells

                boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
                LabelCells(cell_population, p_state,  x1,  y1,  x2,  y2);

                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME);

                MAKE_PTR(TwoPopulationLineModifyer<2>, p_line_modifier);
                p_line_modifier->SetOutputDirectory(output_directory);
                p_line_modifier->SetCutoff(1.5);
                simulator.AddSimulationModifier(p_line_modifier);
                
                simulator.Solve();

                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

            }


            /* 
            * == Mesh Infinite ==
            */
            {
                std::cout << "Mesh Infinite \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Mesh/Infinite/Run_" +  out.str();

                // Create a simple 2D MutableVertexMesh
                CylindricalHoneycombMeshGenerator generator(M_INITIAL_WIDTH, M_INITIAL_LENGTH, 0);
                Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

                std::vector<CellPtr> cells;
                
                GenerateCells(p_mesh->GetNumNodes(), cells, sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: M_PI*0.25 as r=0.5

                MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

                // Set population to output all data to results files
                cell_population.AddPopulationWriter<VoronoiDataWriter>();
                cell_population.SetWriteVtkAsPoints(true);


                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellMutationStatesWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();

                double x1 = 0; double y1 = (0.5/6.0)*M_DOMAIN_HEIGHT;
                double x2 = 0; double y2 = (5.0/6.0)*M_DOMAIN_HEIGHT;

                RoundOutCellPopulation(cell_population, x1,  y1,  x2,  y2);

                // Set up cell-based simulation and output directory
                OffLatticeSimulation<2> simulator(cell_population);

                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Set up force law and pass it to the simulation
                MAKE_PTR(DifferentialAdhesionGeneralisedLinearSpringForce<2>, p_differential_adhesion_force);
                p_differential_adhesion_force->SetMeinekeSpringStiffness(50.0);
                p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
                p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.1);
                p_differential_adhesion_force->SetCutOffLength(1.5);
                simulator.AddForce(p_differential_adhesion_force);

                simulator.SetOutputDirectory(output_directory);

                // Now label some cells

                boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
                LabelCells(cell_population, p_state,  x1,  y1,  x2,  y2);

                // Create some boundary conditions and pass them to the simulation
                c_vector<double,2> point = zero_vector<double>(2);
                c_vector<double,2> normal = zero_vector<double>(2);
                normal(1) = -1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal)); // y>0
                simulator.AddCellPopulationBoundaryCondition(p_bc1);
                point(1) = M_DOMAIN_HEIGHT;
                normal(1) = 1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                simulator.AddCellPopulationBoundaryCondition(p_bc2);
                // normal(0) = 1.0; normal(1) = 0.0;
                // point(0) = M_DOMAIN_WIDTH; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal)); // y>0
                // simulator.AddCellPopulationBoundaryCondition(p_bc3);
                // normal(0) = -1.0; normal(1) = 0.0;
                // point(0) = 0.0; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                // simulator.AddCellPopulationBoundaryCondition(p_bc4);

                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME);

                MAKE_PTR(TwoPopulationLineModifyer<2>, p_line_modifier);
                p_line_modifier->SetOutputDirectory(output_directory);
                p_line_modifier->SetCutoff(1.5);
                simulator.AddSimulationModifier(p_line_modifier);
                
                simulator.Solve();

                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

            }
            
            /* 
            * == Node Default ==
            */
            {
                std::cout << "Node Default \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Node/Default/Run_" +  out.str();

                // Create a simple 2D MutableVertexMesh
                HoneycombMeshGenerator generator(M_INITIAL_WIDTH, M_INITIAL_LENGTH);
                TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

                Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(M_INITIAL_WIDTH);
                p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh,2.0);

                // Slows things down but can use a larger timestep and diffusion forces

                // NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
                // p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0);

                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumNodes(), cells, M_PI*0.25,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: M_PI*0.25 as r=0.5

                NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);

                // Set population to output all data to results files
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellMutationStatesWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();

                double x1 = 0; double y1 = (0.5/6.0)*M_DOMAIN_HEIGHT;
                double x2 = 0; double y2 = (5.0/6.0)*M_DOMAIN_HEIGHT;

                RoundOutCellPopulation(cell_population, x1,  y1,  x2,  y2);

                // Set up cell-based simulation and output directory
                OffLatticeSimulation<2> simulator(cell_population);

                // Set up force law and pass it to the simulation
                MAKE_PTR(DifferentialAdhesionGeneralisedLinearSpringForce<2>, p_differential_adhesion_force);
                p_differential_adhesion_force->SetMeinekeSpringStiffness(50.0);
                p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
                p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.1);
                p_differential_adhesion_force->SetCutOffLength(1.5);
                simulator.AddForce(p_differential_adhesion_force);

                simulator.SetOutputDirectory(output_directory);

                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Create some boundary conditions and pass them to the simulation
                c_vector<double,2> point = zero_vector<double>(2);
                c_vector<double,2> normal = zero_vector<double>(2);
                normal(1) = -1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal)); // y>0
                simulator.AddCellPopulationBoundaryCondition(p_bc1);
                point(1) = M_DOMAIN_HEIGHT;
                normal(1) = 1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                simulator.AddCellPopulationBoundaryCondition(p_bc2);
                // normal(0) = 1.0; normal(1) = 0.0;
                // point(0) = M_DOMAIN_WIDTH; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal)); // y>0
                // simulator.AddCellPopulationBoundaryCondition(p_bc3);
                // normal(0) = -1.0; normal(1) = 0.0;
                // point(0) = 0.0; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                // simulator.AddCellPopulationBoundaryCondition(p_bc4);

                // Now label some cells
                boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
                LabelCells(cell_population, p_state,  x1,  y1,  x2,  y2);

                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME);

                MAKE_PTR(TwoPopulationLineModifyer<2>, p_line_modifier);
                p_line_modifier->SetOutputDirectory(output_directory);
                p_line_modifier->SetCutoff(1.5);
                simulator.AddSimulationModifier(p_line_modifier);
                
                simulator.Solve();

                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

                // Clear memory
                delete p_mesh;

            }

            /* 
            * == Node Small ==
            */
            {
                std::cout << "Node Small \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Node/Small/Run_" +  out.str();
                // Create a simple 2D MutableVertexMesh
                HoneycombMeshGenerator generator(M_INITIAL_WIDTH, M_INITIAL_LENGTH);
                TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

                Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(M_INITIAL_WIDTH);
                p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh,2.0);

                std::vector<CellPtr> cells;
                
                GenerateCells(p_mesh->GetNumNodes(), cells, M_PI*0.25,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: M_PI*0.25 as r=0.5

                NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);

                // Set population to output all data to results files
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellMutationStatesWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();


                double x1 = 0; double y1 = (0.5/6.0)*M_DOMAIN_HEIGHT;
                double x2 = 0; double y2 = (5.0/6.0)*M_DOMAIN_HEIGHT;

                RoundOutCellPopulation(cell_population, x1,  y1,  x2,  y2);

                // Set up cell-based simulation and output directory
                OffLatticeSimulation<2> simulator(cell_population);

                // Set up force law and pass it to the simulation
                MAKE_PTR(DifferentialAdhesionGeneralisedLinearSpringForce<2>, p_differential_adhesion_force);
                p_differential_adhesion_force->SetMeinekeSpringStiffness(50.0);
                p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
                p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.1);
                p_differential_adhesion_force->SetCutOffLength(1.0);
                simulator.AddForce(p_differential_adhesion_force);

                simulator.SetOutputDirectory(output_directory);

                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Create some boundary conditions and pass them to the simulation
                c_vector<double,2> point = zero_vector<double>(2);
                c_vector<double,2> normal = zero_vector<double>(2);
                normal(1) = -1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal)); // y>0
                simulator.AddCellPopulationBoundaryCondition(p_bc1);
                point(1) = M_DOMAIN_HEIGHT;
                normal(1) = 1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                simulator.AddCellPopulationBoundaryCondition(p_bc2);
                // normal(0) = 1.0; normal(1) = 0.0;
                // point(0) = M_DOMAIN_WIDTH; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal)); // y>0
                // simulator.AddCellPopulationBoundaryCondition(p_bc3);
                // normal(0) = -1.0; normal(1) = 0.0;
                // point(0) = 0.0; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                // simulator.AddCellPopulationBoundaryCondition(p_bc4);

                // Now label some cells

                boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
                LabelCells(cell_population, p_state,  x1,  y1,  x2,  y2);

                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME);

                MAKE_PTR(TwoPopulationLineModifyer<2>, p_line_modifier);
                p_line_modifier->SetOutputDirectory(output_directory);
                p_line_modifier->SetCutoff(1.5);
                simulator.AddSimulationModifier(p_line_modifier);
                
                simulator.Solve();

                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

                // Clear memory
                delete p_mesh;

            }

            /* 
            * == Node Large ==
            */
            {
                std::cout << "Node Large \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Node/Large/Run_" +  out.str();
                // Create a simple 2D MutableVertexMesh
                HoneycombMeshGenerator generator(M_INITIAL_WIDTH, M_INITIAL_LENGTH);
                TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

                Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(M_INITIAL_WIDTH);
                p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh,2.0);

                std::vector<CellPtr> cells;
                
                GenerateCells(p_mesh->GetNumNodes(), cells, M_PI*0.25,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: M_PI*0.25 as r=0.5

                NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);

                // Set population to output all data to results files
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellMutationStatesWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();


                double x1 = 0; double y1 = (0.5/6.0)*M_DOMAIN_HEIGHT;
                double x2 = 0; double y2 = (5.0/6.0)*M_DOMAIN_HEIGHT;

                RoundOutCellPopulation(cell_population, x1,  y1,  x2,  y2);

                // Set up cell-based simulation and output directory
                OffLatticeSimulation<2> simulator(cell_population);

                // Set up force law and pass it to the simulation
                MAKE_PTR(DifferentialAdhesionGeneralisedLinearSpringForce<2>, p_differential_adhesion_force);
                p_differential_adhesion_force->SetMeinekeSpringStiffness(50.0);
                p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
                p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.1);
                p_differential_adhesion_force->SetCutOffLength(2.0);
                simulator.AddForce(p_differential_adhesion_force);

                simulator.SetOutputDirectory(output_directory);

                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Create some boundary conditions and pass them to the simulation
                c_vector<double,2> point = zero_vector<double>(2);
                c_vector<double,2> normal = zero_vector<double>(2);
                normal(1) = -1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal)); // y>0
                simulator.AddCellPopulationBoundaryCondition(p_bc1);
                point(1) = M_DOMAIN_HEIGHT;
                normal(1) = 1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                simulator.AddCellPopulationBoundaryCondition(p_bc2);
                // normal(0) = 1.0; normal(1) = 0.0;
                // point(0) = M_DOMAIN_WIDTH; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal)); // y>0
                // simulator.AddCellPopulationBoundaryCondition(p_bc3);
                // normal(0) = -1.0; normal(1) = 0.0;
                // point(0) = 0.0; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                // simulator.AddCellPopulationBoundaryCondition(p_bc4);

                // Now label some cells

                boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
                LabelCells(cell_population, p_state,  x1,  y1,  x2,  y2);

                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME);

                MAKE_PTR(TwoPopulationLineModifyer<2>, p_line_modifier);
                p_line_modifier->SetOutputDirectory(output_directory);
                p_line_modifier->SetCutoff(1.5);
                simulator.AddSimulationModifier(p_line_modifier);
                
                simulator.Solve();

                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

                // Clear memory
                delete p_mesh;

            }

            /* 
            * == Vertex Jagged ==
            */
            {
                std::cout << "Vertex Jagged \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Vertex/Jagged/Run_" +  out.str();
                // Create a simple 2D MutableVertexMesh
                CylindricalHoneycombVertexMeshGenerator generator(M_INITIAL_WIDTH, M_INITIAL_LENGTH, false);
                Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

                p_mesh->Translate( -0.5, -0.575);

                p_mesh->SetCellRearrangementThreshold(0.15);

                // Slows things down but can use a larger timestep and diffusion forces
                //p_mesh->SetCheckForInternalIntersections(true);

                // Set up cells, one for each VertexElement
                std::vector<CellPtr> cells;
                
                GenerateCells(p_mesh->GetNumElements(), cells, sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: M_PI*0.25 as r=0.5


                // Create cell population
                VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

                // Set population to output all data to results files
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellMutationStatesWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();


                double x1 = 0; double y1 = (0.5/6.0)*M_DOMAIN_HEIGHT;
                double x2 = 0; double y2 = (5.0/6.0)*M_DOMAIN_HEIGHT;

                RoundOutCellPopulation(cell_population, x1,  y1,  x2,  y2);

                // Set up cell-based simulation and output directory
                OffLatticeSimulation<2> simulator(cell_population);

                // Set up force law and pass it to the simulation
                MAKE_PTR(NagaiHondaDifferentialAdhesionForce<2>, p_force);


                double alpha = 50.0;
                double beta =  1.0;
                double gamma_int = 1.0;//1/sqrt(3);
                double gamma_free = 2*gamma_int;
                p_force->SetNagaiHondaDeformationEnergyParameter(alpha);
                p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(beta);

                p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(gamma_int);
                p_force->SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(2*gamma_int);
                p_force->SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(gamma_int);

                p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(gamma_free);
                p_force->SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(gamma_free);
                simulator.AddForce(p_force);

                simulator.SetOutputDirectory(output_directory);

                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                MAKE_PTR(JaggedVertexEdgesModifierSimplified<2>, jagged_edge_modifier);
                simulator.AddSimulationModifier(jagged_edge_modifier);

                // Create some boundary conditions and pass them to the simulation
                c_vector<double,2> point = zero_vector<double>(2);
                c_vector<double,2> normal = zero_vector<double>(2);
                point(1) = -0.575;
                normal(1) = -1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal)); // y>0
                simulator.AddCellPopulationBoundaryCondition(p_bc1);
                point(1) = M_DOMAIN_HEIGHT -0.575/2.0;
                normal(1) = 1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                simulator.AddCellPopulationBoundaryCondition(p_bc2);
                // normal(0) = 1.0; normal(1) = 0.0;
                // point(0) = M_DOMAIN_WIDTH; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal)); // y>0
                // simulator.AddCellPopulationBoundaryCondition(p_bc3);
                // normal(0) = -1.0; normal(1) = 0.0;
                // point(0) = 0.0; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                // simulator.AddCellPopulationBoundaryCondition(p_bc4);

                // Now label some cells

                boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
                LabelCells(cell_population, p_state,  x1,  y1,  x2,  y2);

                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME);

                MAKE_PTR(NodeBoundryWriterModifier<2>, p_boundary_modifier);
                p_boundary_modifier->SetOutputDirectory(output_directory + "/results_from_time_0");
                simulator.AddSimulationModifier(p_boundary_modifier);

                MAKE_PTR(TwoPopulationLineModifyer<2>, p_line_modifier);
                p_line_modifier->SetOutputDirectory(output_directory);
                simulator.AddSimulationModifier(p_line_modifier);
                
                simulator.Solve();

                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

                // Clear singletons
                Warnings::Instance()->QuietDestroy();

            }

            /* 
            * == Vertex Smooth ==
            */
            {
                std::cout << "Vertex Smooth \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Vertex/Smooth/Run_" +  out.str();
                // Create a simple 2D MutableVertexMesh
                CylindricalHoneycombVertexMeshGenerator generator(M_INITIAL_WIDTH, M_INITIAL_LENGTH, false);
                Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

                p_mesh->Translate( -0.5, -0.575);

                p_mesh->SetCellRearrangementThreshold(0.1);
                // Slows things down but can use a larger timestep and diffusion forces
                //p_mesh->SetCheckForInternalIntersections(true);

                // Set up cells, one for each VertexElement
                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumElements(), cells, sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: M_PI*0.25 as r=0.5

                // Create cell population
                VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
                SmoothVertexMeshEdges(cell_population);

                // Set population to output all data to results files
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellMutationStatesWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();


                double x1 = 0; double y1 = (0.5/6.0)*M_DOMAIN_HEIGHT;
                double x2 = 0; double y2 = (5.0/6.0)*M_DOMAIN_HEIGHT;

                RoundOutCellPopulation(cell_population, x1,  y1,  x2,  y2);

                // Set up cell-based simulation and output directory
                OffLatticeSimulation<2> simulator(cell_population);

                // Set up force law and pass it to the simulation
                MAKE_PTR(NagaiHondaDifferentialAdhesionForce<2>, p_force);


                double alpha = 50.0;
                double beta =  1.0;
                double gamma_int = 1.0;//1/sqrt(3);
                double gamma_free = 2*gamma_int;
                p_force->SetNagaiHondaDeformationEnergyParameter(alpha);
                p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(beta);

                p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(gamma_int);
                p_force->SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(2*gamma_int);
                p_force->SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(gamma_int);

                p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(gamma_free);
                p_force->SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(gamma_free);
                simulator.AddForce(p_force);

                MAKE_PTR(SmoothVertexEdgesModifierSimplified<2>, smooth_edge_modifier);
                simulator.AddSimulationModifier(smooth_edge_modifier);

                simulator.SetOutputDirectory(output_directory);

                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Create some boundary conditions and pass them to the simulation
                c_vector<double,2> point = zero_vector<double>(2);
                c_vector<double,2> normal = zero_vector<double>(2);
                point(1) = -0.575;
                normal(1) = -1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal)); // y>0
                simulator.AddCellPopulationBoundaryCondition(p_bc1);
                point(1) = M_DOMAIN_HEIGHT -0.575/2.0;
                normal(1) = 1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                simulator.AddCellPopulationBoundaryCondition(p_bc2);
                // normal(0) = 1.0; normal(1) = 0.0;
                // point(0) = M_DOMAIN_WIDTH; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal)); // y>0
                // simulator.AddCellPopulationBoundaryCondition(p_bc3);
                // normal(0) = -1.0; normal(1) = 0.0;
                // point(0) = 0.0; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                // simulator.AddCellPopulationBoundaryCondition(p_bc4);

                // Now label some cells

                boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
                LabelCells(cell_population, p_state,  x1,  y1,  x2,  y2);

                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME);

                MAKE_PTR(NodeBoundryWriterModifier<2>, p_boundary_modifier);
                p_boundary_modifier->SetOutputDirectory(output_directory + "/results_from_time_0");
                simulator.AddSimulationModifier(p_boundary_modifier);

                MAKE_PTR(TwoPopulationLineModifyer<2>, p_line_modifier);
                p_line_modifier->SetOutputDirectory(output_directory);
                simulator.AddSimulationModifier(p_line_modifier);
                
                simulator.Solve();

                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

                // Clear singletons
                Warnings::Instance()->QuietDestroy();

            }

            /* 
            * == Vertex Curved ==
            */
            {
                std::cout << "Vertex Curved \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Vertex/Curved/Run_" +  out.str();
                // Create a simple 2D MutableVertexMesh
                CylindricalHoneycombVertexMeshGenerator generator(M_INITIAL_WIDTH, M_INITIAL_LENGTH, false);
                Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

                p_mesh->Translate( -0.5, -0.575);

                p_mesh->SetCellRearrangementThreshold(0.01);

                // Slows things down but can use a larger timestep and diffusion forces
                //p_mesh->SetCheckForInternalIntersections(true);

                // Set up cells, one for each VertexElement
                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumElements(), cells, sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: M_PI*0.25 as r=0.5

                // Create cell population
                VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

                // Set population to output all data to results files
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellMutationStatesWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();


                double x1 = 0; double y1 = (0.5/6.0)*M_DOMAIN_HEIGHT;
                double x2 = 0; double y2 = (5.0/6.0)*M_DOMAIN_HEIGHT;

                RoundOutCellPopulation(cell_population, x1,  y1,  x2,  y2);

                // Set up cell-based simulation and output directory
                OffLatticeSimulation<2> simulator(cell_population);

                // Set up force law and pass it to the simulation
                MAKE_PTR(NagaiHondaDifferentialAdhesionForce<2>, p_force);


                double alpha = 50.0;
                double beta =  1.0;
                double gamma_int = 1.0;//1/sqrt(3);
                double gamma_free = 2*gamma_int;
                p_force->SetNagaiHondaDeformationEnergyParameter(alpha);
                p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(beta);

                p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(gamma_int);
                p_force->SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(2*gamma_int);
                p_force->SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(gamma_int);

                p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(gamma_free);
                p_force->SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(gamma_free);
                simulator.AddForce(p_force);

                MAKE_PTR(CurvedVertexEdgesModifierSimplified<2>, refinement_modifier);
                simulator.AddSimulationModifier(refinement_modifier);

                simulator.SetOutputDirectory(output_directory);

                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Create some boundary conditions and pass them to the simulation
                c_vector<double,2> point = zero_vector<double>(2);
                c_vector<double,2> normal = zero_vector<double>(2);
                point(1) = -0.575;
                normal(1) = -1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal)); // y>0
                simulator.AddCellPopulationBoundaryCondition(p_bc1);
                point(1) = M_DOMAIN_HEIGHT -0.575/2.0;
                normal(1) = 1.0;
                MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                simulator.AddCellPopulationBoundaryCondition(p_bc2);
                // normal(0) = 1.0; normal(1) = 0.0;
                // point(0) = M_DOMAIN_WIDTH; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal)); // y>0
                // simulator.AddCellPopulationBoundaryCondition(p_bc3);
                // normal(0) = -1.0; normal(1) = 0.0;
                // point(0) = 0.0; point(1) = 0.0;
                // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal)); // y<M_DOMAIN_HEIGHT
                // simulator.AddCellPopulationBoundaryCondition(p_bc4);

                // Now label some cells

                boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
                LabelCells(cell_population, p_state,  x1,  y1,  x2,  y2);

                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME);

                MAKE_PTR(NodeBoundryWriterModifier<2>, p_boundary_modifier);
                p_boundary_modifier->SetOutputDirectory(output_directory + "/results_from_time_0");
                simulator.AddSimulationModifier(p_boundary_modifier);

                MAKE_PTR(TwoPopulationLineModifyer<2>, p_line_modifier);
                p_line_modifier->SetOutputDirectory(output_directory);
                simulator.AddSimulationModifier(p_line_modifier);
                
                simulator.Solve();

                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

                // Clear singletons
                Warnings::Instance()->QuietDestroy();

            }
            
        }

    }






};

#endif /* TestCompetingPopulations_Loop_HPP_ */
