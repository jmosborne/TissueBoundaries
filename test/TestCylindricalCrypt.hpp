#ifndef TESTCYLINDRICALCRYPT_HPP_
#define TESTCYLINDRICALCRYPT_HPP_


#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "SmartPointers.hpp"

#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "Cylindrical2dNodesOnlyMesh.hpp"

#include "CellsGenerator.hpp"

#include "SimpleCellDataWntCellCycleModel.hpp"

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "CellProliferativeTypesCountWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "VoronoiDataWriter.hpp"

#include "OffLatticeSimulationWithMonoclonalStoppingEvent.hpp"

#include "NagaiHondaForce.hpp"
#include "BoundaryForce.hpp"
#include "RepulsionForce.hpp"
#include "DiffusionCaUpdateRule.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"

#include "SimpleTargetAreaModifier.hpp"
#include "VolumeTrackingModifier.hpp"
#include "WntConcentrationModifier.hpp"
#include "SmoothVertexEdgesModifier.hpp"
#include "JaggedVertexEdgesModifier.hpp"
#include "CurvedVertexEdgesModifier.hpp"
#include "RandomDirectionVertexBasedDivisionRule.hpp"


#include "PlaneBasedCellKiller.hpp"

#include "PlaneBoundaryCondition.hpp"

#include "CryptShovingCaBasedDivisionRule.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Warnings.hpp"
#include "Debug.hpp"

/*
 *  This is where you can set parameters to be used in all the simulations.
 */

static const double M_END_STEADY_STATE = 1; //100
static const double M_END_TIME = 10; //1100
static const double M_DT_TIME = 0.001;//0.001;
static const double M_SAMPLE_TIME = 100;
static const double M_CRYPT_DIAMETER = 6; //16
static const double M_CRYPT_LENGTH = 6;
static const double M_CONTACT_INHIBITION_LEVEL = 0.8;
static const double M_BOUNDARY_FORCE_STRENGTH = 50.0;
static const double M_BOUNDARY_FORCE_CUTOFF = 1.0;


static const std::string M_HEAD_FOLDER = "CylindricalCrypt";

class TestCylindricalCrypt : public AbstractCellBasedWithTimingsTestSuite
{
private:

    /*
     * This is a helper method to generate cells and is used in all simulations.
     */ 

    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells, double equilibriumVolume, double quiescentVolumeFraction)
    {
        double typical_cell_cycle_duration = 12.0;

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        for (unsigned i=0; i<num_cells; i++)
        {
            SimpleCellDataWntCellCycleModel* p_model = new SimpleCellDataWntCellCycleModel();
            p_model->SetDimension(2);
            // p_model->SetEquilibriumVolume(equilibriumVolume);
            // p_model->SetQuiescentVolumeFraction(quiescentVolumeFraction);
            // p_model->SetWntThreshold(0.5);


            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_cell_type);
            double birth_time = - RandomNumberGenerator::Instance()->ranf() * typical_cell_cycle_duration;
            p_cell->SetBirthTime(birth_time);

            // Set Target Area so dont need to use a growth model in vertex simulations
            p_cell->GetCellData()->SetItem("target area", 1.0);

            rCells.push_back(p_cell);
        }
    }

public:

    
    void TestCylindricalCrypts()
    {
        
        //Command line argument stuff: Get the seed parameters  
		// TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-run_index"));
        // unsigned start_index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-run_index");
		
		// TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-num_runs"));
        // unsigned num_runs = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-num_runs");
		
        unsigned start_index = 0;
        unsigned num_runs = 2;

        // Loop over the random seed.
		for(unsigned sim_index=start_index; sim_index < start_index + num_runs; sim_index++)
		{
			std::cout << " Run number " << sim_index << "... \n" << std::flush;

			// Reseed the random number generator
			RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
			
			std::stringstream out;
			out << sim_index;

            double interaction_radi[3] = {1.0,1.5,2.0};
            for(unsigned interation_index=0; interation_index<3; interation_index++)
            {
                std::stringstream out_2;
			    out_2 << interaction_radi[interation_index];

                std::cout << "OS - " << interaction_radi[interation_index] << "\n" << std::flush;

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Node/Radi_"  + out_2.str() + "/Run_" +  out.str();

                // Create a simple mesh
                HoneycombMeshGenerator generator(M_CRYPT_DIAMETER, M_CRYPT_LENGTH, 0);
                TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

                double cut_off_length = interaction_radi[interation_index];

                // Convert this to a Cylindrical2dNodesOnlyMesh
                Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(M_CRYPT_DIAMETER);
                p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh,2.0); // So factor of 16

                // Create cells
                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumNodes(), cells, M_PI*0.25,M_CONTACT_INHIBITION_LEVEL); // mature volume: M_PI*0.25 as r=0.5

                // Create a node-based cell population
                NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
                cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAncestorWriter>();

                for (unsigned index = 0; index < cell_population.rGetMesh().GetNumNodes(); index++)
                {
                    cell_population.rGetMesh().GetNode(index)->SetRadius(0.5);
                }

                // Create simulation from cell population
                OffLatticeSimulationWithMonoclonalStoppingEvent simulator(cell_population);
                simulator.SetStartTime(M_END_STEADY_STATE);
                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_STEADY_STATE);
                simulator.SetOutputDirectory(output_directory);
                simulator.SetOutputDivisionLocations(true);
                simulator.SetOutputCellVelocities(true);

                //Add Wnt concentration modifier
                MAKE_PTR(WntConcentrationModifier<2>, p_wnt_modifier);
                p_wnt_modifier->SetType(LINEAR);
                p_wnt_modifier->SetCryptLength(M_CRYPT_LENGTH);
                simulator.AddSimulationModifier(p_wnt_modifier);

                // Add volume tracking modifier
                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Create a force law and pass it to the simulation
                MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
                p_linear_force->SetMeinekeSpringStiffness(50.0);
                p_linear_force->SetCutOffLength(cut_off_length);
                simulator.AddForce(p_linear_force);

                MAKE_PTR(BoundaryForce<2>, p_boundary_force);
                p_boundary_force->SetForceStrength(M_BOUNDARY_FORCE_STRENGTH);
                p_boundary_force->SetCutOffHeight(M_BOUNDARY_FORCE_CUTOFF);
                simulator.AddForce(p_boundary_force);

                // Sloughing killer
                MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, (M_CRYPT_LENGTH-0.5)*unit_vector<double>(2,1), unit_vector<double>(2,1)));
                simulator.AddCellKiller(p_killer);

                // Solve to Steady State
                try
                {
                    simulator.Solve();
                }
                catch (Exception& e)
                {
                    std::cout << "\n Node Run " << sim_index << " Aborted \n" << std::flush;
                }

                // Mark Ancestors
                simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

                // Reset end time for simulation and run for a further duration
                simulator.SetEndTime(M_END_TIME);

                try
                {
                    simulator.Solve();
                }
                catch (Exception& e)
                {
                    std::cout << "\n Node Run " << sim_index << " Aborted \n" << std::flush;
                }
                
                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

                // Clear memory
                delete p_mesh;
            }

            /*
             * == No ghosts Infinite VT == 
             */
            {
                std::cout << "Infinite VT \n" << std::flush;

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Mesh/NoGhosts/Infinite/Run_" +  out.str();

                // Create mesh (no ghosts)
                CylindricalHoneycombMeshGenerator generator(M_CRYPT_DIAMETER, M_CRYPT_LENGTH, 0); // No Ghosts
                Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

                /////// Change the halo parameters so still get longer edges.  //////
                p_mesh->SetHaloScalingFactor(0.5);
                p_mesh->SetHaloOffset(2.0);

                // Create cells
                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumNodes(),cells,sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL);  //mature_volume = sqrt(3.0)/2.0

                // Create tissue
                MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
                cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAncestorWriter>();

                // Output Voroni for visualisation
                cell_population.AddPopulationWriter<VoronoiDataWriter>();
                cell_population.SetWriteVtkAsPoints(true);

                // Create simulation from cell population
                OffLatticeSimulationWithMonoclonalStoppingEvent simulator(cell_population);
                simulator.SetStartTime(M_END_STEADY_STATE);
                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_STEADY_STATE);
                simulator.SetOutputDirectory(output_directory);
                simulator.SetOutputDivisionLocations(true);
                simulator.SetOutputCellVelocities(true);

                //Add Wnt concentration modifier
                MAKE_PTR(WntConcentrationModifier<2>, p_wnt_modifier);
                p_wnt_modifier->SetType(LINEAR);
                p_wnt_modifier->SetCryptLength(M_CRYPT_LENGTH);
                simulator.AddSimulationModifier(p_wnt_modifier);

                // Add volume tracking Modifier
                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Create a force law and pass it to the simulation
                MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
                p_linear_force->SetMeinekeSpringStiffness(50.0);
                simulator.AddForce(p_linear_force);

                MAKE_PTR(BoundaryForce<2>, p_boundary_force);
                p_boundary_force->SetForceStrength(M_BOUNDARY_FORCE_STRENGTH);
                p_boundary_force->SetCutOffHeight(M_BOUNDARY_FORCE_CUTOFF);
                simulator.AddForce(p_boundary_force);

                // Sloughing killer
                MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, (M_CRYPT_LENGTH-0.5)*unit_vector<double>(2,1), unit_vector<double>(2,1)));
                simulator.AddCellKiller(p_killer);

                // Solve to Steady State
                try
                {
                    simulator.Solve();
                }
                catch (Exception& e)
                {
                    std::cout << "\n Infinite VT Run " << sim_index << " Aborted \n" << std::flush;
                }

                // Mark Ancestors
                simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

                // Reset end time for simulation and run for a further duration
                simulator.SetEndTime(M_END_TIME);

                try
                {
                    simulator.Solve();
                }
                catch (Exception& e)
                {
                    std::cout << "\n Infinite VT Run " << sim_index << " Aborted \n" << std::flush;
                }
                
                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
            }

            /*
             * == No ghosts Finite VT == 
             */
            {
                std::cout << "Finite VT \n" << std::flush;

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Mesh/NoGhosts/Finite/Run_" +  out.str();

                // Create mesh (no ghosts)
                CylindricalHoneycombMeshGenerator generator(M_CRYPT_DIAMETER, M_CRYPT_LENGTH, 0); // No Ghosts

                Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

                // Create cells
                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumNodes(),cells,sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL);  //mature_volume = sqrt(3.0)/2.0

                // Create tissue
                MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
                cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAncestorWriter>();

                // Output Voroni for visualisation
                cell_population.AddPopulationWriter<VoronoiDataWriter>();
                cell_population.SetWriteVtkAsPoints(true);

                ////// Bound VT //////
                cell_population.SetBoundVoronoiTessellation(true);

                // Create simulation from cell population
                OffLatticeSimulationWithMonoclonalStoppingEvent simulator(cell_population);
                simulator.SetStartTime(M_END_STEADY_STATE);
                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_STEADY_STATE);
                simulator.SetOutputDirectory(output_directory);
                simulator.SetOutputDivisionLocations(true);
                simulator.SetOutputCellVelocities(true);

                //Add Wnt concentration modifier
                MAKE_PTR(WntConcentrationModifier<2>, p_wnt_modifier);
                p_wnt_modifier->SetType(LINEAR);
                p_wnt_modifier->SetCryptLength(M_CRYPT_LENGTH);
                simulator.AddSimulationModifier(p_wnt_modifier);

                // Add volume tracking Modifier
                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Create a force law and pass it to the simulation
                MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
                p_linear_force->SetMeinekeSpringStiffness(50.0);

                ///////  Add Cut off length to force //////
                p_linear_force->SetCutOffLength(1.5);

                simulator.AddForce(p_linear_force);
            
                MAKE_PTR(BoundaryForce<2>, p_boundary_force);
                p_boundary_force->SetForceStrength(M_BOUNDARY_FORCE_STRENGTH);
                p_boundary_force->SetCutOffHeight(M_BOUNDARY_FORCE_CUTOFF);
                simulator.AddForce(p_boundary_force);

                // Sloughing killer
                MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, (M_CRYPT_LENGTH-0.5)*unit_vector<double>(2,1), unit_vector<double>(2,1)));
                simulator.AddCellKiller(p_killer);

                // Solve to Steady State
                try
                {
                    simulator.Solve();
                }
                catch (Exception& e)
                {
                    std::cout << "\n Finite VT Run " << sim_index << " Aborted \n" << std::flush;
                }

                // Mark Ancestors
                simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

                // Reset end time for simulation and run for a further duration
                simulator.SetEndTime(M_END_TIME);

                try
                {
                    simulator.Solve();
                }
                catch (Exception& e)
                {
                    std::cout << "\n Finite VT Run " << sim_index << " Aborted \n" << std::flush;
                }
                
                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
            }

            /*
             * == Ghosts VT == 
             */
            {
                std::cout << "Ghost VT \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Mesh/Ghosts/Run_" +  out.str();

                // Create mesh
                unsigned thickness_of_ghost_layer = 2;

                CylindricalHoneycombMeshGenerator generator(M_CRYPT_DIAMETER, M_CRYPT_LENGTH, thickness_of_ghost_layer);
                Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

                // Get location indices corresponding to real cells
                std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

                // Create cells
                std::vector<CellPtr> cells;
                GenerateCells(location_indices.size(),cells,sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL);  //mature_volume = sqrt(3.0)/2.0

                // Create tissue
                MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
                cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAncestorWriter>();

                // Output Voroni for visualisation
                cell_population.AddPopulationWriter<VoronoiDataWriter>();
                cell_population.SetWriteVtkAsPoints(true);

                // Create simulation from cell population
                OffLatticeSimulationWithMonoclonalStoppingEvent simulator(cell_population);
                simulator.SetStartTime(M_END_STEADY_STATE);
                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_STEADY_STATE);
                simulator.SetOutputDirectory(output_directory);
                simulator.SetOutputDivisionLocations(true);
                simulator.SetOutputCellVelocities(true);

                //Add Wnt concentration modifier
                MAKE_PTR(WntConcentrationModifier<2>, p_wnt_modifier);
                p_wnt_modifier->SetType(LINEAR);
                p_wnt_modifier->SetCryptLength(M_CRYPT_LENGTH);
                simulator.AddSimulationModifier(p_wnt_modifier);

                // Add volume tracking Modifier
                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Create a force law and pass it to the simulation
                MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
                p_linear_force->SetMeinekeSpringStiffness(50.0);
                simulator.AddForce(p_linear_force);

                MAKE_PTR(BoundaryForce<2>, p_boundary_force);
                p_boundary_force->SetForceStrength(M_BOUNDARY_FORCE_STRENGTH);
                p_boundary_force->SetCutOffHeight(M_BOUNDARY_FORCE_CUTOFF);
                simulator.AddForce(p_boundary_force);

                // Sloughing killer
                MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, (M_CRYPT_LENGTH-0.5)*unit_vector<double>(2,1), unit_vector<double>(2,1)));
                simulator.AddCellKiller(p_killer);

                // Run simulation
                try
                {
                    simulator.Solve();
                }
                catch (Exception& e)
                {
                    std::cout << "\n Ghost VT Run " << sim_index << " Aborted \n" << std::flush;
                }

                // Mark Ancestors
                simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

                // Reset end time for simulation and run for a further duration
                simulator.SetEndTime(M_END_TIME);

                try
                {
                    simulator.Solve();
                }
                catch (Exception& e)
                {
                    std::cout << "\n Ghost VT Run " << sim_index << " Aborted \n" << std::flush;
                }
                
                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
            }

            /*
            * == VM Smooth ==
            */
            {
                std::cout << "Vertex Smooth \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Vertex/Smooth/Run_" +  out.str();

                // Create mesh
                bool is_flat_bottom = true; // only different here with jagged is this.
                CylindricalHoneycombVertexMeshGenerator generator(M_CRYPT_DIAMETER, M_CRYPT_LENGTH, is_flat_bottom);
                Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();
                p_mesh->SetCellRearrangementThreshold(0.05);
                p_mesh->SetCellRearrangementRatio(1.5);

                // Create cells
                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumElements(),cells,1.0,M_CONTACT_INHIBITION_LEVEL); //mature_volume = 1.0

                // Create tissue
                VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
                cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellIdWriter>();

                // Create crypt simulation from cell population
                OffLatticeSimulationWithMonoclonalStoppingEvent simulator(cell_population);
                simulator.SetStartTime(M_END_STEADY_STATE);
                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_STEADY_STATE);
                simulator.SetOutputDirectory(output_directory);
                simulator.SetOutputDivisionLocations(true);
                simulator.SetOutputCellVelocities(true);
                cell_population.AddCellWriter<CellAncestorWriter>();

                //Add Wnt concentration modifier
                MAKE_PTR(WntConcentrationModifier<2>, p_wnt_modifier);
                p_wnt_modifier->SetType(LINEAR);
                p_wnt_modifier->SetCryptLength(M_CRYPT_LENGTH);
                simulator.AddSimulationModifier(p_wnt_modifier);

                // Add volume tracking modifier
                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Refine the edges on boundary to get smooth edges
                MAKE_PTR(SmoothVertexEdgesModifier<2>, smooth_edge_modifier);
                simulator.AddSimulationModifier(smooth_edge_modifier);

                // Create Forces and pass to simulation NOTE : these are not the default ones and chosen to give a stable growing monolayer
                MAKE_PTR(NagaiHondaForce<2>, p_force);
                p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
                p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
                p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
                p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0);
                simulator.AddForce(p_force);
                
                MAKE_PTR(BoundaryForce<2>, p_boundary_force);
                p_boundary_force->SetForceStrength(M_BOUNDARY_FORCE_STRENGTH);
                p_boundary_force->SetCutOffHeight(M_BOUNDARY_FORCE_CUTOFF);
                simulator.AddForce(p_boundary_force);

                // Sloughing killer
                MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, M_CRYPT_LENGTH*unit_vector<double>(2,1), unit_vector<double>(2,1)));
                simulator.AddCellKiller(p_killer);

                // Run simulation
                try
                {
                    simulator.Solve();
                }
                catch (Exception& e)
                {
                    std::cout << "\n Vertex Smooth Run " << sim_index << " Aborted \n" << std::flush;
                }

                // Mark Ancestors
                simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

                // Reset end time for simulation and run for a further duration
                simulator.SetEndTime(M_END_TIME);

                try
                {
                    simulator.Solve();
                }
                catch (Exception& e)
                {
                    std::cout << "\n Vertex Smooth Run " << sim_index << " Aborted \n" << std::flush;
                }
                
                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
                
                // Clear singletons
                Warnings::Instance()->QuietDestroy();
            }

            /*
            * == VM Jagged ==
            */
            {
                std::cout << "Vertex Jagged \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Vertex/Jagged/Run_" +  out.str();
      
                // Create mesh
                bool is_flat_bottom = false; // only different here with smoothed is this.
                CylindricalHoneycombVertexMeshGenerator generator(M_CRYPT_DIAMETER, M_CRYPT_LENGTH, is_flat_bottom);
                Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();
                p_mesh->SetCellRearrangementThreshold(0.05);
                p_mesh->SetCellRearrangementRatio(1.5);

                // Create cells
                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumElements(),cells,1.0,M_CONTACT_INHIBITION_LEVEL); //mature_volume = 1.0

                // Create tissue
                VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
                cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellIdWriter>();

                // Create crypt simulation from cell population
                OffLatticeSimulationWithMonoclonalStoppingEvent simulator(cell_population);
                simulator.SetStartTime(M_END_STEADY_STATE);
                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_STEADY_STATE);
                simulator.SetOutputDirectory(output_directory);
                simulator.SetOutputDivisionLocations(true);
                simulator.SetOutputCellVelocities(true);
                cell_population.AddCellWriter<CellAncestorWriter>();

                //Add Wnt concentration modifier
                MAKE_PTR(WntConcentrationModifier<2>, p_wnt_modifier);
                p_wnt_modifier->SetType(LINEAR);
                p_wnt_modifier->SetCryptLength(M_CRYPT_LENGTH);
                simulator.AddSimulationModifier(p_wnt_modifier);

                // Add volume tracking modifier
                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Refine the edges on boundary to get accurate edges
                MAKE_PTR(JaggedVertexEdgesModifier<2>, refinement_modifier);
                simulator.AddSimulationModifier(refinement_modifier);

                // Create Forces and pass to simulation NOTE : these are not the default ones and chosen to give a stable growing monolayer
                MAKE_PTR(NagaiHondaForce<2>, p_force);
                p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
                p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
                p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
                p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0);
                simulator.AddForce(p_force);
                
                MAKE_PTR(BoundaryForce<2>, p_boundary_force);
                p_boundary_force->SetForceStrength(M_BOUNDARY_FORCE_STRENGTH);
                p_boundary_force->SetCutOffHeight(M_BOUNDARY_FORCE_CUTOFF);
                simulator.AddForce(p_boundary_force);

                // Sloughing killer
                MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, M_CRYPT_LENGTH*unit_vector<double>(2,1), unit_vector<double>(2,1)));
                simulator.AddCellKiller(p_killer);

                // Run simulation
                try
                {
                    simulator.Solve();
                }
                catch (Exception& e)
                {
                    std::cout << "\n Vertex Jagged Run " << sim_index << " Aborted \n" << std::flush;
                }

                // Mark Ancestors
                simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

                // Reset end time for simulation and run for a further duration
                simulator.SetEndTime(M_END_TIME);

                try
                {
                    simulator.Solve();
                }
                catch (Exception& e)
                {
                    std::cout << "\n Vertex Jagged Run " << sim_index << " Aborted \n" << std::flush;
                }
                
                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
                
                // Clear singletons
                Warnings::Instance()->QuietDestroy();
            }

            /*
            * == VM Curved ==
            */
            {
                std::cout << "Vertex Curved \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Vertex/Curved/Run_" +  out.str();

                // Create mesh
                CylindricalHoneycombVertexMeshGenerator generator(M_CRYPT_DIAMETER, M_CRYPT_LENGTH, false);
                Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();
                p_mesh->SetCellRearrangementThreshold(0.05);
                p_mesh->SetCellRearrangementRatio(1.5);
                
                // Create cells
                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumElements(),cells,1.0,M_CONTACT_INHIBITION_LEVEL); //mature_volume = 1.0

                // Create tissue
                VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
                cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellIdWriter>();

                // Create crypt simulation from cell population
                OffLatticeSimulationWithMonoclonalStoppingEvent simulator(cell_population);
                simulator.SetStartTime(M_END_STEADY_STATE);
                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_STEADY_STATE);
                simulator.SetOutputDirectory(output_directory);
                simulator.SetOutputDivisionLocations(true);
                simulator.SetOutputCellVelocities(true);
                cell_population.AddCellWriter<CellAncestorWriter>();

                //Add Wnt concentration modifier
                MAKE_PTR(WntConcentrationModifier<2>, p_wnt_modifier);
                p_wnt_modifier->SetType(LINEAR);
                p_wnt_modifier->SetCryptLength(M_CRYPT_LENGTH);
                simulator.AddSimulationModifier(p_wnt_modifier);

                // Add volume tracking modifier
                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Refine the edges on boundary to get smooth edges, and fix T1, T2 and T3 swaps:
                MAKE_PTR(CurvedVertexEdgesModifier<2>, refinement_modifier);
                simulator.AddSimulationModifier(refinement_modifier);

                // Create Forces and pass to simulation NOTE : these are not the default ones and chosen to give a stable growing monolayer
                MAKE_PTR(NagaiHondaForce<2>, p_force);
                p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
                p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
                p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
                p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0);
                simulator.AddForce(p_force);
                
                MAKE_PTR(BoundaryForce<2>, p_boundary_force);
                p_boundary_force->SetForceStrength(M_BOUNDARY_FORCE_STRENGTH);
                p_boundary_force->SetCutOffHeight(M_BOUNDARY_FORCE_CUTOFF);
                simulator.AddForce(p_boundary_force);

                // Sloughing killer
                MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, M_CRYPT_LENGTH*unit_vector<double>(2,1), unit_vector<double>(2,1)));
                simulator.AddCellKiller(p_killer);

                boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_division_rule_to_set(new RandomDirectionVertexBasedDivisionRule<2>());
                cell_population.SetVertexBasedDivisionRule(p_division_rule_to_set);

                // Run simulation
                try
                {
                    simulator.Solve();
                }
                catch (Exception& e)
                {
                    std::cout << "\n Vertex Curved Run " << sim_index << " Aborted \n" << std::flush;
                }

                // Mark Ancestors
                simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

                // Reset end time for simulation and run for a further duration
                simulator.SetEndTime(M_END_TIME);

                try
                {
                    simulator.Solve();
                }
                catch (Exception& e)
                {
                    std::cout << "\n Vertex Curved Run " << sim_index << " Aborted \n" << std::flush;
                }
                
                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
                
                // Clear singletons
                Warnings::Instance()->QuietDestroy();
            }
        }
    }
};

#endif /* TESTCYLINDRICALCRYPT_HPP_ */

