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

#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"

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
#include "VertexBoundaryRefinementModifier.hpp"
#include "SmoothVertexEdgesModifier.hpp"
#include "VertexEdgesModifier.hpp"
#include "RandomDirectionVertexBasedDivisionRule.hpp"


#include "PlaneBasedCellKiller.hpp"

#include "PlaneBoundaryCondition.hpp"

#include "CryptShovingCaBasedDivisionRule.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Warnings.hpp"

/*
 *  This is where you can set parameters to be used in all the simulations.
 */

static const double M_END_STEADY_STATE = 100; //100
static const double M_END_TIME = 100; //1100
static const double M_DT_TIME = 0.001;
static const double M_SAMPLE_TIME = 1;
static const double M_CRYPT_DIAMETER = 6; //16
static const double M_CRYPT_LENGTH = 10;
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

    /*
     * == OS ==
     *
     * Simulate cell proliferation in the colorectal crypt using the
     * Overlapping Spheres model.
     */
    void xTestNodeBasedCrypt()
    {
        std::string output_directory = M_HEAD_FOLDER + "/Node/Default";

        // Create a simple mesh
        HoneycombMeshGenerator generator(M_CRYPT_DIAMETER, M_CRYPT_LENGTH, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        double cut_off_length = 1.5; //this is the default

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
        OffLatticeSimulation<2> simulator(cell_population);
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

        // Run simulation
        simulator.Solve();

        // Mark Ancestors
        simulator.SetEndTime(M_END_TIME);
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Run simulation to new end time
        simulator.Solve();

        // Clear memory
        delete p_mesh;
    }

    void noTestNodeBasedLargeCutoffCrypt()
    {
        std::string output_directory = M_HEAD_FOLDER + "/Node/LargeCutoff";

        // Create a simple mesh
        HoneycombMeshGenerator generator(M_CRYPT_DIAMETER, M_CRYPT_LENGTH, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        double cut_off_length = 2.0; //this is larger than the default

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
        OffLatticeSimulation<2> simulator(cell_population);
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

        // Run simulation
        simulator.Solve();

        // Mark Ancestors
        simulator.SetEndTime(M_END_TIME);
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Run simulation to new end time
        simulator.Solve();

        // Clear memory
        delete p_mesh;
    }

    void noTestNodeBasedSmallCutoffCrypt()
    {
        std::string output_directory = M_HEAD_FOLDER + "/Node/SmallCutoff";

        // Create a simple mesh
        HoneycombMeshGenerator generator(M_CRYPT_DIAMETER, M_CRYPT_LENGTH, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        double cut_off_length = 1.0; //this is smaller than the default

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
        OffLatticeSimulation<2> simulator(cell_population);
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

        // Run simulation
        simulator.Solve();

        // Mark Ancestors
        simulator.SetEndTime(M_END_TIME);
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Run simulation to new end time
        simulator.Solve();

        // Clear memory
        delete p_mesh;
    }

    /*
     * == VT ==
     *
     * Simulate cell proliferation in the colorectal crypt using the
     * Voronoi Tesselation model.
     */


    /*
     * == No ghosts Ininite VT == 
     */
    void noTestMeshBasedNoGhostsInfiniteVtCrypt()
    {
        std::string output_directory = M_HEAD_FOLDER + "/Mesh/NoGhosts/Infinite";

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

        // static_cast<MeshBasedCellPopulation<2>*>(&(p_simulator_1->rGetCellPopulation()))->SetBoundVoronoiTessellation(true);


        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
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
        simulator.Solve();

        // Mark Ancestors
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Reset end time for simulation and run for a further duration
        simulator.SetEndTime(M_END_TIME);

        simulator.Solve();
    }

    /*
     * == No ghosts Finite VT == 
     */
    void noTestMeshBasedNoGhostsFiniteVtCrypt()
    {
        std::string output_directory = M_HEAD_FOLDER + "/Mesh/NoGhosts/Finite";

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

        // Bound VT
        cell_population.SetBoundVoronoiTessellation(true);

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
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
        simulator.Solve();

        // Mark Ancestors
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Reset end time for simulation and run for a further duration
        simulator.SetEndTime(M_END_TIME);

        simulator.Solve();
    }

    void xTestMeshBasedGhostsCrypt()
    {
        std::string output_directory = M_HEAD_FOLDER + "/Mesh/Ghosts";

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
        OffLatticeSimulation<2> simulator(cell_population);
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
        simulator.Solve();

        // Mark Ancestors
        simulator.SetEndTime(M_END_TIME);
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Run simulation to new end time
        simulator.Solve();
    }

    /*
     * == VM ==
     *
     * Simulate cell proliferation in the colorectal crypt using the
     * Cell Vertex model. With smoothed edges.
     */
    void xTestVertexBasedSmoothCrypt()
    {
        std::string output_directory = M_HEAD_FOLDER + "/Vertex/Smooth";

        // Create mesh
        bool is_flat_bottom = true; // only different here with jagged is this.
        CylindricalHoneycombVertexMeshGenerator generator(M_CRYPT_DIAMETER, M_CRYPT_LENGTH, is_flat_bottom);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Create cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells,1.0,M_CONTACT_INHIBITION_LEVEL); //mature_volume = 1.0

        // Create tissue
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();

        // Create crypt simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
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
        simulator.Solve();

        // Mark Ancestors
        simulator.SetEndTime(M_END_TIME);
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Run simulation to new end time
        simulator.Solve();

        // Clear singletons
        Warnings::Instance()->QuietDestroy();
    }

    void xTestVertexBasedJaggedCrypt()
    {
        std::string output_directory = M_HEAD_FOLDER + "/Vertex/Jagged";

        // Create mesh
        bool is_flat_bottom = false; // only different here with smoothed is this.
        CylindricalHoneycombVertexMeshGenerator generator(M_CRYPT_DIAMETER, M_CRYPT_LENGTH, is_flat_bottom);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Create cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells,1.0,M_CONTACT_INHIBITION_LEVEL); //mature_volume = 1.0

        // Create tissue
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();

        // Create crypt simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
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
        simulator.Solve();

        // Mark Ancestors
        simulator.SetEndTime(M_END_TIME);
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Run simulation to new end time
        simulator.Solve();

        // Clear singletons
        Warnings::Instance()->QuietDestroy();
    }

    /*
     * == VM ==
     *
     * Simulate cell proliferation in the colorectal crypt using the
     * Cell Vertex model with curved edges.
     */
    void TestVertexBasedCurvedCrypt()
    {

        RandomNumberGenerator::Instance()->Reseed(3);

        std::string output_directory = M_HEAD_FOLDER + "/Vertex/Curved_test_2";

        // Create mesh
        CylindricalHoneycombVertexMeshGenerator generator(M_CRYPT_DIAMETER, M_CRYPT_LENGTH, false);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Create cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells,1.0,M_CONTACT_INHIBITION_LEVEL); //mature_volume = 1.0

        // Create tissue
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();

        // Create crypt simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(M_DT_TIME);
        simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
        // simulator.SetEndTime(M_END_STEADY_STATE);
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
        MAKE_PTR(VertexBoundaryRefinementModifier<2>, refinement_modifier);
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
        // simulator.Solve();

        MAKE_PTR(VertexEdgesModifier<2>, edge_modifier);
        simulator.AddSimulationModifier(edge_modifier);

        // Mark Ancestors
        simulator.SetEndTime(M_END_TIME);
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Run simulation to new end time
        simulator.Solve();

        // Clear singletons
        Warnings::Instance()->QuietDestroy();
    }
};

#endif /* TESTCYLINDRICALCRYPT_HPP_ */

