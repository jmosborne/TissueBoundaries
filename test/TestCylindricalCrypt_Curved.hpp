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
#include "BoundaryNodeWriter.hpp"
#include "BoundaryCellWriter.hpp"

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
     * == VM ==
     *
     * Simulate cell proliferation in the colorectal crypt using the
     * Cell Vertex model with curved edges.
     */
    void TestVertexBasedCurvedCrypt()
    {

        RandomNumberGenerator::Instance()->Reseed(1);

        std::string output_directory = M_HEAD_FOLDER + "/Vertex/Curved_test_1";

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
        cell_population.AddPopulationWriter<BoundaryNodeWriter>();

        // cell_population.AddCellWriter<BoundaryCellWriter>();


        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();

        // Create crypt simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(M_DT_TIME);
        // simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
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

        MAKE_PTR(VertexEdgesModifier<2>, edge_modifier);
        simulator.AddSimulationModifier(edge_modifier);


        simulator.SetSamplingTimestepMultiple(200);

        simulator.SetEndTime(480);
        
        // Mark Ancestors
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Run simulation to new end time
        simulator.Solve();

        // Clear singletons
        Warnings::Instance()->QuietDestroy();
    }
};

#endif /* TESTCYLINDRICALCRYPT_HPP_ */

