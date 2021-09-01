#ifndef TESTGROWINGMONOLAYER_HPP_
#define TESTGROWINGMONOLAYER_HPP_


#include <cxxtest/TestSuite.h> 

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp" 

#include "SmartPointers.hpp"

#include "HoneycombMeshGenerator.hpp"
// #include "MutableMesh.hpp"
// #include "NodesOnlyMesh.hpp"
#include "HoneycombVertexMeshGenerator.hpp"

#include "CellsGenerator.hpp"

#include "ContactInhibitionCellCycleModel.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "CellProliferativeTypesCountWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellVolumesWriter.hpp"

#include "OffLatticeSimulation.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"

#include "VolumeTrackingModifier.hpp"
#include "SimpleTargetAreaModifier.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Warnings.hpp"

/*
 *  This is where you can set parameters to be used in all the simulations.
 */

static const double M_END_TIME = 100; //100
static const double M_CONTACT_INHIBITION_LEVEL = 0.8;

class TestGrowingMonolayer : public AbstractCellBasedWithTimingsTestSuite
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
            ContactInhibitionCellCycleModel* p_model = new ContactInhibitionCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetEquilibriumVolume(equilibriumVolume);
            p_model->SetQuiescentVolumeFraction(quiescentVolumeFraction);


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
     * Simulate a growing monolayer using the
     * Overlapping Spheres model.
     */
    void TestNodeBasedGrowingMonolayer()
    {
        /* 
         * == Default cut-off ==
         */
         // Create a simple mesh
        HoneycombMeshGenerator generator(2, 2);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        double cut_off_length = 1.5; //this is the default

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0);

        // Create cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(), cells, M_PI*0.25,M_CONTACT_INHIBITION_LEVEL); // mature volume: M_PI*0.25 as r=0.5

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(0.005);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_END_TIME); //50
        simulator.SetOutputDirectory("GrowingSpheroid/Node/Default");
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        p_linear_force->SetCutOffLength(cut_off_length);
        simulator.AddForce(p_linear_force);

        // Run simulation
        simulator.Solve();

        // Clear memory
        delete p_mesh;

        /*
         * == Larger cut-off ==
         */

    }

    /* 
     * == VT ==
     *
     * Simulate a growing monolayer using the
     * Voronoi Tessellation model.
     */
    void TestMeshBasedGrowingMonolayer()
    {
        /* 
         * == No ghosts == 
         */
         // Create mesh
        HoneycombMeshGenerator generator(2, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        GenerateCells(location_indices.size(),cells,sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL);  //mature_volume = sqrt(3.0)/2.0

        // Create tissue
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(0.005);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_END_TIME); //50
        simulator.SetOutputDirectory("GrowingSpheroid/Mesh/NoGhosts");
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        simulator.AddForce(p_linear_force);

        // Run simulation
        simulator.Solve();

        /*
         * == Ghosts == 
         */

    }

    /* 
     * == VM ==
     *
     * Simulate a growing monolayer using the
     * Cell Vertex model.
     */
    void TestVertexBasedGrowingMonolayer() 
    {
        /* 
         * == Jagged edge == 
         */
        // Create cells
        HoneycombVertexMeshGenerator mesh_generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = mesh_generator.GetMesh();
        // p_mesh->SetCellRearrangementThreshold(0.1);

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
        simulator.SetDt(0.005);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_END_TIME); //50
        simulator.SetOutputDirectory("GrowingSpheroid/Vertex/Jagged");
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create Forces and pass to simulation NOTE : these are not the default ones and chosen to give a stable growing monolayer
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
        simulator.AddForce(p_force);

        // Add target area modifier
        // MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        // p_growth_modifier->SetGrowthDuration(0);
        // simulator.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulator.Solve();

        /*
         * == Smooth edge == 
         */

        /*
         * == Curved edge == 
         */

    }
};

#endif /* TESTGROWINGMONOLAYER_HPP_ */
