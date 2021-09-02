#ifndef TESTINTERNALVOID_HPP_
#define TESTINTERNALVOID_HPP_


#include <cxxtest/TestSuite.h> 

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp" 
#include "CheckpointArchiveTypes.hpp"

#include "SmartPointers.hpp"

#include "HoneycombMeshGenerator.hpp"
#include "PeriodicNodesOnlyMesh.hpp"
#include "MutableMesh.hpp"
#include "ToroidalHoneycombVertexMeshGenerator.hpp"

#include "CellsGenerator.hpp"

#include "NoCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "CellVolumesWriter.hpp"

#include "OffLatticeSimulation.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"

#include "SimpleTargetAreaModifier.hpp"
#include "VolumeTrackingModifier.hpp"

#include "PlaneBoundaryCondition.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp" 
#include "PetscSetupAndFinalize.hpp"
#include "Warnings.hpp"

/*
 *  This is where you can set parameters to be used in all the simulations.
 */

static const double M_END_STEADY_STATE = 10; //25
static const double M_END_TIME = 60; //50
static const double M_DOMAIN_WIDTH = 12;
static const double M_DOMAIN_LENGTH = 12;

class TestInternalVoid : public AbstractCellBasedWithTimingsTestSuite
{
private:
    /**
    * Helper method. Smooth out edges of a vertex mesh.
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
    * Helper method. Iterate over all cells and define the 'hole' by
    * killing those cells whose centres are located in a given region.
    * 
    * @param rCellPopulation a cell population
    * @param holeWidth the width of the hole
    * @param xMin the left boundary of the hole
    * @param xMax the right boundary of the hole
    * @param yMin the bottom boundary of the hole
    * @param yMax the top boundary of the hole
    */
    void CreateHoleInCellPopulation(AbstractCellPopulation<2>& rCellPopulation,
                                    double holeWidth,
                                    double xMin,
                                    double xMax,
                                    double yMin,
                                    double yMax)
    {
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
        {
            // Get the coordinates of this cell centre
            c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            double x = centre_of_cell[0];
            double y = centre_of_cell[1];

            if ((fabs(y-x)<holeWidth) && (x>xMin) && (x<xMax) && (y>yMin) && (y<yMax))
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

public:

    /* 
     * == OS ==
     *
     * Simulate an internal void using the
     * Overlapping Spheres model.
     */
    void TestNodeBasedInternalVoid()
    {
        /* 
         * == Pre-void == 
         */
         // Create simple mesh
        HoneycombMeshGenerator generator(M_DOMAIN_WIDTH, M_DOMAIN_LENGTH, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();
        p_generating_mesh->Scale(0.8, 0.8);

        double cut_off_length = 1.5; //this is the default

        // Convert this to a PeriodicNodesOnlyMesh
        c_vector<double,2> periodic_width = zero_vector<double>(2);
        periodic_width[0] = 10.0;
        periodic_width[1] = 10.0;
        PeriodicNodesOnlyMesh<2>* p_mesh = new PeriodicNodesOnlyMesh<2>(periodic_width);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(0.005);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory("InternalVoid/Node/Pre-void");
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

        // Save simulation in steady state
		CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        // Clear memory
        delete p_mesh;

        /*
         * == Void default cut-off == 
         */
        // Load steady state
        OffLatticeSimulation<2>* p_simulator1 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("InternalVoid/Node/Pre-void",M_END_STEADY_STATE);
        NodeBasedCellPopulation<2>* cell_population1 = static_cast<NodeBasedCellPopulation<2>*>(&(p_simulator1->rGetCellPopulation()));

        // Now remove cells in a given region using a helper method
        CreateHoleInCellPopulation(*cell_population1, 1.75, 2.0, 8.0, 2.0, 8.0);

        // Reset timestep, sampling timestep and end time for simulation and run for a further duration
        p_simulator1->SetDt(0.005);
        p_simulator1->SetSamplingTimestepMultiple(200);
        p_simulator1->SetEndTime(M_END_TIME);
        p_simulator1->SetOutputDirectory("InternalVoid/Node/DefaultCutOff");
        p_simulator1->Solve();

        // Tidy up
        delete p_simulator1;

        /*
         * == Larger cut-off ==
         */

    }

    /* 
     * == VT ==
     * 
     * Simulate internal voide using the
     * Voronoi Tesselation model.
     */
    void TestMeshBasedInternalVoid()
    {
        // Create mesh
        unsigned thickness_of_ghost_layer = 2;

        HoneycombMeshGenerator generator(M_DOMAIN_WIDTH, M_DOMAIN_LENGTH, thickness_of_ghost_layer);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->Scale(0.8, 0.8);

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_differentiated_type);

        // Create tissue
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(0.005);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory("InternalVoid/Mesh/Pre-void");
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking Modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        simulator.AddForce(p_linear_force);

        // Create boundary condition y > 0
        c_vector<double, 2> point1 = zero_vector<double>(2);
        c_vector<double, 2> normal1 = zero_vector<double>(2);
        normal1(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_condition1, (&cell_population, point1, normal1));

        // Create boundary condition x > 0
        c_vector<double, 2> point2 = zero_vector<double>(2);
        c_vector<double, 2> normal2 = zero_vector<double>(2);
        normal2(0) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_condition2, (&cell_population, point2, normal2));

        // Create boundary condition y < 10
        c_vector<double, 2> point3 = zero_vector<double>(2);
        point3(1) = 10.0;
        normal1(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_condition3, (&cell_population, point3, normal1));

        // Create boundary condition x < 10
        c_vector<double, 2> point4 = zero_vector<double>(2);
        point4(0) = 10.0;
        normal2(0) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_condition4, (&cell_population, point4, normal2));

        // Pass the boundary conditions to the simulation
        simulator.AddCellPopulationBoundaryCondition(p_condition1);
        simulator.AddCellPopulationBoundaryCondition(p_condition2);
        simulator.AddCellPopulationBoundaryCondition(p_condition3);
        simulator.AddCellPopulationBoundaryCondition(p_condition4);

        // Run simulation
        simulator.Solve();

        // Save simulation in steady state
		CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        /*
         * == No ghosts == 
         */
         // Load steady state
        OffLatticeSimulation<2>* p_simulator1 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("InternalVoid/Mesh/Pre-void",M_END_STEADY_STATE);
        MeshBasedCellPopulation<2>* cell_population1 = static_cast<MeshBasedCellPopulation<2>*>(&(p_simulator1->rGetCellPopulation()));

        // Now remove cells in a given region using a helper method
        CreateHoleInCellPopulation(*cell_population1, 1.75, 2.0, 8.0, 2.0, 8.0);

        // Reset timestep, sampling timestep and end time for simulation and run for a further duration
        p_simulator1->SetDt(0.005);
        p_simulator1->SetOutputDirectory("InternalVoid/Mesh/NoGhosts");
        p_simulator1->SetSamplingTimestepMultiple(200);
        p_simulator1->SetEndTime(M_END_TIME);
        p_simulator1->Solve();

        // Tidy up
        delete p_simulator1;

        /*
         * == Ghosts ==
         */

    }

    /* 
     * == VM ==
     * 
     * Simulation internal void using the
     * Cell Vertex model.
     */
    void TestVertexBasedInternalVoid()
    {
        /* 
         * == Pre-void == 
         */
         // Create mesh
        ToroidalHoneycombVertexMeshGenerator generator(M_DOMAIN_WIDTH, M_DOMAIN_LENGTH);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();
        p_mesh->Scale(0.8, 0.8);
        p_mesh->SetHeight(10);
        p_mesh->SetWidth(10);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        // Create tissue
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(0.005);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory("InternalVoid/Vertex/Pre-void");
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
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->SetGrowthDuration(0);
        simulator.AddSimulationModifier(p_growth_modifier);
        
        // Smooth out edges to get nice box domain
        SmoothVertexMeshEdges(cell_population);

        // Run simulation
        simulator.Solve();

        // Save simulation in steady state
		CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        /*
         * == Smooth void == 
         */
         // Load steady state
        OffLatticeSimulation<2>* p_simulator1 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("InternalVoid/Vertex/Pre-void",M_END_STEADY_STATE);
        VertexBasedCellPopulation<2>* cell_population1 = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator1->rGetCellPopulation()));

        // Now remove cells in a given region using a helper method
        CreateHoleInCellPopulation(*cell_population1, 1.75, 2.0, 8.0, 2.0, 8.0);
        SmoothVertexMeshEdges(*cell_population1);

        // Reset timestep, sampling timestep and end time for simulation and run for a further duration
        p_simulator1->SetDt(0.005);
        p_simulator1->SetSamplingTimestepMultiple(200);
        p_simulator1->SetEndTime(M_END_TIME);
        p_simulator1->SetOutputDirectory("InternalVoid/Vertex/Smooth");
        p_simulator1->Solve();

        // Tidy up
        delete p_simulator1;

        /*
         * == Jagged void ==
         */
         // Load steady state
        OffLatticeSimulation<2>* p_simulator2 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("InternalVoid/Vertex/Pre-void",M_END_STEADY_STATE);
        VertexBasedCellPopulation<2>* cell_population2 = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator2->rGetCellPopulation()));
        // Now remove cells in a given region using a helper method
        CreateHoleInCellPopulation(*cell_population2, 1.75, 2.0, 8.0, 2.0, 8.0);

        // Reset timestep, sampling timestep and end time for simulation and run for a further duration
        p_simulator2->SetDt(0.005);
        p_simulator2->SetSamplingTimestepMultiple(200);
        p_simulator2->SetEndTime(M_END_TIME);
        p_simulator2->SetOutputDirectory("InternalVoid/Vertex/Jagged");
        p_simulator2->Solve();

        // Tidy up
        delete p_simulator2;

        /*
         * == Curved edges ==
         */

    }
};

#endif /* TESTINTERNALVOID_HPP_ */