#ifndef TESTGROWINGMONOLAYER_HPP_
#define TESTGROWINGMONOLAYER_HPP_


#include <cxxtest/TestSuite.h> 

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp" 

#include "SmartPointers.hpp"

#include "HoneycombMeshGenerator.hpp"
#include "MutableMesh.hpp"
#include "NodesOnlyMesh.hpp"
#include "HoneycombVertexMeshGenerator.hpp"

#include "CellsGenerator.hpp"

#include "BernoulliTrialWithContactInhibitionCellCycleModel.hpp"
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
#include "VertexBoundaryRefinementModifier.hpp"
#include "SmoothVertexEdgesModifier.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Warnings.hpp"
#include "Debug.hpp"

#include "CicularityCalcModifier.hpp"
#include "BoundaryCellWriter.hpp"

#include "GhostNodeRemovalModifier.hpp"

/*
 *  This is where you can set parameters to be used in all the simulations.
 */

static const double M_END_TIME = 40; //100
static const double M_DT_TIME = 0.005;
static const double M_SAMPLE_TIME = 50;

static const double M_CONTACT_INHIBITION_LEVEL = 0.8;
static const double M_STEM_CELL_DIVISION_PROBABILITY = 0.1;
static const double M_STEM_CELL_MINIMUM_DIVISION_AGE = 1.0;

static const double M_INITIAL_WIDTH = 10;
static const double M_INITIAL_LENGTH = 10;

static const double M_DOMAIN_WIDTH = 5.0;
static const double M_DOMAIN_X_MIN = -5.0;
static const double M_DOMAIN_X_MAX = 5.0;
static const double M_DOMAIN_Y_MIN = -5.0;
static const double M_DOMAIN_Y_MAX = 5.0;

static const std::string M_HEAD_FOLDER = "GrowingMonolayerSweep_40hrs";

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

    /**
    * Helper method. Iterate over all cells and 'round out' the domain by
    * killing those cells whose centres are located outside a given region.
    * 
    * @param rCellPopulation a cell population
    */
    void RoundOutCellPopulation(AbstractCellPopulation<2>& rCellPopulation)
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

                unsigned location_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
 
                // if ((fabs(y-x)>M_DOMAINWIDTH) && (x<M_DOMAIN_X_MIN) && (x>M_DOMAIN_X_MAX) && (y<M_DOMAIN_Y_MIN) && (y>M_DOMAIN_Y_MAX))
                if ( norm_2(centre_of_cell) > M_DOMAIN_WIDTH)
                {   
                    // Delete cell and store it as a ghost node
                    rCellPopulation.RemoveCellUsingLocationIndex(location_index, (*cell_iter));
                    // Update vector of cells
                    cell_iter = rCellPopulation.rGetCells().erase(cell_iter);
 
                    // Change to chost node            
                    ghost_node_indices.insert(location_index);
                }
                else
                {
                    ++cell_iter;
                }
            }
            dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>* >(&rCellPopulation)->SetGhostNodes(ghost_node_indices);
            //dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>* >(&rCellPopulation)->RemoveDeadCells();
            //dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>* >(&rCellPopulation)->Update();
        }
        else
        {
            for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
                    cell_iter != rCellPopulation.End();
                    ++cell_iter)
            {
                // Get the coordinates of this cell centre
                c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                // if ((fabs(y-x)>M_DOMAINWIDTH) && (x<M_DOMAIN_X_MIN) && (x>M_DOMAIN_X_MAX) && (y<M_DOMAIN_Y_MIN) && (y>M_DOMAIN_Y_MAX))
                if ( norm_2(centre_of_cell) > M_DOMAIN_WIDTH)
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

    /*
     * This is a helper method to generate cells and is used in all simulations.
     */ 

    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells, double equilibriumVolume, double quiescentVolumeFraction, 
            double stemCellDivisionProbability, double stemCellMinimumDivisionAge)
    {

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_transit_cell_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_stem_cell_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        for (unsigned i=0; i<num_cells; i++)
        {
            BernoulliTrialWithContactInhibitionCellCycleModel* p_model = new BernoulliTrialWithContactInhibitionCellCycleModel();
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

            rCells.push_back(p_cell);
        }
    }

public:

    void TestGrowingMonolayers()
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

            /* 
            * == OS ==
            *
            * Simulate a growing monolayer using the
            * Overlapping Spheres model.
            */
            {
                std::cout << "Node Default \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Node/Default/Run_" +  out.str();

                /* 
                * == Default cut-off ==
                */
                // Create a simple mesh
                HoneycombMeshGenerator generator(2.0*M_INITIAL_WIDTH, 3.0*M_INITIAL_LENGTH,0);
                boost::shared_ptr<MutableMesh<2, 2> > p_generating_mesh = generator.GetMesh();

                p_generating_mesh->Translate(-M_INITIAL_WIDTH,-0.75*sqrt(3.0)*M_INITIAL_LENGTH);

                double cut_off_length = 1.5; //this is the default

                // Convert this to a NodesOnlyMesh
                NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
                p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0);

                // Create cells
                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumNodes(), cells, M_PI*0.25,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: M_PI*0.25 as r=0.5
                /* If want to make default node the same size as mesh, need to set the mature volume to around
                * 0.92, r\approx 0.52.
                */

                // Create a node-based cell population
                NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
                cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAncestorWriter>();

                cell_population.AddCellWriter<BoundaryCellWriter>();

                RoundOutCellPopulation(cell_population);

                // Create simulation from cell population
                OffLatticeSimulation<2> simulator(cell_population);
                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME); //50
                simulator.SetOutputDirectory(output_directory);
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

                MAKE_PTR(CicularityCalcModifier<2>, circularity_modifier);
                circularity_modifier->SetOutputDirectory(output_directory);
                simulator.AddSimulationModifier(circularity_modifier);

                // Mark Ancestors
                simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

                // Run simulation
                simulator.Solve();

                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

                // Clear memory
                delete p_mesh;
            }

            /*
            * == Larger cut-off ==
            * Cut-off = 2.0
            */
            {
                std::cout << "Node Large Cutoff \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Node/LargeCutoff/Run_" +  out.str();

                // Create a simple mesh
                HoneycombMeshGenerator generator(2.0*M_INITIAL_WIDTH, 3.0*M_INITIAL_LENGTH,0);
                boost::shared_ptr<MutableMesh<2, 2> > p_generating_mesh = generator.GetMesh();

                p_generating_mesh->Translate(-M_INITIAL_WIDTH,-0.75*sqrt(3.0)*M_INITIAL_LENGTH);

                double cut_off_length = 2.0; //this is larger than the default

                // Convert this to a NodesOnlyMesh
                NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
                p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0);

                // Create cells
                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumNodes(), cells, M_PI*0.25,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: M_PI*0.25 as r=0.5

                // Create a node-based cell population
                NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
                cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAncestorWriter>();

                cell_population.AddCellWriter<BoundaryCellWriter>();

                RoundOutCellPopulation(cell_population);

                // Create simulation from cell population
                OffLatticeSimulation<2> simulator(cell_population);
                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME); //50
                simulator.SetOutputDirectory(output_directory);
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

                MAKE_PTR(CicularityCalcModifier<2>, circularity_modifier);
                circularity_modifier->SetOutputDirectory(output_directory);
                simulator.AddSimulationModifier(circularity_modifier);

                // Mark Ancestors
                simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

                // Run simulation
                simulator.Solve();

                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

                // Clear memory
                delete p_mesh;
            }

            /*
            * == Smaller cut-off ==
            * Cut-off = 1.0
            */
            {
                std::cout << "Node Small Cutoff \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Node/SmallCutoff/Run_" +  out.str();

                // Create a simple mesh
                HoneycombMeshGenerator generator(2.0*M_INITIAL_WIDTH, 3.0*M_INITIAL_LENGTH,0);
                boost::shared_ptr<MutableMesh<2, 2> > p_generating_mesh = generator.GetMesh();

                p_generating_mesh->Translate(-M_INITIAL_WIDTH,-0.75*sqrt(3.0)*M_INITIAL_LENGTH);

                double cut_off_length = 1.0; //this is smaller than the default

                // Convert this to a NodesOnlyMesh
                NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
                p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 2.0);

                // Create cells
                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumNodes(), cells, M_PI*0.25,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: M_PI*0.25 as r=0.5

                // Create a node-based cell population
                NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
                cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAncestorWriter>();

                cell_population.AddCellWriter<BoundaryCellWriter>();

                RoundOutCellPopulation(cell_population);

                // Create simulation from cell population
                OffLatticeSimulation<2> simulator(cell_population);
                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME); //50
                simulator.SetOutputDirectory(output_directory);
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

                MAKE_PTR(CicularityCalcModifier<2>, circularity_modifier);
                circularity_modifier->SetOutputDirectory(output_directory);
                simulator.AddSimulationModifier(circularity_modifier);

                // Mark Ancestors
                simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

                // Run simulation
                simulator.Solve();

                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

                // Clear memory
                delete p_mesh;
            }

            /* 
            * == VT ==
            *
            * Simulate a growing monolayer using the
            * Voronoi Tessellation model.
            */
            {
                std::cout << "Ghost VT \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Mesh/Ghosts/Run_" +  out.str();

                /* 
                * == Ghosts == 
                */
                // Create mesh

                unsigned thickness_of_ghost_layer = 8;
                HoneycombMeshGenerator generator(2.0*M_INITIAL_WIDTH, 3.0*M_INITIAL_LENGTH, thickness_of_ghost_layer);
                boost::shared_ptr<MutableMesh<2, 2> > p_mesh = generator.GetMesh();

                p_mesh->Translate(-M_INITIAL_WIDTH,-0.75*sqrt(3.0)*M_INITIAL_LENGTH);

                // Get location indices corresponding to real cells
                std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

                // Create cells
                std::vector<CellPtr> cells;
                GenerateCells(location_indices.size(),cells,sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: sqrt(3.0)/2.0
                // To make same size as node default, decrease mature volume to around 0.8

                // Create tissue
                MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
                // MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, std::vector<unsigned>(), false, 15 , 1, 1);

                // Output Voroni for visualisation
                cell_population.AddPopulationWriter<VoronoiDataWriter>();
                cell_population.SetWriteVtkAsPoints(true);

                cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAncestorWriter>();

                cell_population.AddCellWriter<BoundaryCellWriter>();

                RoundOutCellPopulation(cell_population);

                // Create simulation from cell population
                OffLatticeSimulation<2> simulator(cell_population);
                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME); //50
                simulator.SetOutputDirectory(output_directory);
                simulator.SetOutputDivisionLocations(true);
                simulator.SetOutputCellVelocities(true);

                // Add volume tracking modifier
                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Create a force law and pass it to the simulation
                MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
                p_linear_force->SetMeinekeSpringStiffness(50.0);
                simulator.AddForce(p_linear_force);

                MAKE_PTR(CicularityCalcModifier<2>, circularity_modifier);
                circularity_modifier->SetOutputDirectory(output_directory);
                simulator.AddSimulationModifier(circularity_modifier);

                MAKE_PTR(GhostNodeRemovalModifier<2>, GNR_modifier);
                simulator.AddSimulationModifier(GNR_modifier);   

                // Mark Ancestors
                simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

                // Run simulation
                simulator.Solve();

                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
            }

            /* 
            * == No Ghosts infinite VT == 
            */
            // Create mesh
            {
                std::cout << "Infinite VT \n" << std::flush;

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Mesh/NoGhosts/Infinite/Run_" +  out.str();

                HoneycombMeshGenerator generator(2.0*M_INITIAL_WIDTH, 3.0*M_INITIAL_LENGTH, 0);
                boost::shared_ptr<MutableMesh<2, 2> > p_mesh = generator.GetMesh();

                p_mesh->Translate(-M_INITIAL_WIDTH,-0.75*sqrt(3.0)*M_INITIAL_LENGTH);

                // Create cells
                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumNodes(), cells,sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: sqrt(3.0)/2.0
                // Note: Tissue shrinks without proliferation due to boundary effects

                // Create tissue
                MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

                // Output Voroni for visualisation
                cell_population.AddPopulationWriter<VoronoiDataWriter>();
                cell_population.SetWriteVtkAsPoints(true);

                cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAncestorWriter>();

                cell_population.AddCellWriter<BoundaryCellWriter>();

                RoundOutCellPopulation(cell_population);

                // Create simulation from cell population
                OffLatticeSimulation<2> simulator(cell_population);
                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME); //50
                simulator.SetOutputDirectory(output_directory);
                simulator.SetOutputDivisionLocations(true);
                simulator.SetOutputCellVelocities(true);

                // Add volume tracking modifier
                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Create a force law and pass it to the simulation
                MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
                p_linear_force->SetMeinekeSpringStiffness(50.0);
                simulator.AddForce(p_linear_force);


                MAKE_PTR(CicularityCalcModifier<2>, circularity_modifier);
                circularity_modifier->SetOutputDirectory(output_directory);
                simulator.AddSimulationModifier(circularity_modifier);

                // Mark Ancestors
                simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

                // Run simulation
                simulator.Solve();

                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
            }

            /* 
            * == No Ghosts finite VT == 
            */
            // Create mesh
            {
                std::cout << "Finite VT \n" << std::flush;

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Mesh/NoGhosts/Finite/Run_" +  out.str();

                
                HoneycombMeshGenerator generator(2.0*M_INITIAL_WIDTH, 3.0*M_INITIAL_LENGTH);
                boost::shared_ptr<MutableMesh<2, 2> > p_mesh = generator.GetMesh();

                p_mesh->Translate(-M_INITIAL_WIDTH,-0.75*sqrt(3.0)*M_INITIAL_LENGTH);

                // Get location indices corresponding to real cells
                std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

                // Create cells
                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumNodes(), cells,sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); // mature volume: sqrt(3.0)/2.0
                // Note: Tissue shrinks without proliferation due to boundary effects

                // Create tissue
                MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
                cell_population.SetBoundVoronoiTessellation(true);
                cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAncestorWriter>();

                cell_population.AddCellWriter<BoundaryCellWriter>();

                RoundOutCellPopulation(cell_population);

                // Create simulation from cell population
                OffLatticeSimulation<2> simulator(cell_population);
                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME); //50
                simulator.SetOutputDirectory(output_directory);
                simulator.SetOutputDivisionLocations(true);
                simulator.SetOutputCellVelocities(true);

                // Add volume tracking modifier
                MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                simulator.AddSimulationModifier(p_modifier);

                // Create a force law and pass it to the simulation
                MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
                p_linear_force->SetMeinekeSpringStiffness(50.0);
                simulator.AddForce(p_linear_force);

                MAKE_PTR(CicularityCalcModifier<2>, circularity_modifier);
                circularity_modifier->SetOutputDirectory(output_directory);
                simulator.AddSimulationModifier(circularity_modifier);

                // Mark Ancestors
                simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

                // Run simulation
                simulator.Solve();

                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
            }

            /* 
            * == VM ==
            *
            * Simulate a growing monolayer using the
            * Cell Vertex model.
            */
            {
                std::cout << "Vertex Jagged \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Vertex/Jagged/Run_" +  out.str();

                /* 
                * == Jagged edge == 
                */
                // Create cells
                HoneycombVertexMeshGenerator mesh_generator(2.0*M_INITIAL_WIDTH, 3.0*M_INITIAL_LENGTH);
                boost::shared_ptr<MutableVertexMesh<2, 2> > p_mesh = mesh_generator.GetMesh();
                // p_mesh->SetCellRearrangementThreshold(0.1);

                p_mesh->Translate(-M_INITIAL_WIDTH,-0.75*sqrt(3.0)*M_INITIAL_LENGTH+0.25);

                // Create cells
                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumElements(), cells,sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); //mature_volume = 1.0

                // Create tissue
                VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
                cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAncestorWriter>();

                cell_population.AddCellWriter<BoundaryCellWriter>();

                RoundOutCellPopulation(cell_population);

                // Create crypt simulation from cell population
                OffLatticeSimulation<2> simulator(cell_population);
                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME); //50
                simulator.SetOutputDirectory(output_directory);
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
                p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(2.0);
                simulator.AddForce(p_force);

                // Mark Ancestors
                simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

                MAKE_PTR(CicularityCalcModifier<2>, circularity_modifier);
                circularity_modifier->SetOutputDirectory(output_directory);
                simulator.AddSimulationModifier(circularity_modifier);

                // Add target area modifier
                // MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
                // p_growth_modifier->SetGrowthDuration(0);
                // simulator.AddSimulationModifier(p_growth_modifier);

                // Run simulation
                simulator.Solve();

                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

                // Clear singletons
                Warnings::Instance()->QuietDestroy();
            }

            /*
            * == Smooth edge == 
            */
            {
                std::cout << "Vertex Smooth \n" << std::flush; 

                p_gen->Reseed(50);
                std::string output_directory = M_HEAD_FOLDER + "/Vertex/Smooth/Run_" +  out.str();

                // Create cells
                HoneycombVertexMeshGenerator mesh_generator(2.0*M_INITIAL_WIDTH, 3.0*M_INITIAL_LENGTH);
                boost::shared_ptr<MutableVertexMesh<2, 2> > p_mesh = mesh_generator.GetMesh();
                // p_mesh->SetCellRearrangementThreshold(0.1);

                p_mesh->Translate(-M_INITIAL_WIDTH,-0.75*sqrt(3.0)*M_INITIAL_LENGTH+0.25);

                // Create cells
                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumElements(), cells,sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); //mature_volume = 1.0

                // Create tissue
                VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
                SmoothVertexMeshEdges(cell_population);
                cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAncestorWriter>();

                cell_population.AddCellWriter<BoundaryCellWriter>();


                RoundOutCellPopulation(cell_population);

                // Create crypt simulation from cell population
                OffLatticeSimulation<2> simulator(cell_population);
                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME); //50
                simulator.SetOutputDirectory(output_directory);
                simulator.SetOutputDivisionLocations(true);
                simulator.SetOutputCellVelocities(true);

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
                p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(2.0);
                simulator.AddForce(p_force);

                MAKE_PTR(CicularityCalcModifier<2>, circularity_modifier);
                circularity_modifier->SetOutputDirectory(output_directory);
                simulator.AddSimulationModifier(circularity_modifier);

                // Mark Ancestors
                simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

                // Add target area modifier
                // MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
                // p_growth_modifier->SetGrowthDuration(0);
                // simulator.AddSimulationModifier(p_growth_modifier);

                // Run simulation
                simulator.Solve();

                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

                // Clear singletons
                Warnings::Instance()->QuietDestroy();
            }

            /*
            * == Curved edge == 
            */
            {
                std::cout << "Vertex Curved \n" << std::flush; 

                p_gen->Reseed(sim_index);
                std::string output_directory = M_HEAD_FOLDER + "/Vertex/Curved/Run_" +  out.str();

                
                // Create cells
                HoneycombVertexMeshGenerator mesh_generator(2.0*M_INITIAL_WIDTH, 3.0*M_INITIAL_LENGTH);
                boost::shared_ptr<MutableVertexMesh<2, 2> > p_mesh = mesh_generator.GetMesh();
                // p_mesh->SetCellRearrangementThreshold(0.1);

                p_mesh->Translate(-M_INITIAL_WIDTH,-0.75*sqrt(3.0)*M_INITIAL_LENGTH+0.25);

                // Create cells
                std::vector<CellPtr> cells;
                GenerateCells(p_mesh->GetNumElements(), cells,sqrt(3.0)/2.0,M_CONTACT_INHIBITION_LEVEL,
                        M_STEM_CELL_DIVISION_PROBABILITY,M_STEM_CELL_MINIMUM_DIVISION_AGE); //mature_volume = 1.0

                // Create tissue
                VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
                cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellAncestorWriter>();

                cell_population.AddCellWriter<BoundaryCellWriter>();

                RoundOutCellPopulation(cell_population);

                // Create crypt simulation from cell population
                OffLatticeSimulation<2> simulator(cell_population);
                simulator.SetDt(M_DT_TIME);
                simulator.SetSamplingTimestepMultiple(M_SAMPLE_TIME);
                simulator.SetEndTime(M_END_TIME); //50
                simulator.SetOutputDirectory(output_directory);
                simulator.SetOutputDivisionLocations(true);
                simulator.SetOutputCellVelocities(true);

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
                p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(2.0);
                simulator.AddForce(p_force);

                MAKE_PTR(CicularityCalcModifier<2>, circularity_modifier);
                circularity_modifier->SetOutputDirectory(output_directory);
                simulator.AddSimulationModifier(circularity_modifier);

                // Mark Ancestors
                simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

                // Add target area modifier
                // MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
                // p_growth_modifier->SetGrowthDuration(0);
                // simulator.AddSimulationModifier(p_growth_modifier);

                // Run simulation
                simulator.Solve();

                // Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);

                // Clear singletons
                Warnings::Instance()->QuietDestroy();
            }
        }
    }
};

#endif /* TESTGROWINGMONOLAYER_HPP_ */
