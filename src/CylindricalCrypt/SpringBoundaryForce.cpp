#include "SpringBoundaryForce.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "MutableVertexMesh.hpp"

template<unsigned DIM>
SpringBoundaryForce<DIM>::SpringBoundaryForce()
    : AbstractForce<DIM>(),
	  mMinYvalue(0.0)
{
}

template<unsigned DIM>
SpringBoundaryForce<DIM>::~SpringBoundaryForce()
{
}

template<unsigned DIM>
void SpringBoundaryForce<DIM>::SetMinimumYValue(double minYvalue)
{
    mMinYvalue = minYvalue;
}

template<unsigned DIM>
double SpringBoundaryForce<DIM>::GetMinimumYValue()
{
    return mMinYvalue;
}

template<unsigned DIM>
void SpringBoundaryForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    double dt = SimulationTime::Instance()->GetTimeStep();
    double turn_on_spring_distance = mMinYvalue + 1.0;

    //  (dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation))
    if(bool(dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&rCellPopulation))) 
    {
        // Iterate over the nodes
        for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End();
                 ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<DIM>* p_node = rCellPopulation.GetNode(node_index);
            c_vector<double, DIM> node_location = p_node->rGetLocation();

            c_vector<double, DIM> force_contribution;
            force_contribution[0] = 0.0;
            force_contribution[1] = 0.0;

            if(node_location[1] < turn_on_spring_distance)
            {
                force_contribution[1] = 10.0*(turn_on_spring_distance - node_location[1]);            
            }

            p_node->AddAppliedForceContribution(force_contribution);

        }

    }
    else if(bool(dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation))) 
    {
        VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        MutableVertexMesh<DIM,DIM>* p_mesh = static_cast<MutableVertexMesh<DIM,DIM>*>(&(p_cell_population->rGetMesh()));

        // Iterate over the nodes
        for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End();
                 ++cell_iter)
        {
            c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

            c_vector<double, DIM> force_contribution;
            force_contribution[0] = 0.0;
            force_contribution[1] = 0.0;

            if(real_node_location[1] < turn_on_spring_distance)
            {
                force_contribution[1] = 0.1*(turn_on_spring_distance - real_node_location[1]);            
            }

            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(*cell_iter);
            for (unsigned i=0; i<p_element->GetNumNodes(); i++)
            {
                unsigned node_global_index = p_element->GetNodeGlobalIndex(i);
                Node<DIM>* p_node_a = p_mesh->GetNode(node_global_index);

                p_node_a->AddAppliedForceContribution(force_contribution);   
            }

        }

    }
        
    
}

template<unsigned DIM>
void SpringBoundaryForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<mMinYvalue>" << mMinYvalue << "</mMinYvalue> \n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class SpringBoundaryForce<1>;
template class SpringBoundaryForce<2>;
template class SpringBoundaryForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SpringBoundaryForce)
