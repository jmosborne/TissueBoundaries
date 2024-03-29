/*

Copyright (c) 2005-2021, University of Oxford.
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

#ifndef JaggedVertexEdgesModifierSimplified_HPP_
#define JaggedVertexEdgesModifierSimplified_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MutableMesh.hpp"

/**
 * A modifier class which at each simulation time step removes nodes
 * that are only contained in one element.
 * 
 * Only for use with VertexBasedCellPopulation.
 */
template<unsigned DIM>
class JaggedVertexEdgesModifierSimplified : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
        archive & mMaxEdgeLength;
        archive & mMinEdgeLength;
    }
    /**
     * The longest an edge can be.
     */
    double mMaxEdgeLength;

    /**
     * The smallest an edge can be.
     */
    double mMinEdgeLength;

public:

    /**
     * Default constructor.
     */
    JaggedVertexEdgesModifierSimplified();

    /**
     * Destructor.
     */
    virtual ~JaggedVertexEdgesModifierSimplified();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    double DistanceToEdgeFromPoint( c_vector<double, DIM> P1,  c_vector<double, DIM> P2,  c_vector<double, DIM> P0);

    bool IsNodeABNeighbours( unsigned node_a,unsigned node_b, AbstractCellPopulation<DIM,DIM>& rCellPopulation );

    bool Remove3NodeVoid(bool performed_edge_modifier, unsigned node_index, Node<DIM>* p_neighbour_node_1, Node<DIM>* p_neighbour_node_2, VertexBasedCellPopulation<DIM>* p_cell_population, MutableVertexMesh<DIM,DIM>* p_mesh);

    bool Remove4NodeVoid(bool performed_edge_modifier, unsigned node_index, unsigned node_neighbour_1, unsigned node_neighbour_2, VertexBasedCellPopulation<DIM>* p_cell_population, MutableVertexMesh<DIM,DIM>* p_mesh);

    /**
     * Helper method to Do the remshing.
     *
     * @param rCellPopulation reference to the cell population
     */
    // bool SmoothEdges(AbstractCellPopulation<DIM,DIM>& rCellPopulation);
    bool CheckForVConfig(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here);
    bool CheckFor3NodeVoids(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here);
    bool CleanUpRogueNodes(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here);
    bool CheckForInternalNode(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here);
    bool CheckForNodeEdgeInt(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here);
    bool CheckT3Distrupt(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here);
    bool DirtyHack(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here);
    bool ClearFreeInternalNodes(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here);

    

    // bool RefineEdges(AbstractCellPopulation<DIM,DIM>& rCellPopulation);
    bool RefineEdges(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned numb_times_here);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(JaggedVertexEdgesModifierSimplified)

#endif /*JaggedVertexEdgesModifierSimplified_HPP_*/
