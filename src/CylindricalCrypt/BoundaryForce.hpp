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

#ifndef BOUNDARYFORCE_HPP_
#define BOUNDARYFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

/**
 * A boundary force class for use in crypt simulations
 * to prevent cells flattening along the bottom of the crypt (y=0).
 *
 * The boundary force is taken to be proportional to 1/y^2
 * The constant of proportionality is given by the parameter mForceStrength.
 */
template<unsigned DIM>
class BoundaryForce  : public AbstractForce<DIM>
{

private:

    /** Parameters determining the strength of the force acting on nodes below y=CutOffHeight. */
    double mForceStrength;
    double mCutOffHeight;

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
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mForceStrength;
        archive & mCutOffHeight;
    }

public:

    /**
     * Constructor.
     */
    BoundaryForce();

    /**
     * Destructor.
     */
    ~BoundaryForce();

    /**
     * Strength of the repulsive force at y=1.
     * (Force goes to infinity as y goes to 0.)
     *
     * @param forceStrength the force strength
     */
    void SetForceStrength(double forceStrength);

    /**
     * Get the force strength.
     *
     * @return mForceStrength
     */
    double GetForceStrength();

    /**
     * Use a cutoff point, don't add boundary force contribution if 
     * a node is too far from the boundary
     *
     * @param cutOffHeight the cutoff to use
     */
    void SetCutOffHeight(double cutOffHeight);

    /**
     * Get the cutoff height.
     *
     * @return mCutOffHeight
     */
    double GetCutOffHeight();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the boundary force on each node in the cell population.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BoundaryForce)

#endif /*BOUNDARYFORCE_HPP_*/
