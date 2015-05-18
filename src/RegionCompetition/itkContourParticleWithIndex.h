/*=========================================================================
 *
 *  Copyright MOSAIC Group, ETHZ and MPI-CBG Dresden
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

/*=========================================================================
 * This file is part of the Region Competition segmentation algorithm 
 * described in:
 *
 * Cardinale J, Paul G, Sbalzarini I (2012), "Discrete region com-
 * petition for unknown numbers of connected regions". 
 * Image Processing, IEEE Transactions on, 21(8):3531--3545
 *
 * In order to ensure financial support and allow 
 * further development of this software, please cite above publication in
 * all your documents and manuscripts that made use of this software. 
 * Thanks a lot!
 *
 * Authors:
 * Janick Cardinale
 *=========================================================================*/
/*
 * ContourPointWithIndex.h
 *
 */

#ifndef ITKCONTOURPARTICLEWITHINDEX_H_
#define ITKCONTOURPARTICLEWITHINDEX_H_

#include "itkContourIndex.h"
#include "itkContourParticle.h"

namespace itk {

    template <class TContourIndex, class TContourParticle >
    class ContourParticleWithIndex {
    public:
        typedef TContourParticle ContourParticleType;
        typedef TContourIndex ContourIndexType;
        ContourIndexType m_ContourIndex;
        ContourParticleType m_ContourParticle;

        inline ContourParticleWithIndex(TContourIndex aIndex, TContourParticle aPoint)
        : m_ContourIndex(aIndex), m_ContourParticle(aPoint) {
        };
    private:

    };

    template <class TContourIndex, class TContourParticle >
    bool operator<(const ContourParticleWithIndex<TContourIndex, TContourParticle>& left,
    const ContourParticleWithIndex<TContourIndex, TContourParticle>& right) {
        return left.m_ContourParticle.m_energyDifference < right.m_ContourParticle.m_energyDifference;
    }




}

#endif /* ITKCONTOURPARTICLEWITHINDEX_H_ */
