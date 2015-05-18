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
 * File:   itkMinimalParticle.h
 * Author: janickc
 *
 * Created on January 24, 2012, 4:37 PM
 */

#ifndef ITKMINIMALPARTICLE_H
#define	ITKMINIMALPARTICLE_H

#include "itkIndex.h"

namespace itk {

    template <unsigned int VDimension = 2, typename TLabel = unsigned int >
            class MinimalParticle {
    public:
        typedef MinimalParticle<VDimension> Self;
        typedef Index<VDimension> IndexType;
        typedef TLabel LabelType;

        IndexType m_Index;
//        LabelType m_Label;
        LabelType m_CandidateLabel;
        float m_Proposal;

        MinimalParticle(){}; // default constructor

//        inline MinimalParticle(IndexType aIndex, LabelType aLabel, LabelType aCandLabel)//, unsigned int aTopoFG)
//        : m_Index(aIndex), m_Label(aLabel), m_CandidateLabel(aCandLabel) {//, m_TopoFG(aTopoFG) {
//        };
        inline MinimalParticle(IndexType aIndex,  LabelType aCandLabel, float aProposal)
        : m_Index(aIndex),  m_CandidateLabel(aCandLabel), m_Proposal(aProposal) {
        };
    private:
    };

    template <unsigned int VDimension, typename TLabel>
    bool operator<(
            const MinimalParticle<VDimension, TLabel>& left,
            const MinimalParticle<VDimension, TLabel>& right) {

        for (unsigned int i = 0; i < VDimension; i++) {
            if (left.m_Index[i] < right.m_Index[i]) return true;
            else if (left.m_Index[i] == right.m_Index[i]) continue;
            else return false;
        }

//        if (left.m_Label < right.m_Label) return true;
//        else if (left.m_Label > right.m_Label) return false;

        if (left.m_CandidateLabel < right.m_CandidateLabel) return true;
        /// in all other cases we can return false 
        return false;
    }

    template <unsigned int VDimension>
    struct MinimalParticleHashFunctor {

        size_t operator()(const itk::MinimalParticle<VDimension> & aKey) const {
            // Computation of the hash function:
            int vVal = 0;
            unsigned int nbits = 32 / VDimension;
            for (unsigned int vI = 0; vI < VDimension; vI++) {
                if (vI != 0) {
                    vVal = vVal << nbits;
                }
                vVal = vVal | (aKey.m_Index[vI] & (unsigned int) pow(2, nbits));
            }
            /// we just add up and hope:
//            vVal += 27 * aKey.m_Label + 31 * aKey.m_CandidateLabel;
            vVal += 31 * aKey.m_CandidateLabel;
            
            return vVal;
        }

        bool operator()(
                const itk::MinimalParticle<VDimension>& aLeft,
                const itk::MinimalParticle<VDimension>& aRight) const {
            //here should be the code to compare two keys
            for (unsigned int vI = 0; vI < VDimension; vI++) {
                if (aLeft[vI] != aRight[vI]) return false;
            }
//            if(aLeft.m_Label != aRight.m_Label) return false;
            if(aLeft.m_CandidateLabel != aRight.m_CandidateLabel) return false;
            return true;
        }
    };

    template <>
    struct MinimalParticleHashFunctor < 2 > {

        size_t operator()(const itk::MinimalParticle < 2 > & aKey) const {
            // here is the computation of the hash function
            return (int) (aKey.m_Index[0] << 15 |
                    ((aKey.m_Index[1] & 0x7FFF) << 2) |
//                    ((2 * aKey.m_Label + 3 * aKey.m_CandidateLabel) & 0x0003));
                    ((aKey.m_CandidateLabel) & 0x0003));
        }

        bool operator()(
                const itk::MinimalParticle < 2 > & aLeft,
                const itk::MinimalParticle < 2 > & aRight) const {
            // compare two keys
            return (aLeft.m_Index[0] == aRight.m_Index[0]) &&
                    (aLeft.m_Index[1] == aRight.m_Index[1]) &&
//                    (aLeft.m_Label == aRight.m_Label) &&
                    (aLeft.m_CandidateLabel == aRight.m_CandidateLabel);
        }
    };

    template <>
    struct MinimalParticleHashFunctor < 3 > {

        size_t operator()(const itk::MinimalParticle < 3 > & aKey) const {
            // here is the computation of the hash function
            return (int) (aKey.m_Index[2] << 22 |
                    ((aKey.m_Index[1] & 0x03FF) << 12) |
                    ((aKey.m_Index[0] & 0x03FF) << 2) |
//                    ((2 * aKey.m_Label + 3 * aKey.m_CandidateLabel) & 0x0003));
                    ((aKey.m_CandidateLabel) & 0x0003));
        }

        bool operator()(
                const itk::MinimalParticle < 3 > & aLeft,
                const itk::MinimalParticle < 3 > & aRight) const {
            // compare two keys
            return (aLeft.m_Index[0] == aRight.m_Index[0]) &&
                    (aLeft.m_Index[1] == aRight.m_Index[1]) &&
                    (aLeft.m_Index[2] == aRight.m_Index[2]) &&
//                    (aLeft.m_Label == aRight.m_Label) &&
                    (aLeft.m_CandidateLabel == aRight.m_CandidateLabel);
        }
    };
    
    template <unsigned int VDimension, typename TLabel >
    std::ostream & operator<<(std::ostream & os, const MinimalParticle< VDimension, TLabel > & ind)
    {
        os << "(";
//        os << ind.m_Index << ", " << ind.m_Label << ", " << ind.m_CandidateLabel;
        os << ind.m_Index << ", " << ind.m_CandidateLabel;
        os << ")";
        return os;
    }
} // end namespace

#endif	/* ITKMINIMALPARTICLE_H */

