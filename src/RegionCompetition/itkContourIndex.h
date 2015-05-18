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
 * ContourIndex.h
 *
 *  Created on: Nov 16, 2009
 *      Author: janickc
 */

#ifndef ITKCONTOURINDEX_H_
#define ITKCONTOURINDEX_H_

#include "itkIndex.h"

namespace itk {

    template <unsigned int VDimension = 2 >
            class ContourIndex : public itk::Index < VDimension > {
        //    typedef ContourIndex Self;

    public:
        static const unsigned int IndexDimension = VDimension;
    };

    template <unsigned int VDimension>
    bool operator<(
            const ContourIndex<VDimension>& left,
            const ContourIndex<VDimension>& right) {

        for (unsigned int i = 0; i < VDimension; i++) {
            if (left[i] < right[i]) return true;
            if (left[i] == right[i]) continue;
            return false;
        }
        return false;
    }




    /**
     * Returns the squared distance to another ContourIndex.
     */
    //    template <unsigned int VDimension>
    //    float distTo2(const ContourIndex<VDimension>& aCIndex) {
    //        float vDist2 = 0.0f;
    //        for(int vD = 0; vD < VDimension; vD++) {
    //            vDist2 += (m_Index[vD] - aCIndex.m_Index[vD]) * (m_Index[vD] - aCIndex.m_Index[vD]);
    //        }
    //    }

    //    bool operator<(const ContourIndex<3>& left, const ContourIndex<3>& right) {
    //    	if(left.x < right.x) return true;
    //    	if(left.x == right.x){
    //    		if(left.y < right.y) return true;
    //    		if(left.y == right.y) {
    //    			return left.z < right.z;
    //    		}
    //    	}
    //    	return false;
    //    }

    //bool operator==(const ContourIndex<2>& left, const ContourIndex<2>& right) {
    //	return (left.x == right.x && left.y == right.y);
    //}
    //
    //bool operator==(const ContourIndex<3>& left, const ContourIndex<3>& right) {
    //	return (left.x == right.x && left.y == right.y && left.z == right.z);
    //}

    //} // end namespace itk

    //namespace Functor {

    template <unsigned int VDimension>
    struct ContourIndexHashFunctor {

        size_t operator()(const itk::ContourIndex<VDimension> & aKey) const {
            // here is the computation of the hash function
            int vVal = 0;
            unsigned int nbits = 32 / VDimension;
            for (unsigned int vI = 0; vI < VDimension; vI++) {
                if (vI != 0) {
                    vVal = vVal << nbits;
                }
                vVal = vVal | (aKey[vI] & (unsigned int) pow(2, nbits));
            }
            return vVal;
        }

        bool operator()(
                const itk::ContourIndex<VDimension>& aLeft,
                const itk::ContourIndex<VDimension>& aRight) const {
            //here should be the code to compare two keys
            for (unsigned int vI = 0; vI < VDimension; vI++) {
                if (aLeft[vI] != aRight[vI]) return false;
            }
            return true;
        }
    };

    template <>
    struct ContourIndexHashFunctor < 2 > {

        size_t operator()(const itk::ContourIndex < 2 > & aKey) const {
            // here is the computation of the hash function
            return (int) (aKey[0] << 16 | ((aKey[1] & 0xFFFF)));
        }

        bool operator()(
                const itk::ContourIndex < 2 > & aLeft,
                const itk::ContourIndex < 2 > & aRight) const {
            // compare two keys
            return (aLeft[0] == aRight[0]) && (aLeft[1] == aRight[1]);
        }
    };

    template <>
    struct ContourIndexHashFunctor < 3 > {

        size_t operator()(const itk::ContourIndex < 3 > & aKey) const {
            // here is the computation of the hash function
            return (int) (aKey[2] << 22 |
                    ((aKey[1] & 0x07FF) << 11) |
                    (aKey[0] & 0x07FF));
        }

        bool operator()(
                const itk::ContourIndex < 3 > & aLeft,
                const itk::ContourIndex < 3 > & aRight) const {
            // compare two keys
            return (aLeft[0] == aRight[0]) &&
                    (aLeft[1] == aRight[1]) &&
                    (aLeft[2] == aRight[2]);
        }
    };
}

#endif /* ITKCONTOURINDEX_H_ */
