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
 * ContourContainer.h
 *
 *  Created on: Nov 16, 2009
 *      Author: janickc
 */

#ifndef ITKCONTOURCONTAINER_H_
#define ITKCONTOURCONTAINER_H_ 

#include "itkContourIndex.h"
#include <map>

namespace itk {

    template <class TIndex, class TValue>
    class ContourContainer : public std::map<TIndex, TValue > {
    public:
        typedef ContourContainer<TIndex, TValue> Self;
        typedef std::map<TIndex, TValue > SuperclassType;
        typedef typename SuperclassType::iterator IteratorType; //TODO test Self::iterator IteratorType;
        typedef TIndex IndexType;
        typedef TValue ValueType;

//        void RemoveSet(Self* aSetToRemove) {
//            //remove points from the container
//            IteratorType vRemoveIterator = aSetToRemove->begin();
//            IteratorType vRemoveMapEnd = aSetToRemove->end();
//            for (; vRemoveIterator != vRemoveMapEnd; ++vRemoveIterator) {
//                this->erase(vRemoveIterator->first); //remove by key
//            }
//        };
    };

} //end namespace itk

#endif /* CONTOURCONTAINER_H_ */
