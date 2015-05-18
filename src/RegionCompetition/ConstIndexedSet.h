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
 * File:   ConstIndexedSet.h
 * Author: janickc
 *
 * Created on February 21, 2012, 3:19 PM
 */

#ifndef CONSTINDEXEDSET_H
#define	CONSTINDEXEDSET_H

#include <map> // TODO: make a hashed version
#include <vector>
#include <utility>

template <class TElement>
class ConstIndexedSet {

    typedef TElement ElementType;

    typedef std::vector<const ElementType*> PtrVecType;

    typedef unsigned int IndexingType;
//    typedef PtrVecType::size_type IndexingType;

    typedef std::map<ElementType, IndexingType> MapType;

    typedef typename MapType::value_type MapValueType;
    typedef typename MapType::iterator MapIteratorType;

    MapType m_map;
    PtrVecType m_ptr_vec;

public:

    ConstIndexedSet() {
    };

    unsigned long size() const {
        return m_map.size();
    }

    IndexingType insert(const ElementType &aE) {

        MapIteratorType it = m_map.find(aE);
        if (it == m_map.end()) {
            it = m_map.insert(MapValueType(aE, m_map.size())).first;
            m_ptr_vec.push_back(&it->first); 
        }
        return it->second; // return the index
    }

    const ElementType& elementAt(unsigned int aIndex) {
        return *(m_ptr_vec[aIndex]);
    }
    
    ElementType& operator[](const unsigned int aIndex) {
        return *(m_ptr_vec[aIndex]);
    }
    
    
    /**
     * @param aE Element to erase from the container.
     * @return true if an element was deleted; false if there has not been such
     * an element in the container.
     */
    bool erase(const ElementType &aE) {
        MapIteratorType it = m_map.find(aE);

        if (it == m_map.end()) {
            return false;
        } else {
            /// Get the index of the element to delete:
            IndexingType vIndex = it->second;

            /// We move the last element:
            m_map[*m_ptr_vec.back()] = vIndex;

            /// Delete it in the map
            m_map.erase(it);
            
            /// Update the vector: move the last element to the element to 
            /// delete and remove the last element.
            m_ptr_vec[vIndex] = m_ptr_vec.back();
            m_ptr_vec.pop_back();

        }
        return true;
    }
};
#endif	/* CONSTINDEXEDSET_H */

