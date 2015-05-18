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

#ifndef CONSTINDEXEDHASHSET_H
#define	CONSTINDEXEDHASHSET_H

#include "itksys/hash_map.hxx"
#include <vector>
#include <utility>

template <class TElement, class THashFcn, class TEqualKey>
class ConstIndexedHashSet {

    typedef ConstIndexedHashSet<TElement, THashFcn, TEqualKey> SelfType;
    
    typedef TElement ElementType;

    typedef std::vector<const ElementType*> PtrVecType;

    typedef unsigned int IndexingType;
//    typedef PtrVecType::size_type IndexingType;

    typedef THashFcn HashFunctionType;
    typedef TEqualKey EqualKeyType;
    typedef itksys::hash_map<ElementType, IndexingType, HashFunctionType, EqualKeyType> MapType;

    typedef typename MapType::value_type MapValueType;
    typedef typename MapType::iterator MapIteratorType;

    MapType m_map;
    PtrVecType m_ptr_vec;
    ElementType m_last_del_element;

public:

    ConstIndexedHashSet() {
    };

    const unsigned long size() const {
        return m_map.size();
    }

    /**
     * Finds an element in (amortized) constant time.
     * @param aE Element to find.
     * @return The index of the element. size() is returned if the element is 
     * not part of the set. 
     */
    IndexingType find(const ElementType &aE) {
        MapIteratorType it = m_map.find(aE);
        if (it == m_map.end()) {
            return m_map.size();
        } 
        return it->second;
    }
    
    bool contains(const ElementType &aE) {
        return !(find(aE) == m_map.size());
    }
    
    /**
     * Can be used also to replace elements.
     * @param aE
     * @return 
     */
    IndexingType insert(const ElementType &aE) {

        MapIteratorType it = m_map.find(aE);
        if (it == m_map.end()) {
            it = m_map.insert(MapValueType(aE, m_map.size())).first;
            m_ptr_vec.push_back(&it->first);
        } else {
            /// replace the entry  
            IndexingType vec_i = it->second;
            m_last_del_element = it->first;
            m_map.erase(it);
            it = m_map.insert(MapValueType(aE, vec_i)).first;
            m_ptr_vec[vec_i] = &it->first;
        }
        return it->second; // return the index
    }
    
    void join(const SelfType &aSet) {
        for(unsigned int vI = 0; vI < aSet.size; vI++) {
            this->insert(aSet[vI]);
        }
    }
    
    const ElementType& elementAt(unsigned int aIndex) const {
        return *(m_ptr_vec[aIndex]);
    }
    
    const ElementType& operator[](const unsigned int aIndex) const {
        return *(m_ptr_vec[aIndex]);
    }

    bool erase(const ElementType &aE) {
        MapIteratorType it = m_map.find(aE);
        if (it == m_map.end()) {
            return false;
        } else {
            /// Get the index of the element to delete:
            IndexingType vIndex = it->second;
            m_last_del_element = it->first;
            
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
    
    const ElementType& getLastDeletedElement() const {
        return m_last_del_element;
    }
};

#endif	/* CONSTINDEXEDHASHSET_H */

