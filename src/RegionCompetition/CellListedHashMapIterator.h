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
#ifndef __CellListedHashMapIterator_H
#define __CellListedHashMapIterator_H


//#include <iostream>
#include <utility>
//#include <string>

template <class TCellListedHashMap>
class CellListedHashMapIterator {
public:
    typedef TCellListedHashMap CellListedHashMapType;
    typedef CellListedHashMapIterator<CellListedHashMapType> Self;

    typedef typename CellListedHashMapType::HashMapType HashMapType;
    typedef typename CellListedHashMapType::HashMapContainerType ContainerType;
    typedef typename CellListedHashMapType::IndexType CLHMKeyType;
    typedef typename CellListedHashMapType::ValueType CLHMValueType;
    typedef typename CellListedHashMapType::HashFunctionType CLHMHashFunctionType;
    typedef typename CellListedHashMapType::EqualKeyType CLHMEqualKeyType;

    typedef typename HashMapType::iterator HMIteratorType;
    typedef typename HashMapType::key_type KeyType;
    typedef typename HashMapType::data_type DataType;
    typedef typename HashMapType::value_type ValueType;

    template <
    unsigned int,
    class CLHMKeyType,
    class CLHMValueType,
    class CLHMHashFunctionType,
    class CLHMEqualKeyType> friend class CellListedHashMap;

private:
    ContainerType* m_container;

    unsigned int m_currentCell;
    HMIteratorType m_HMIt;

    unsigned int m_nbCells;

public:

    CellListedHashMapIterator(CellListedHashMapType* aCLHM, bool aBegin) {
        m_container = &(aCLHM->m_Container);
        m_nbCells = m_container->size();

        //    std::cout << "check for the end iterators:" << std::endl;
        // if( ((*aC)[0]).end() ==  ((*aC)[1]).end()) {
        //   std::cout << "the same end!" <<  std::endl;
        // }

        if (aBegin) {
            m_currentCell = 0;
            m_HMIt = ((*m_container)[m_currentCell]).begin();
//            std::cout << "size: " << ((*aC)[m_currentCell]).size() << std::endl;
            while (m_currentCell < m_nbCells-1 &&
                    (((*m_container)[m_currentCell]).begin() ==
                    ((*m_container)[m_currentCell]).end())) {
                m_currentCell++;
                m_HMIt = ((*m_container)[m_currentCell]).begin();
            }
        } else {
            m_currentCell = (*m_container).size() - 1;
            m_HMIt = ((*m_container)[m_nbCells - 1]).end();
        }
    };

    CellListedHashMapIterator(CellListedHashMapType* aCLHM,
            KeyType aIndex,
            unsigned int aCellIndex) {

        m_container = &(aCLHM->m_Container);
        m_nbCells = m_container->size();

        m_HMIt = (*m_container)[aCellIndex].find(aIndex);
        if (m_HMIt != (*m_container)[aCellIndex].end()) {
            m_currentCell = aCellIndex;
            //      std::cout << "found!" << std::endl;
        } else {
            //      std::cout << "not found!" << std::endl;
            m_currentCell = (*m_container).size() - 1;
            m_HMIt = ((*m_container)[m_nbCells - 1]).end();
        }
    };

    CellListedHashMapIterator & operator++() {
        ++m_HMIt;
        if (m_HMIt == (*m_container)[m_currentCell].end()) {
            m_currentCell++;

            //feed forward the empty cells
            while (m_currentCell < m_nbCells &&
                    (*m_container)[m_currentCell].begin() ==
                    (*m_container)[m_currentCell].end()) {
                m_currentCell++;
            }
            if (m_currentCell < m_nbCells) {
                // set the iterator to the first element in the new cell
                m_HMIt = (*m_container)[m_currentCell].begin();
            } else {
                // cleanup: set the counter to the last cell
                m_currentCell = m_nbCells - 1;
                // give back the end of the last cell
                m_HMIt = (*m_container)[m_currentCell].end();
            }
        }
        return *this;
    };

    std::pair<const KeyType, DataType> & operator*() {
        return *m_HMIt;
    };

    std::pair<const KeyType, DataType> * operator->() {
        return &(*m_HMIt);
    };

    bool operator!=(const CellListedHashMapIterator& aIt) {
        return (m_HMIt != aIt.m_HMIt || m_currentCell != aIt.m_currentCell);
    };

    bool operator==(const CellListedHashMapIterator& aIt) {
        return (m_HMIt == aIt.m_HMIt && m_currentCell == aIt.m_currentCell);
    };

    unsigned int getCellID() {
        return m_currentCell;
    }
};

#endif /*__CellListedHashMapIterator_H */
