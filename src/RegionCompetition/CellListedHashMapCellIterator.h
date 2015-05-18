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
#ifndef CELLLISTEDHASHMAPCELLITERATOR_H
#define	CELLLISTEDHASHMAPCELLITERATOR_H

template <class TCellListedHashMap>
class CellListedHashMapCellIterator {
    /**
     * An iterator enabling (asymetric) particle interactions using a cell list
     * datastructure.
     * Elements stored in the 9 and 27 cells are iterated for 2D and 3D,
     * respectivly (for boundary cells, there are less cells in the neighb.).
     * The iterator position consists of 2 position variables:
     * - m_HMIt: the iterator within the current cell.
     * - m_currentIndexInNeigh: the index(or ID) of the cell within the neighborhood
     * In addition the vector m_Neigh stores all indices (or IDs) of the the
     * cells in the neighborhood.
     */
public:
    typedef TCellListedHashMap CellListedHashMapType;
    typedef CellListedHashMapCellIterator<CellListedHashMapType> Self;

    typedef typename CellListedHashMapType::HashMapType HashMapType;
    typedef typename CellListedHashMapType::HashMapContainerType ContainerType;
    typedef typename CellListedHashMapType::IndexType CLHMKeyType;
    typedef typename CellListedHashMapType::IndexType CellIndexType;
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
    /// the vector m_NeighborCells stores the IDs(linear index) of the neighbors.
    std::vector<unsigned int> m_NeighborCells;
    /// m_currentIndexInNeigh is the index of the current
    /// index in m_Neighborcells corresponding to the iterator position.
    unsigned int m_currentIndexInNeigh;
    /// m_HMIt is the original iterator of the hash-map in the current cell.
    HMIteratorType m_HMIt;
    /// a reference to the nb of cells array in CellListedHashMap used for bound
    /// checking.
    unsigned int* m_NbCells;


public:

      CellListedHashMapCellIterator(CellListedHashMapType* aCLHM,
              CellIndexType aCellIndex, bool aBegin) {
        m_container = &(aCLHM->m_Container);
        m_NbCells = aCLHM->m_NbCells;

        // find the neighboring cells and store their ID (linear index):
        if(TCellListedHashMap::m_DIM== 2) {
//            std::cout << "cell field size: " << m_NbCells[0] << " x " << m_NbCells[1] << std::endl;
//            std::cout << "neighboring cells of [" << aCellIndex[0] << ", " << aCellIndex[1] << "]\n";

            /// iterat the neighborhood
            for (int vY = -1; vY <= 1; vY++) {
                for (int vX = -1; vX <= 1; vX++) {
                    CellIndexType vNeighIndex = aCellIndex;
                    vNeighIndex[0] = aCellIndex[0] + vX;
                    vNeighIndex[1] = aCellIndex[1] + vY;
                    // check the bounds:
                    if(vNeighIndex[0] >= 0 && vNeighIndex[1] >= 0 &&
                            vNeighIndex[0] < m_NbCells[0] && vNeighIndex[1] < m_NbCells[1]) {

                        // add the (linear) index of the neighbor cell:
                        m_NeighborCells.push_back(
                                aCLHM->getLinearizedCellIndexFromCellIndex(vNeighIndex));
//                        std::cout << aCLHM->getLinearizedCellIndexFromCellIndex(vNeighIndex)<< ", ";
                    }
                }
            }
//            std::cout << std::endl;
        } else if (TCellListedHashMap::m_DIM == 3) {
            /// iterat the 1-neighborhood
            for (int vZ = -1; vZ <= 1; vZ++) {
                for (int vY = -1; vY <= 1; vY++) {
                    for (int vX = -1; vX <= 1; vX++) {
                        CellIndexType vNeighIndex = aCellIndex;
                        vNeighIndex[0] = aCellIndex[0] + vX;
                        vNeighIndex[1] = aCellIndex[1] + vY;
                        vNeighIndex[2] = aCellIndex[2] + vZ;
                        /// check the  bounds
                        if (vNeighIndex[0] >= 0 &&
                                vNeighIndex[1] >= 0 &&
                                vNeighIndex[2] >= 0 &&
                                vNeighIndex[0] < m_NbCells[0] &&
                                vNeighIndex[1] < m_NbCells[1] &&
                                vNeighIndex[2] < m_NbCells[2]) {
                            /// add the ID of the neighbor cell.
                            m_NeighborCells.push_back(
                                    aCLHM->getLinearizedCellIndexFromCellIndex(vNeighIndex));
                        }
                    }
                }
            }
        } else {
            assert("CellListedHashMapCellIterator only implemented for 2D and 3D.");
        }


        if (aBegin) {
            m_currentIndexInNeigh = 0;
            unsigned int vCurrentCell = m_NeighborCells[m_currentIndexInNeigh];
            m_HMIt = ((*m_container)[vCurrentCell]).begin();

            // feed forward empty cells:
            while (m_currentIndexInNeigh < m_NeighborCells.size() - 1 &&
                    (((*m_container)[vCurrentCell]).begin() ==
                    ((*m_container)[vCurrentCell]).end())) {

                vCurrentCell = m_NeighborCells[++m_currentIndexInNeigh];
                m_HMIt = ((*m_container)[vCurrentCell]).begin();
            }
        } else {
            /// return the end of the last cell
            m_currentIndexInNeigh = m_NeighborCells.size() - 1;
            unsigned int vCurrentCell = m_NeighborCells.back();
            m_HMIt = ((*m_container)[vCurrentCell]).end();
        }
    };

    //    CellListedHashMapCellIterator(CellListedHashMapType* aCLHM,
    //            KeyType aIndex,
    //            unsigned int aCellIndex) {
//
//        m_container = &(aCLHM->m_Container);
//        m_nbCells = m_container->size();
//
//        m_HMIt = (*m_container)[aCellIndex].find(aIndex);
//        if (m_HMIt != (*m_container)[aCellIndex].end()) {
//            m_currentCell = aCellIndex;
//            //      std::cout << "found!" << std::endl;
//        } else {
//            //      std::cout << "not found!" << std::endl;
//            m_currentCell = (*m_container).size() - 1;
//            m_HMIt = ((*m_container)[m_nbCells - 1]).end();
//        }
//    };

    CellListedHashMapCellIterator & operator++() {

        // Increment the iterator of the original container:
        ++m_HMIt;
        // Get the corresponding cell
        unsigned int vCurrentCell = m_NeighborCells[m_currentIndexInNeigh];

        if (m_HMIt == (*m_container)[vCurrentCell].end()) {
            // jump to the next cell
            vCurrentCell = m_NeighborCells[++m_currentIndexInNeigh];
            
            // feed forward the empty cells
            while (m_currentIndexInNeigh < m_NeighborCells.size() &&
                    (*m_container)[vCurrentCell].begin() ==
                    (*m_container)[vCurrentCell].end()) {
                vCurrentCell = m_NeighborCells[++m_currentIndexInNeigh];
            }
            /// In case of a non-empty cell we stop; else we set the iterator
            /// to the last position of the last cell.
            if (m_currentIndexInNeigh < m_NeighborCells.size()) {
                // set the iterator to the first element in the new cell
                m_HMIt = (*m_container)[vCurrentCell].begin();
            } else {
                // cleanup: set the counter to the last cell
                m_currentIndexInNeigh--; // go back to the last element (equals to size - 1)
                vCurrentCell = m_NeighborCells.back();
                // give back the end of the last cell
                m_HMIt = (*m_container)[vCurrentCell].end();
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

    bool operator!=(const CellListedHashMapCellIterator& aIt) {
        return (m_HMIt != aIt.m_HMIt ||
                m_currentIndexInNeigh != aIt.m_currentIndexInNeigh);
    };

    bool operator==(const CellListedHashMapCellIterator& aIt) {
        return (m_HMIt == aIt.m_HMIt &&
                m_currentIndexInNeigh == aIt.m_currentIndexInNeigh);
    };

    unsigned int getCellID() {
        return m_NeighborCells[m_currentIndexInNeigh];
    }
};


#endif	/* CELLLISTEDHASHMAPCELLITERATOR_H */