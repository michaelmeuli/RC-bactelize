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
#ifndef __CellListedHashMap_H
#define __CellListedHashMap_H


#include <vector>
#include "itksys/hash_map.hxx"
#include "CellListedHashMapIterator.h"
#include "CellListedHashMapCellIterator.h"

template <unsigned int VDim, class TIndex> class CellListedHashMapHelper;

template <unsigned int VDim, class TIndex, class TValue, class THashFcn, class TEqualKey>
class CellListedHashMap {
public:
    typedef CellListedHashMap <VDim, TIndex, TValue, THashFcn, TEqualKey> Self;
    typedef TIndex IndexType;
    typedef TValue ValueType;
    typedef THashFcn HashFunctionType;
    typedef TEqualKey EqualKeyType;
    typedef itksys::hash_map<IndexType, ValueType, HashFunctionType, EqualKeyType> HashMapType;

    typedef std::vector<HashMapType> HashMapContainerType;
    typedef CellListedHashMapIterator<Self> iterator;
    typedef CellListedHashMapCellIterator<Self> cell_iterator;
    // we should be able to define teh const_cell_iterator as follows:
//    typedef CellListedHashMapCellIterator<const Self> const_cell_iterator;
    typedef typename HashMapType::size_type size_type;
//    typedef typename HashMapType::iterator cell_iterator;

private:
    typedef CellListedHashMapHelper<VDim, IndexType> IndexCalculatorType;


    template <class T> friend class CellListedHashMapIterator;
    template <class T> friend class CellListedHashMapCellIterator;
//    friend std::ostream & operator<<(std::ostream& aOS, const Self &aHashMap);

    HashMapContainerType m_Container;
    unsigned int m_CellSize[3];
    unsigned int m_RegionSize[3];
    unsigned int m_NbCells[3];
    int m_RegionIndex[3];
    unsigned int m_LinNbCells;
    static const unsigned int m_DIM = TIndex::IndexDimension;

//    IndexCalculatorType m_IndexCalculator;

public:
    CellListedHashMap() {
        unsigned int vDefaultCellSize[3] = {20, 20, 10};
        unsigned int vDefaultRegionSize[3] = {200, 200, 10};
        int vDefaultRegionIndex[3] = {0, 0, 0};
        ClearAndConstruct(vDefaultCellSize, vDefaultRegionSize, vDefaultRegionIndex);
    };

    CellListedHashMap(unsigned int* aCellSize, unsigned int* aRegionSize, int* aRegionIndex) {
        ClearAndConstruct(aCellSize, aRegionSize, aRegionIndex);
    };

private:
    void ClearAndConstruct(unsigned int* aCellSize, unsigned int* aRegionSize, int* aRegionIndex){
        m_Container.clear();
        for (unsigned int vD = 0; vD < m_DIM; vD++) {
            m_CellSize[vD] = aCellSize[vD];
            m_RegionSize[vD] = aRegionSize[vD];
            m_RegionIndex[vD] = aRegionIndex[vD];
            m_NbCells[vD] = (float) m_RegionSize[vD] / (float) m_CellSize[vD];

            if (m_RegionSize[vD] % m_CellSize[vD] != 0) {
                m_NbCells[vD]++;
            }
        }
        m_LinNbCells = 1;
        for (unsigned int vD = 0; vD < m_DIM; vD++) {
            m_LinNbCells *= m_NbCells[vD];
        }

        /// Construct the hashmaps:
        for (unsigned int vI = 0; vI < m_LinNbCells; vI++) {
            HashMapType vHM;
            m_Container.push_back(vHM);
        }
//        std::cout << "container size when init:" << m_Container.size() << std::endl;
        
    };

public:

    void SetCellAndRegionSize(unsigned int* aCellSize, unsigned int* aRegionSize,
            int* aRegionIndex) {

        ClearAndConstruct(aCellSize, aRegionSize, aRegionIndex);

//        std::cout << "DEBUGGING output for CellLists: \n"
//                << "CellSize:   [" << m_CellSize[0]   << "," << m_CellSize[1]   << "," << m_CellSize[2] << "]\n"
//                << "RegionSize: [" << m_RegionSize[0] << "," << m_RegionSize[1] << "," << m_RegionSize[2] << "]\n"
//                << "nb_cells:   [" << m_NbCells[0]    << "," << m_NbCells[1]    << "," << m_NbCells[2] << "]\n";


        //TODO: the next block segfaults if the container is non-empty
        /*
        /// Copy the elements into arrays
        unsigned int vSize = this->size();
        IndexType* vIndices;
        ValueType* vValues;
        if (vSize > 0) {
            vIndices = new IndexType[vSize];
            vValues = new ValueType[vSize];
        }
        iterator vIt = this->begin();
        iterator vItEnd = this->end();

        unsigned int vI = 0;
        for(; vIt != vItEnd; ++vIt) {
            vI++;
            vIndices[vI] = vIt->first;
            vValues[vI] = vIt->second;
        }

        /// Clear the container and reset members:
        this->ClearAndConstruct(aCellSize, aRegionSize, aRegionIndex);

        /// Refill the container:
        for(vI = 0; vI < vSize; vI++) {
            (*this)[vIndices[vI]] = vValues[vI];
        }

        if (vSize > 0) {
            delete vIndices;
            delete vValues;
        }
         * */
    };
    
    iterator begin() {
        return iterator(this, true);
    };

    iterator end() {
        return iterator(this, false);
    };

    iterator find(const IndexType& aKey) {
        return iterator(this, aKey, getLinearizedCellIndexFromParticle(aKey));
    };

    cell_iterator cell_begin(const IndexType& aParticleIndex) {
        return cell_iterator(this, getCellIndexFromParticleIndex(aParticleIndex), true);
    };

    cell_iterator cell_end(const IndexType& aParticleIndex) {
        return cell_iterator(this, getCellIndexFromParticleIndex(aParticleIndex), false);
    };

    ValueType & operator[](const IndexType& aKey) {
//        if(getLinearizedCellIndex(aKey) >= m_Container.size()) {
//            std::cout << "Warning!: out of bounds CellListedHashMap: Index: " <<
//                    aKey <<
//                    "; index: " << getLinearizedCellIndex(aKey) << std::endl;
//            std::cout << "container size: " << m_Container.size();
//        }
        return (m_Container[getLinearizedCellIndexFromParticle(aKey)])[aKey];
    };
    

    size_type erase(const IndexType& aKey) {
        return m_Container[getLinearizedCellIndexFromParticle(aKey)].erase(aKey);
    };

    void erase(const iterator aIt) {
        m_Container[aIt.m_currentCell].erase(aIt.m_HMIt);
    };


    void clear() {
        for (unsigned int vI = 0; vI < m_LinNbCells; vI++) {
            m_Container[vI].clear();
        }
    };

    size_type size() const {
        size_type vCount = 0;
        for (unsigned int vI = 0; vI < m_LinNbCells; vI++) {
            vCount += m_Container[vI].size();
        }
        return vCount;
    }

    bool empty() {
        for (unsigned int vI = 0; vI < m_LinNbCells; vI++) {
            if (!m_Container[vI].empty()) {
                return false;
            }
        }
        return true;
    };

    inline unsigned int getLinearizedCellIndexFromParticle(const IndexType& aParticleIndex) {
        return IndexCalculatorType::getLinearizedCellIndexFromParticle(
                aParticleIndex, m_CellSize, m_NbCells, m_RegionIndex);
    }

    inline IndexType getCellIndexFromParticleIndex(const IndexType& aParticleIndex) {
        return IndexCalculatorType::getCellIndexFromParticleIndex(
                aParticleIndex, m_CellSize, m_RegionIndex);
    }

    inline unsigned int getLinearizedCellIndexFromCellIndex(const IndexType& aCellIndex) {
        return IndexCalculatorType::getLinearizedCellIndexFromCellIndex(aCellIndex, m_NbCells);
    }

};

template <unsigned int VDim, class TIndex, class TValue, class THashFcn, class TEqualKey>
std::ostream & operator<<(std::ostream& aOS,
 CellListedHashMap<VDim, TIndex, TValue, THashFcn, TEqualKey> &aHashMap) {
    typedef CellListedHashMap<VDim, TIndex, TValue, THashFcn, TEqualKey> HashMapType;
    typedef typename HashMapType::iterator HashMapIteratorType;
    HashMapIteratorType vIt = aHashMap.begin();
    HashMapIteratorType vEnd = aHashMap.end();

    for (; vIt != vEnd; ++vIt) {

        aOS << vIt->first << " " << vIt.getCellNb() << std::endl;
    }
    return aOS;
}

template <unsigned int VDim, class TIndex>
class CellListedHashMapHelper{
    typedef TIndex IndexType;
public:
    static inline unsigned int getLinearizedCellIndexFromParticle(
            const IndexType& aKey,
            unsigned int* a_CellSize,
            unsigned int* a_NbCells){
        assert("CellListedHashMap not implemented for this dimension.");
        return 0;
    }

    static inline unsigned int getLinearizedCellIndexFromCellIndex(
            const IndexType& aCellIndex,
            unsigned int* a_NbCells) {
        assert("CellListedHashMap not implemented for this dimension.");
        return 0;
    }

    static inline IndexType getCellIndexFromParticleIndex(
            const IndexType& aParticleIndex,
            unsigned int* a_CellSize,
            int* a_RegionIndex) {
        assert("CellListedHashMap not implemented for this dimension.");
        IndexType vI;
        return vI;
    }
};

template <class TIndex>
class CellListedHashMapHelper < 3, TIndex> {
    typedef TIndex IndexType;
public:

    static inline unsigned int getLinearizedCellIndexFromParticle(
            const IndexType& aKey,
            unsigned int* a_CellSize,
            unsigned int* a_NbCells,
            int* a_RegionIndex) {
        return ((aKey[0] - a_RegionIndex[0]) / a_CellSize[0]) // x-coordinate
                + ((aKey[1] - a_RegionIndex[1]) / a_CellSize[1]) * a_NbCells[0] // y-coordinate
                + ((aKey[2] - a_RegionIndex[2]) / a_CellSize[2]) * a_NbCells[0] * a_NbCells[1];
    }

    static inline unsigned int getLinearizedCellIndexFromCellIndex(
            const IndexType& aCellIndex,
            unsigned int* a_NbCells) {

        return ((aCellIndex[0])) // x-coordinate
                + (aCellIndex[1]) * a_NbCells[0] // y-coordinate
                + (aCellIndex[2]) * a_NbCells[0] * a_NbCells[1];
    }

    static inline IndexType getCellIndexFromParticleIndex(
            const IndexType& aParticleIndex,
            unsigned int* a_CellSize,
            int* a_RegionIndex) {
        IndexType vCellIndex;
        vCellIndex[0] = (aParticleIndex[0] - a_RegionIndex[0]) / a_CellSize[0];
        vCellIndex[1] = (aParticleIndex[1] - a_RegionIndex[1]) / a_CellSize[1];
        vCellIndex[2] = (aParticleIndex[2] - a_RegionIndex[2]) / a_CellSize[2];
        return vCellIndex;
    }

};

template <class TIndex>
class CellListedHashMapHelper < 2, TIndex> {
    typedef TIndex IndexType;
public:

    static inline unsigned int getLinearizedCellIndexFromParticle(
            const IndexType& aKey,
            unsigned int* a_CellSize,
            unsigned int* a_NbCells,
            int* a_RegionIndex) {
        return ((aKey[0] - a_RegionIndex[0]) / a_CellSize[0]) // x-coordinate
                + ((aKey[1] - a_RegionIndex[1]) / a_CellSize[1]) * a_NbCells[0];
    }

    static inline unsigned int getLinearizedCellIndexFromCellIndex(
            const IndexType& aCellIndex,
            unsigned int* a_NbCells) {

        return ((aCellIndex[0])) // x-coordinate
                + (aCellIndex[1]) * a_NbCells[0]; // y-coordinate
    }

    static inline IndexType getCellIndexFromParticleIndex(
            const IndexType& aParticleIndex,
            unsigned int* a_CellSize,
            int* a_RegionIndex) {
        IndexType vCellIndex;
        vCellIndex[0] = (aParticleIndex[0] - a_RegionIndex[0]) / a_CellSize[0];
        vCellIndex[1] = (aParticleIndex[1] - a_RegionIndex[1]) / a_CellSize[1];
        return vCellIndex;
    }
};


//#include "CellListedHashMap.hxx"
#endif /* __CellListedHashMap_H */
