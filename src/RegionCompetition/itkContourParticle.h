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
 * File:   itkContourPoint.h
 * Author: janickc
 *
 * Created on April 20, 2010, 1:32 PM
 */

#ifndef _ITKCONTOURPOINT_H
#define	_ITKCONTOURPOINT_H

#include <list>
#include <utility>
#include "itkImage.h"
#include "itkContourIndex.h"
//#include "itkContourPointCandidateElement.h"

namespace itk {

    template <class TDataImage, class TLabelImage>
    class ContourParticle {
    public:
        typedef TDataImage DataImageType;
        typedef TLabelImage LabelImageType;
        typedef typename DataImageType::PixelType DataPixelType;
        typedef typename LabelImageType::PixelType LabelPixelType;
        typedef float EnergyDifferenceType;
        typedef std::set<LabelPixelType>  TestedLabelsListType;
//        typedef ContourPointCandidateElementType ContourPointCandidateElement;
//        typedef std::set<ContourPointCandidateElement> ContourPointCandidateListType;


        // TODO: change list types to std::list<unsigned int> and the candidate list to a deque.
        //        typedef std::list<unsigned int > MotherIndexListType;
        //        typedef std::list<unsigned int > DaughterIndexListType;
        typedef std::list<ContourIndex<TDataImage::ImageDimension> > MotherIndexListType;
        typedef std::list<ContourIndex<TDataImage::ImageDimension> > DaughterIndexListType;

        inline ContourParticle() :
        m_data(0),
        m_label(0),
        m_candidateLabel(0),
        m_energyDifference(NumericTraits<EnergyDifferenceType>::max()),
        m_externalEnergyDifference(NumericTraits<EnergyDifferenceType>::max()),
        m_modifiedCounter(1),
        m_isDaughter(false),
        m_isMother(false){
        };

        inline ContourParticle(
                DataPixelType aD,
                LabelPixelType aL,
                LabelPixelType aCandL) :
        m_data(aD),
        m_label(aL),
        m_candidateLabel(aCandL),
        m_energyDifference(NumericTraits<EnergyDifferenceType>::max()),
        m_externalEnergyDifference(NumericTraits<EnergyDifferenceType>::max()),
        m_modifiedCounter(1),
        m_isDaughter(false),
        m_isMother(false) {
        };

//        void SortCandidateList(){
//            m_candidates.sort();
//        }

//        inline bool HasLabelBeenTested(LabelPixelType aLabel) {
//            return (m_testedList.find(aLabel) != m_testedList.end());
//        }
//
//        inline void LabelHasBeenTested(LabelPixelType aLabel) {
//            m_testedList.insert(aLabel);
//        }

        DataPixelType m_data;
        LabelPixelType m_label;

        LabelPixelType m_candidateLabel;
        EnergyDifferenceType m_energyDifference;
        EnergyDifferenceType m_externalEnergyDifference;
        EnergyDifferenceType m_filteredEnergyDifference;

        
        unsigned int m_referenceCount;
        unsigned int m_modifiedCounter;
//        TestedLabelsListType m_testedList;
        //        ContourPointCandidateListType m_candidates;

        MotherIndexListType m_motherIndices;
        DaughterIndexListType m_daughterIndices;

        bool m_processed;
        bool m_isMother;
        bool m_isDaughter;
      //  bool m_merge;
    };



    template <class TDataImage, class TLabelImage>
    std::ostream & operator<<(
    std::ostream &os,
    const ContourParticle <TDataImage, TLabelImage> &aCPoint) {
        os << "label: \t" << aCPoint.m_label << std::endl;
        os << "data: \t" << aCPoint.m_data << std::endl;
        os << "candidate label: \t" << aCPoint.m_candidateLabel << std::endl;
        os << "energyDifference: \t" << aCPoint.m_energyDifference << std::endl;
        os << "is mother: \t" << aCPoint.m_isMother << std::endl;
        os << "is daughter: \t" << aCPoint.m_isDaughter << std::endl;
        os << "reference count: \t" << aCPoint.m_referenceCount << std::endl;
        return os;
    }


    

} //end namespace itk

#endif	/* _ITKCONTOURPOINT_H */

