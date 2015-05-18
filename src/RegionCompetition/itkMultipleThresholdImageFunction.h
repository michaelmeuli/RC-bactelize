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
#ifndef __itkMultipleThresholdImageFunction_h
#define __itkMultipleThresholdImageFunction_h

#include "itkImageFunction.h"
#include <utility>
#include <vector>

namespace itk {

    /** \class MultipleThresholdImageFunction
     * \brief Returns true is the value of an image lies within a range
     *        of thresholds
     * This ImageFunction returns true (or false) if the pixel value lies
     * within (outside) a a set of lower and upper threshold values. The threshold
     * range can be set with the ThresholdBelow, ThresholdBetween or
     * ThresholdAbove methods.  The input image is set via method
     * SetInputImage().
     *
     * Methods Evaluate, EvaluateAtIndex and EvaluateAtContinuousIndex
     * respectively evaluate the function at an geometric point, image index
     * and continuous image index.
     *
     * \ingroup ImageFunctions
     *
     */
    template <class TInputImage, class TCoordRep = float>
    class ITK_EXPORT MultipleThresholdImageFunction :
    public ImageFunction<TInputImage, bool, TCoordRep> {
    public:
        /** Standard class typedefs. */
        typedef MultipleThresholdImageFunction Self;
        typedef ImageFunction<TInputImage, bool, TCoordRep> Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        /** Run-time type information (and related methods). */
        itkTypeMacro(MultipleThresholdImageFunction, ImageFunction);

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** InputImageType typedef support. */
        typedef typename Superclass::InputImageType InputImageType;

        /** Typedef to describe the type of pixel. */
        typedef typename TInputImage::PixelType PixelType;

        /** Dimension underlying input image. */
        itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

        /** Point typedef support. */
        typedef typename Superclass::PointType PointType;

        /** Index typedef support. */
        typedef typename Superclass::IndexType IndexType;

        /** ContinuousIndex typedef support. */
        typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

        typedef typename std::pair<PixelType, PixelType> ThresholdPairType;

        typedef typename std::vector<ThresholdPairType> ThresholdPairVectorType;

        /** MultipleThreshold the image at a point position
         *
         * Returns true if the image intensity at the specified point position
         * satisfies the threshold criteria.  The point is assumed to lie within
         * the image buffer.
         *
         * ImageFunction::IsInsideBuffer() can be used to check bounds before
         * calling the method. */

        virtual bool Evaluate(const PointType& point) const {
            IndexType index;
            this->ConvertPointToNearestIndex(point, index);
            return ( this->EvaluateAtIndex(index));
        }

        /** MultipleThreshold the image at a continuous index position
         *
         * Returns true if the image intensity at the specified point position
         * satisfies the threshold criteria.  The point is assumed to lie within
         * the image buffer.
         *
         * ImageFunction::IsInsideBuffer() can be used to check bounds before
         * calling the method. */
        virtual bool EvaluateAtContinuousIndex(
                const ContinuousIndexType & index) const {
            IndexType nindex;

            this->ConvertContinuousIndexToNearestIndex(index, nindex);
            return this->EvaluateAtIndex(nindex);
        }

        /** MultipleThreshold the image at an index position.
         *
         * Returns true if the image intensity at the specified point position
         * satisfies the threshold criteria.  The point is assumed to lie within
         * the image buffer.
         *
         * ImageFunction::IsInsideBuffer() can be used to check bounds before
         * calling the method. */
        virtual bool EvaluateAtIndex(const IndexType & index) const {
            PixelType value = this->GetInputImage()->GetPixel(index);
            for (unsigned int vI = 0; vI < m_NThresholds; vI++) {
                if (m_Thresholds[vI].first <= value && value <= m_Thresholds[vI].second) {
                    return true;
                }
            }
            return false;
        }

        /** Get the lower threshold value. */
        //itkGetConstReferenceMacro(Thresholds,std::vector<std::pair<PixelType,PixelType>>);

        /** Values that lie between lower and upper inclusive are inside. */
        void AddThresholdBetween(PixelType lower, PixelType upper);

        void ClearThresholds() {
            m_Thresholds.clear();
            m_NThresholds = 0;
        }

    protected:
        MultipleThresholdImageFunction();

        ~MultipleThresholdImageFunction() {
        };
        void PrintSelf(std::ostream& os, Indent indent) const;

    private:
        MultipleThresholdImageFunction(const Self&); //purposely not implemented
        void operator=(const Self&); //purposely not implemented

        ThresholdPairVectorType m_Thresholds;
        unsigned int m_NThresholds; // to not call size() of the vector at each evaluation.
    };

} // end namespace itk


// Define instantiation macro for this template.
#define ITK_TEMPLATE_MultipleThresholdImageFunction(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT MultipleThresholdImageFunction< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef MultipleThresholdImageFunction< ITK_TEMPLATE_2 x > \
                               MultipleThresholdImageFunction##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkMultipleThresholdImageFunction+-.h"
#endif

#if ITK_TEMPLATE_TXX
#include "itkMultipleThresholdImageFunction.hxx"
#endif

#endif
