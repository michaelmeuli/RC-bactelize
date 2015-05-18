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
#ifndef __itkFrontsCompetitionImageFilter_H_
#define __itkFrontsCompetitionImageFilter_H_

#include <utility>
#include <map>
#include <set>
#include <stack>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include <algorithm>
#include <list>
#include <iomanip>

#include "itkRGBPixel.h"
#include "itkImage.h"
#include "itkImageFunction.h"
#include "itkNumericTraits.h"
#include "itkFixedArray.h"
#include "itksys/hash_map.hxx"
#include "itksys/hash_set.hxx"
#include "CellListedHashMap.h"
#include "CellListedHashMapIterator.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkFloodFilledImageFunctionConditionalIterator.h"
#include "itkFloodFilledImageFunctionConditionalConstIterator.h"
#include "itkConstantBoundaryCondition.h"

#include "itkMultipleThresholdImageFunction.h"

#include "itkBinaryThresholdImageFunction.h"
#include "itkCastImageFilter.h" // write images

#include "itkImageToImageFilter.h" // base class
#include "itkLabelContourImageFilter.h"
#include "itkShiftScaleImageFilter.h" // write images
#include "itkMaskImageFilter.h" 
#include "itkScalarConnectedComponentImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkImageFileWriter.h"

#include "itkConnectivity.h"
#include "itkBackgroundConnectivity.h"
#include "itkTopologicalNumberImageFunction.h"
#include "itkSimplicityByTopologicalNumbersImageFunction.h"

#include "itkContourIndex.h"
#include "itkContourParticle.h"
#include "itkContourParticleWithIndex.h"
#include "itkContourContainer.h"
#include "itkMinimalParticle.h"


#include "itkTimeProbe.h"
#include "itkExceptionObject.h"
#include "itkVector.hxx"

#include "itkRCExternalEnergyBaseClass.h"
#include "itkRCInternalEnergyBaseClass.h"
#include "itkRCParticleInteractionEnergyBaseClass.h"
#include "itkRCGPUEnergyBaseClass.h"


//TODO: insert const wherever it makes sense (where it can be used and makes a difference).
//TODO: template the connectivity types.
//TODO: When smoothing the extenal energies (Sobolev gradients), the performance 
//      will be significantly improved with using symmetric cell lists.


namespace itk {

    /**
     * @brief This filter locally optimizes a energy functional \f$\mathcal{E} \f$
     * by propagating the the fronts of the foreground object.
     *
     * From a given initial front set, fronts are propageted such as to minimize
     * a energy functional \f$\mathcal{E}(\Gamma)\f$ where \f$\Gamma\f$
     * is the front/edge set.
     * The enclosed open set of voxel by a subset of \f$\Gamma\f$ is called a region.
     * 
     * Reference for details and citation:
     * 
     * Cardinale J, Paul G, Sbalzarini I (2012) Discrete region competition for
     * unknown numbers of connected regions. 
     * Image Processing, IEEE Transactions on 21(8):3531 –3545
     *
     * Fusion and fission of regions can be controlled using the concept of
     * digital topology. The digital topology implementation was used from
     * 
     * J. Lamy, “Integrating digital topology in image-processing libraries,”
     * Comput. Methods Programs Biomed., vol. 85, no. 1, pp. 51–58, Jan. 2007.
     * 
     * \ingroup Segmentation
     */
    template <class TInputImage, class TInitImage, class TOutputImage >
    class ITK_EXPORT FrontsCompetitionImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage> {
    public:

        /** Standard class typedefs         */
        typedef FrontsCompetitionImageFilter< TInputImage, TInitImage, TOutputImage> Self;
        typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        
        
        typedef int LabelPixelType; //also need negative numbers for the contour //TODO:template?
        typedef unsigned int LabelAbsPixelType; //TODO:template or cancel
        typedef unsigned char MemoryPixelType;
        typedef itk::Vector<float, TInputImage::ImageDimension> VectorType;

        
        /** Image Types */
        typedef TInputImage InputImageType;
        typedef TInitImage InitImageType;
        typedef TOutputImage OutputImageType;

        typedef Image<LabelPixelType, InputImageType::ImageDimension> LabelImageType;
        typedef Image<LabelAbsPixelType, InputImageType::ImageDimension> LabelAbsImageType;
        typedef Image<MemoryPixelType, InputImageType::ImageDimension> MemoryImageType;


        /** Typedef to describe the output and input image region types. */
        typedef typename InputImageType::RegionType InputImageRegionType;
        typedef typename OutputImageType::RegionType OutputImageRegionType;
        typedef typename InitImageType::RegionType InitImageRegionType;
        typedef typename LabelImageType::RegionType LabelImageRegionType;

        /** Typedef to describe the type of pixel. */
        typedef typename InputImageType::PixelType InputPixelType;
        typedef typename InitImageType::PixelType InitPixelType;
        typedef typename OutputImageType::PixelType OutputPixelType;

        /** Typedef to describe the output and input image index and size types. */
        typedef typename InputImageType::IndexType InputImageIndexType;
        typedef typename InputImageType::SizeType InputImageSizeType;
        typedef typename InputImageType::OffsetType InputImageOffsetType;
        typedef typename OutputImageType::IndexType OutputImageIndexType;
        typedef typename OutputImageType::SizeType OutputImageSizeType;
        typedef typename OutputImageType::OffsetType OutputImageOffsetType;
        typedef typename InitImageType::IndexType InitImageIndexType;
        typedef typename InitImageType::SizeType InitImageSizeType;
        typedef typename InitImageType::OffsetType InitImageOffsetType;
        typedef typename LabelImageType::IndexType LabelImageIndexType;
        typedef typename LabelImageType::SizeType LabelImageSizeType;
        typedef typename LabelImageType::OffsetType LabelImageOffsetType;
        typedef typename LabelImageType::SizeType SizeType;
        
        typedef typename LabelImageType::Pointer LabelImagePointerType;
        typedef typename InitImageType::Pointer InitImagePointerType;
        typedef typename MemoryImageType::Pointer MemoryImagePointerType;

        typedef float InternalPixelType; // float or double; for energy values etc.
        typedef Image<InternalPixelType, InputImageType::ImageDimension> InternalImageType;
        typedef typename InternalImageType::Pointer InternalImagePointerType;
        typedef typename InternalImageType::SizeType InternalImageSizeType;
        typedef typename InternalImageType::IndexType InternalImageIndexType;

        typedef FixedArray<float, InputImageType::ImageDimension> ArrayType;
        typedef FixedArray<float, InputImageType::ImageDimension+1> Centroid_Coeff_Type;

        typedef MultipleThresholdImageFunction<LabelImageType> MultipleThsFunctionType;
        typedef typename MultipleThsFunctionType::Pointer MultipleThsFunctionPointerType;
        
   protected:
        FrontsCompetitionImageFilter();
        virtual ~FrontsCompetitionImageFilter();
        void PrintSelf(std::ostream& aOS, Indent aIndent) const;

        void GenerateData();
        
    public:
        /** Setters and Getters for the parameters */
        itkSetMacro(MaxNbIterations, unsigned int);
        itkGetConstMacro(MaxNbIterations, unsigned int);
        itkGetConstMacro(converged, bool);

        itkSetMacro(MinimalCellSize, SizeType);
        itkGetConstMacro(MinimalCellSize, SizeType);

        itkSetMacro(AllowFusion, bool);
        itkGetConstMacro(AllowFusion, bool);
        itkSetMacro(AllowFission, bool);
        itkGetConstMacro(AllowFission, bool);
        itkSetMacro(AllowHandles, bool);
        itkGetConstMacro(AllowHandles, bool);

        itkSetMacro(UseFastEvolution, bool);
        itkGetConstMacro(UseFastEvolution, bool);
        itkSetMacro(UseSobolevGradients, bool);
        itkGetConstMacro(UseSobolevGradients, bool);
        itkSetMacro(SobolevKernelSigma, ArrayType);
        itkGetConstMacro(SobolevKernelSigma, ArrayType);
        
        itkSetMacro(Verbose, bool);
        itkGetConstMacro(Verbose, bool);

        
        /** Set the initialization image */
        void SetInitImageInput(InitImageType *input) {
            // Process object is not const-correct so the const casting is required.
            this->SetNthInput(1, const_cast<InitImageType *> (input));
        }

        /** Set a forbidden region image */
        void SetForbiddenRegionInput(InitImageType *input) {
            // Process object is not const-correct so the const casting is required.
            this->SetNthInput(2, const_cast<InitImageType *> (input));
            m_UseForbiddenRegion = true;
        }

        /** Get the input image */
        InputImageType * GetDataInput() {
            return static_cast<InputImageType*> (const_cast<DataObject *>
                    (this->ProcessObject::GetInput(0)));
        }

        /** Get the init image */
        InitImageType * GetInitInput() {
            return static_cast<InitImageType*> (const_cast<DataObject *>
                    (this->ProcessObject::GetInput(1)));
        }

        /** Get the forbidden region image */
        InitImageType * GetForbiddenRegionInput() {
            return static_cast<InitImageType*> (const_cast<DataObject *>
                    (this->ProcessObject::GetInput(2)));
        }

    private:
        typedef bool MaskPixelType;
        typedef Image<MaskPixelType, InputImageType::ImageDimension> MaskImageType;
        
        static const unsigned int m_Dim = InputImageType::ImageDimension;

        typedef ContourIndex<InputImageType::ImageDimension> ContourIndexType;
        typedef ContourParticle<InputImageType, LabelImageType> ContourParticleType;
        typedef ContourParticleWithIndex<ContourIndexType, ContourParticleType>
        ContourParticleWithIndexType;
        typedef typename ContourParticleType::EnergyDifferenceType EnergyDifferenceType;
        typedef bool EnergyEventType;

        typedef std::pair<EnergyDifferenceType, EnergyEventType> EnergyDifferenceAndMergePairType;

        //outer contour container types:
        typedef ContourIndexType ContourContainerKeyType;
        typedef ContourParticle<InputImageType, LabelAbsImageType> ContourContainerValueType;

        typedef ContourIndex < 2 > CompetingRegionsPairType;
        typedef std::map<CompetingRegionsPairType, ContourIndexType> CompetingRegionsMapType;
        CompetingRegionsMapType m_CompetingRegionsMap;

        typedef itk::ContourIndexHashFunctor<m_Dim> ContourIndexHashFunctorType;

        // The hash map version
        typedef itksys::hash_map<

        ContourContainerKeyType,
        ContourContainerValueType,
        ContourIndexHashFunctorType,
        ContourIndexHashFunctorType> HashMapType;
        typedef typename HashMapType::iterator HashMapIteratorType;
        typedef typename HashMapType::value_type HashMapValueType;

        // The cell listed hash map version:
        typedef CellListedHashMap<
        m_Dim,
        ContourContainerKeyType,
        ContourContainerValueType,
        ContourIndexHashFunctorType,
        ContourIndexHashFunctorType> CellListedHashMapType;

        // The standard hash map version:
        //        typedef itksys::hash_map<
        //        OuterContourContainerKeyType,
        //        OuterContourContainerValueType,
        //        ContourIndexHashFunctorType,
        //        ContourIndexHashFunctorType> CellListedHashMapType;


        
        // The map iterator:
        typedef typename CellListedHashMapType::iterator CellListedHashMapIteratorType;
        typedef typename CellListedHashMapType::cell_iterator CellListedHashMapCellIteratorType;
    public:
        typedef CellListedHashMapCellIteratorType InteractionParticlesIteratorType;
    private:
        typedef Connectivity<m_Dim, m_Dim - 1 > ForegroundConnectivityType;
        typedef typename BackgroundConnectivity<ForegroundConnectivityType>::Type BackgroundConnectivityType;

     public:
        /* Register energies */    
        typedef RCExternalEnergyBaseClass<LabelImageType, InputImageType> ExternalEnergyType;
        typedef typename ExternalEnergyType::Pointer ExternalEnergyPointerType;
        void RegisterExternalEnergy(ExternalEnergyType* aEnergyPtr){
            ExternalEnergyPointerType vExtEnergy = aEnergyPtr;
            m_ExtEnergies.push_back(vExtEnergy);
        };
        typedef RCInternalEnergyBaseClass<LabelImageType> InternalEnergyType;
        typedef typename InternalEnergyType::Pointer InternalEnergyPointerType;
        void RegisterInternalEnergy(InternalEnergyType* aEnergyPtr){
            InternalEnergyPointerType vIntEnergy = aEnergyPtr;
            m_IntEnergies.push_back(vIntEnergy);
        };
        typedef RCParticleInteractionEnergyBaseClass<LabelImageType, Self> ParticleInteractionEnergyType;
        typedef typename ParticleInteractionEnergyType::Pointer ParticleInteractionEnergyPointerType;
        void RegisterParticleInteractionEnergy(ParticleInteractionEnergyType* aEnergyPtr) {
            ParticleInteractionEnergyPointerType vPIEnergy = aEnergyPtr;
            m_PIEnergies.push_back(vPIEnergy);
        }
        typedef RCGPUEnergyBaseClass<LabelImageType, InputImageType> GPUEnergyType;
        typedef typename GPUEnergyType::GPU_IntType GPU_IntType;
        typedef typename GPUEnergyType::Pointer GPUEnergyPointerType;
        void RegisterGPUEnergy(GPUEnergyType* aEnergyPtr){
            m_GPUEnergy = aEnergyPtr;
        }

    private:
        
        /*
         * Parameters
         */
        
        // Maximum number of iterations (if the optimization doesn't converge)
        unsigned int m_MaxNbIterations;
        
        // Internal parameter, how long should the memory be to detect spatial convergence
        unsigned char m_MemoryLength;
        
        // Determines the Sigma of the Gaussian (if Sobolev gradients are used)
        ArrayType m_SobolevKernelSigma;
        
        // After spatial convergence, how much to reduce the number of accepted points
        float m_AcceptedPointsReductionFactor; 
        
        // Stores the cell size for the cell lists. Only affects performance.
        SizeType m_MinimalCellSize;
        
        // Topology mode
        bool m_AllowFusion;
        bool m_AllowFission;
        bool m_AllowHandles;
        
        /* 
         * modes  
         */
        bool m_UseForbiddenRegion;
        bool m_UseFastEvolution;
        bool m_UseSobolevGradients;
        bool m_Verbose;
        
        /* Flags */
        bool m_converged;
        bool m_spatiallyConverged;
        
        /*
         * Private typedefs
         */
        typedef itk::ConstantBoundaryCondition<LabelImageType> ConstBoundaryConditionType;
        typedef itk::ConstNeighborhoodIterator<LabelImageType, ConstBoundaryConditionType>
        LabelImageConstNeighborhoodIteratorType;
        typedef itk::NeighborhoodIterator<LabelImageType, ConstBoundaryConditionType>
        LabelImageNeighborhoodIteratorType;

        typedef typename ContourIndexType::OffsetType ContourOffsetType;
        typedef typename LabelImageNeighborhoodIteratorType::RadiusType
        LabelImageNeighborhoodIteratorRadiusType;



        /*
         * Private data structures
         */
        CellListedHashMapType m_InnerContourContainer;
        CellListedHashMapType m_Candidates;
        float m_AcceptedPointsFactor; // stores current acceptance reduction factor (while annealing)
        
        typedef std::pair<ContourIndexType, LabelAbsPixelType> SeedType;
        typedef std::set<SeedType> SeedSetType;
        typedef typename SeedSetType::iterator SeedSetIteratorType;
        SeedSetType m_Seeds;

        typedef TopologicalNumberImageFunction<LabelImageType, ForegroundConnectivityType>
        TopologicalNumberCalculatorType;
        typedef typename TopologicalNumberCalculatorType::Pointer TopologicalNumberCalculatorPtrType;
        TopologicalNumberCalculatorPtrType m_TopologicalNumberFunction;


        unsigned int m_NeighborhoodSize_FG_Connectivity;
        unsigned int m_NeighborhoodSize_BG_Connectivity;

        const InputImageOffsetType * m_NeighborsOffsets_BG_Connectivity;
        const InputImageOffsetType * m_NeighborsOffsets_FG_Connectivity;

        LabelAbsPixelType m_MaxNLabels; // past the end counter
        LabelImagePointerType m_LabelImage;
        MemoryImagePointerType m_MemoryImage;

        LabelPixelType m_ForbiddenRegionLabel;

        unsigned int m_iteration_counter; 

        void UpdateStatisticsWhenJump(
                InternalImageIndexType aIndex,
                InputPixelType aData,
                LabelAbsPixelType aFromLabel,
                LabelAbsPixelType aToLabel);

        bool EvaluateCandidatesAndApplyMoves();
        void RebuildAndOptimizeCandidateList(CellListedHashMapType* aCandidateContainer);
        void AddMothersToCandidateList(CellListedHashMapType* aCandidateContainer);
        void AddDaughtersToCandidateList(CellListedHashMapType* aCandidateContainer);
        void GetPropagatingCandidates(
                CellListedHashMapType* aCandidateContainer,
                std::vector<std::vector<int> >* aCoords,
                std::vector<LabelAbsPixelType>* aCandLabels,
                std::vector<unsigned int>* aReferenceCounts);
        void ComputeEnergiesForMothers(CellListedHashMapType* aCandidateContainer);
        void ComputeEnergiesForDaughters(
                CellListedHashMapType* aCandidateContainer,
                std::vector<std::vector<GPU_IntType> >* aCoords,
                std::vector<LabelAbsPixelType>* aCandLabels,
                std::vector<EnergyDifferenceType>* aInternalEnergy,
                std::vector<EnergyDifferenceType>* aExternalEnergy,
                std::vector<EnergyEventType>* aMerge);

        void MinimizeEnergy(
                CellListedHashMapType* aCandidateContainer,
                std::vector<std::vector<GPU_IntType> >* aCoords,
                std::vector<LabelAbsPixelType>* aCandLabels,
                std::vector<EnergyDifferenceType>* aInternalEnergy,
                std::vector<EnergyDifferenceType>* aExternalEnergy,
                std::vector<EnergyEventType>* aMerge,
                std::vector<unsigned int>* aReferenceCounts);

        void AddUpEnergies(CellListedHashMapType* aCandidateContainer);


    
        unsigned int InitializeGPUArrays(
                CellListedHashMapType* aCandidateContainer,
                std::vector<std::vector<GPU_IntType> >* aCoords,
                std::vector<LabelAbsPixelType>* aCandLabels,
                unsigned int ***aGPUCoordinates,
                unsigned int** aGPUToLabel,
                float** aGPUOutInternalEnergy,
                float** aGPUOutExternalEnergy,
                unsigned int** aGPUOutMerge);

        void DistroyGPUArrays(
                unsigned int **aGPUCoordinates,
                unsigned int* aGPUToLabel,
                float* aGPUOutInternalEnergy,
                float* aGPUOutExternalEnergy,
                unsigned int* aGPUMerge);

        void CopyGPUResults(
                CellListedHashMapType* aCandidateContainer,
                std::vector<EnergyDifferenceType>* aInternalEnergy,
                std::vector<EnergyDifferenceType>* aExternalEnergy,
                std::vector<EnergyEventType>* aMerge,
                float* aGPUOutInternalEnergy,
                float* aGPUOutExternalEnergy,
                unsigned int* aGPUOutMerge);


        void SmoothExternalEnergies(CellListedHashMapType* aCandidateContainer);
        void CleanUp();
        void CopyRegionsAndAllocateOutput();
        void ChangeContourParticleLabelToCandidateLabel(
                ContourIndexType aIndex, ContourContainerValueType* aIt);
        void AddNeighborsAtRemove(
                LabelAbsPixelType aAbsLabel,
                ContourContainerKeyType aIndex);
        void MaintainNeighborsAtAdd(
                LabelAbsPixelType aLabelAbs,
                ContourContainerKeyType aIndex);
        void MaintainNeighborsAtAdd2(CellListedHashMapIteratorType aIt);
        void RemoveContourParticle(
                CellListedHashMapType* aTempAddContainer,
                LabelAbsPixelType aAbsLabel,
                CellListedHashMapIteratorType aIt);

        void AddContourParticle(
                CellListedHashMapType* aTempAddContainer,
                CellListedHashMapIteratorType aIt);

        void UpdateMemoryImage();

        void FilterCandidatesContainerUsingRanks(CellListedHashMapType* aContainer);
        float SumAllEnergies(CellListedHashMapType* aContainer);
//        void RemoveSinglePointRegions();

        void WriteLabelImageToOutput();
        

        
        std::pair<EnergyDifferenceType, EnergyEventType> CalculateExternalEnergyDifferenceForLabel(
                ContourContainerKeyType aContourIndex,
                ContourContainerValueType * aPtr,
                LabelAbsPixelType aL);
        EnergyDifferenceType CalculateInternalEnergyDifferenceForLabel(
                ContourContainerKeyType aContourIndex,
                ContourContainerValueType * aPtr,
                LabelAbsPixelType aL);
        EnergyDifferenceType CalculateAverageInternalMotherEnergyForLabel(
                CellListedHashMapType* aCandidateContainer,
                ContourContainerValueType * aPtr,
                LabelAbsPixelType aL);
        bool CalculateMergingEnergyForLabel(
                LabelAbsPixelType aL1,
                LabelAbsPixelType aL2);

        void FreeLabelStatistics(LabelAbsPixelType);


        void InitializeLabelImageAndContourContainer(InitImagePointerType aInitImage);

        bool DoOneIteration();

        
        inline bool IsEnclosedByLabel_FGConnectivity(const LabelImageNeighborhoodIteratorType& aIt) const;
        inline bool IsEnclosedByLabel_FGConnectivity(const LabelImageIndexType& aIndex, LabelAbsPixelType aLabel) const;
        
        inline bool HasBGNeighbor(const LabelImageIndexType& aIndex) const;
        inline bool IsBoundaryPoint(const LabelImageIndexType& aIndex) const;
        inline bool IsSingleFGPoint(const LabelImageIndexType& aIndex, LabelAbsPixelType aAbsLabel) const;

        void RelabelAllAdjacentRegionsAfterTopologicalChange(LabelImagePointerType aLabelImage,
                InputImageIndexType aIndex, LabelAbsPixelType aMaxLabel);
        void RegisterSeedsAfterSplit(
                LabelImagePointerType aLabelImage, ContourIndexType aIndex,
                LabelAbsPixelType aLabel, CellListedHashMapType* aCandidateContainer);
        void RelabelRegionsAfterSplit(
                LabelImagePointerType aLabelImage, InputImageIndexType aIndex,
                LabelAbsPixelType aLabel);
        void RelabelRegionsAfterFusions(LabelImagePointerType aLabelImage,
                InputImageIndexType aIndex, LabelAbsPixelType aL1);
        void ForestFire(LabelImagePointerType aLabelImage,
                InputImageIndexType aIndex,
                MultipleThsFunctionPointerType aMultiThsFunction,
                LabelAbsPixelType);
        void RenewStatistics(LabelImagePointerType aInitImage);


        void PrepareEnergyCaluclation();
        void PrepareEnergyCaluclationForEachIteration();


    private:
        std::vector<ExternalEnergyPointerType> m_ExtEnergies;
        std::vector<InternalEnergyPointerType> m_IntEnergies;
        std::vector<ParticleInteractionEnergyPointerType> m_PIEnergies;
        // Only one gpu energy is allowed. In order to overcome this limitation
        // each energy would have to invalidate the cl state and the whole 
        // schabang with it. The idea here rather is that one energy completely
        // determines the energy differences for all particles in one go. 
        GPUEnergyPointerType m_GPUEnergy; 


    private:
        /** Debugging help variables */
        
        std::vector<double> m_TimingsHist;
        std::vector<unsigned int> m_MergingHist;
        std::vector<unsigned int> m_FramesHist;

        void WriteLabelImageToFile(char* aFilename);
        void WriteLabelImageToFile2(char* aFilename);

        void WriteLabelImageTypeToFile(const char* aFilename, LabelImageType* aImage);
        void WriteLabelAbsImageTypeToFile(const char* aFilename, LabelAbsImageType* aImage);
        void WriteMemoryImageToFile(char* aFilename);

        void WriteInputImageTypeToFile(const char* aFilename, InputImageType* aImage, float aScale);
        void WriteMaskImageTypeToFile(const char* aFilename, MaskImageType* aImage, float aScale);
        void WriteEnergyGradientImageToFile(char* aFilename);
        void WriteEnergyDifferencesToImageFile(CellListedHashMapType * aContainer, const char* aFilename);
        void WriteDebugImages(unsigned int aIterationNumber);
        void WriteContourParticleContainer(const char* aFilename,
                CellListedHashMapType &aContainer);
        void DebugCellListedHashMapIterator(const char* aFilename, CellListedHashMapType &aContainer);

        std::vector<double> m_TimeHistRelabeling;
        std::vector<double> m_TimeHistEnergyComputation;
        std::vector<double> m_TimeHistRebuilding;
        std::vector<double> m_TimeHistMinimization;
        std::vector<double> m_TimeHistFrontPropagation;
    };

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFrontsCompetitionImageFilter.hxx"
#endif

#endif /* __itkFrontsCompetitionImageFilter_H_ */
