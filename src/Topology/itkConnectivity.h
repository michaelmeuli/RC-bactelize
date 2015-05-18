/*=========================================================================
 *
 * This code has been downloaded from 
 * http://www.insight-journal.org/browse/publication/120.
 * 
 * The code was written and contributed by Lamy Julien, 2006.
 * 
 *=========================================================================*/
#ifndef itkConnectivity_h
#define itkConnectivity_h

#include <itkMacro.h>

#include <vnl/vnl_vector_fixed.h>

namespace itk {

    /**
     * @brief Connectivity information.
     *
     * This class describes the k-neighbors of a point in a digital space of
     * dimension n.
     * The definition used for this connectivity is the one using cell
     * decomposition, as defined by Malandain in "On Topology in Multidimensional
     * Discrete Spaces", ftp://ftp.inria.fr/INRIA/publication/RR/RR-2098.ps.gz
     *
     * The advantage of this definition instead of the 4- or 8-connectivity (in 2D)
     * or 6-, 18- and 26-connectivity (in 3D) is that it is consistent and
     * extensible in n-dimensions.
     *
     * In 2D, 4- and 8-connectivity are respectively corresponding to 1- and
     * 0-connectivities. In 3D, 6-, 18- and 26-connectivity are respectively
     * corresponding to 2-, 1- and 0-connectivity.
     */
    template<unsigned int VDim, unsigned int VCellDim>
    class ITK_EXPORT Connectivity {
    public:
        /// @brief Type for a point in n-D.
        typedef vnl_vector_fixed<int, VDim> Point;

        /// @brief Offset in an array.
        typedef unsigned int Offset;

        /// @brief The dimension of the space.
        itkStaticConstMacro(Dimension, unsigned int, VDim);

        /// @brief The dimension of the cell.
        itkStaticConstMacro(CellDimension, unsigned int, VCellDim);


        /// @brief returns the full number of neighbors independent
        ///        on the connectivity type (only dependent on the dimension)
        unsigned int GetNeighborhoodSize() const {
            return m_NeighborhoodSize;
        }

        /// @brief returns the number of neighbors with respect to the
        ///        connectivity type.
        unsigned int GetNumberOfNeighbors() const {
            return m_NumberOfNeighbors;
        }

        /// @brief Accessor to the singleton.
        static Connectivity<VDim, VCellDim> const & GetInstance();

        Point const * const GetNeighborsPoints() const {
            return m_NeighborsPoints;
        }

        Point const * const GetNeighborsOffsets() const {
            m_NeighborsPoints = new Point[m_NumberOfNeighbors];
            return m_NeighborsOffsets;
        }

        typedef typename itk::Point<double, VDim> ItkPointType;

        ItkPointType const * const GetNeighborsITKPoints() const {
            ItkPointType* neighborPoints = new ItkPointType[m_NumberOfNeighbors];
            for(int i = 0; i < m_NumberOfNeighbors; i++) {
                for(unsigned int j = 0; j < VDim; j++) {
                    neighborPoints[i][j] = m_NeighborsPoints[i][j];
                }
            }
            return neighborPoints;
        }

        typedef typename itk::Offset<VDim> ItkOffsetType;

        ItkOffsetType const * const GetNeighborsITKOffsets() const {
            ItkOffsetType* neighborOffsets = new ItkOffsetType[m_NumberOfNeighbors];
            for(unsigned int i = 0; i < m_NumberOfNeighbors; i++) {
                for(unsigned int j = 0; j < VDim; j++) {
                    neighborOffsets[i][j] = m_NeighborsPoints[i][j];
                }
            }
            return neighborOffsets;
        }

        /// @brief Test if two points are neighbors
        bool AreNeighbors(Point const & p1, Point const & p2) const;

        /// @brief Test if two points are neighbors
        bool AreNeighbors(Offset const & o1, Point const & o2) const;

        /// @brief Test if a point is a neighbor of 0
        bool IsInNeighborhood(Point const & p) const;

        /// @brief Test if a point is a neighbor of 0
        bool IsInNeighborhood(Offset const & o) const;

        /// @brief Convert an offset to a point, in a 3x3x3 cube
        Point OffsetToPoint(Offset const offset) const;

        /// @brief Convert a point to an offset, in a 3x3x3 cube
        Offset PointToOffset(Point const p) const;

        /// @brief Compute the number of neighbors within this connectivity.
        static int ComputeNumberOfNeighbors();

    private:
        static Connectivity<VDim, VCellDim> const * m_Instance;

        /// @brief Size of the whole neighborhood (independent of the cell type).
        unsigned int const m_NeighborhoodSize;

        /// @brief Number of neighbors (within the cell).
        unsigned int const m_NumberOfNeighbors;

        /// @brief Neighbors as points.
        Point * const m_NeighborsPoints;

        /// @brief Neighbors as offsets.
        Offset * const m_NeighborsOffsets;

        Connectivity();
        ~Connectivity();

        // Purposedly not implemeted
        Connectivity(Connectivity<VDim, VCellDim> const & other);

        // Purposedly not implemeted
        Connectivity & operator=(Connectivity<VDim, VCellDim> const & other);


        static int factorialrec(int n);
    };

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkConnectivity.hxx"
#endif

#endif // itkConnectivity_h
