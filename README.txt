This program is based on Region Competition:
J. Cardinale, G. Paul, and I. F. Sbalzarini. Discrete region competition for unknown numbers of connected regions. IEEE Trans. Image Process., 21(8):3531–3545, 2012.
http://mosaic.mpi-cbg.de/?q=downloads/RCITK

Lif_extractor.ijm is to be used with Fiji and is taken from Supplementary Data to:
A. Rizk, G. Paul, P. Incardona, M. Bugarski, M. Mansouri, A. Niemann, U. Ziegler, P. Berger, and I. F. Sbalzarini. Segmentation and quantification of subcellular structures in fluorescence microscopy images using Squassh. Nature Protocols, 9(3):586–596, 2014.
http://mosaic.mpi-cbg.de/?q=downloads/imageJ

Region Competition is also available as plugin for Fiji:
http://mosaic.mpi-cbg.de/?q=downloads/imageJ
http://mosaic.mpi-cbg.de/MosaicToolboxSuite/MosaicToolsuiteTutorials.html

This program uses the "Insight Segmentation and Registration Toolkit (ITK) 4.8" with:
Module_ITKReview=ON
Module_LesionSizingToolkit=ON
Module_SCIFIO=ON

boost_1_49_0 also needs to be installed:
http://www.boost.org/users/history/version_1_49_0.html

To avoid SCIFIO throwhing exeption because of memory limits:
export JAVA_FLAGS=-Xmx5400m
To have it permanently add the following line to /etc/environment:
JAVA_FLAGS=-Xmx5400m

Use of itkSCIFIOImageIOTest with ome-tiff files:
./SCIFIOTestDriver itkSCIFIOImageIOTest /path/to/dead-A.ome.tiff /path/to/dead-A.ome.tiff -w -a -d 5

