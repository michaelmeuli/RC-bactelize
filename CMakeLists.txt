CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# The main project, this defines the projects root src path
PROJECT(RC-bactelize_project)

# Find boost
FIND_PACKAGE(Boost REQUIRED)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

# Find ITK.
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

# Cuda OpenCL path
OPTION (USE_GPU OFF)
SET(CUDA_DIR CACHE PATH "Cuda binary directory")
if (USE_GPU)
	if (CUDA_DIR)
		add_definitions(-DUSE_GPU)
	else()
		MESSAGE(FATAL_ERROR "USE_GPU is on. Please choose a cuda dir")
	endif()
endif()

# Add the frontspropagation code
SET(FC_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
SET(SRC_DIRS ${FC_SRC_DIR}/RegionCompetition )
LIST(APPEND SRC_DIRS ${FC_SRC_DIR}/Topology )
LIST(APPEND SRC_DIRS ${FC_SRC_DIR}/Energies )
LIST(APPEND SRC_DIRS ${FC_SRC_DIR}/Initializers)
if (USE_GPU)
	LIST(APPEND SRC_DIRS ${FC_SRC_DIR}/OpenCLSupport )
	INCLUDE_DIRECTORIES(${CUDA_DIR}/include)
endif()
INCLUDE_DIRECTORIES(${SRC_DIRS})

SET(CMAKE_VERBOSE_MAKEFILE ON)

ADD_EXECUTABLE(RC-bactelize main.cxx)

TARGET_LINK_LIBRARIES(RC-bactelize boost_program_options ${ITK_LIBRARIES})
if (USE_GPU)
	TARGET_LINK_LIBRARIES(RC-bactelize OpenCL)
endif()
