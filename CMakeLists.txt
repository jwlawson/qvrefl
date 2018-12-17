cmake_minimum_required(VERSION 3.5)
project(qvrefl)

find_package(Boost REQUIRED)
find_package(BLAS REQUIRED)
find_package(Armadillo REQUIRED)
find_package(QV REQUIRED)
find_package(benchmark)

if(Boost_FOUND AND NOT TARGET Boost::boost)
  add_library(Boost::boost INTERFACE IMPORTED)
  set_target_properties(Boost::boost PROPERTIES
    INTERFACE_INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
  )
endif()

if(Armadillo_FOUND AND NOT TARGET Arma::arma)
  add_library(Arma::arma UNKNOWN IMPORTED)
  set_target_properties(Arma::arma PROPERTIES
    IMPORTED_LOCATION ${ARMADILLO_LIBRARIES}
    INTERFACE_INCLUDE_DIRS ${ARMADILLO_INCLUDE_DIRS}
  )
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
include(QVReflHelpers)

set(_source_dir src)
set(_sources
  ${_source_dir}/cartan.cc
  ${_source_dir}/cartan_exchange_graph.cc
  ${_source_dir}/cartan_iterator.cc
  ${_source_dir}/fully_compatible_check.cc
  ${_source_dir}/mutation_star.cc
  ${_source_dir}/semi_positive_filter.cc
  ${_source_dir}/unique_cartan_filter.cc
  ${_source_dir}/unique_matrix_filter.cc
  ${_source_dir}/util.cc
)

qv_library(TARGET refl
  SOURCES ${_sources}
  PUBLIC_LIBRARIES Boost::boost
                   QV::qv
                   Arma::arma
                   ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES}
)

qv_executable(TARGET qvrefl
  SOURCES ${_source_dir}/qvrefl.cc
  PUBLIC_LIBRARIES refl
)

qv_executable(TARGET cexch
  SOURCES ${_source_dir}/cexch.cc
  PUBLIC_LIBRARIES refl
)

qv_executable(TARGET cmut
  SOURCES ${_source_dir}/cmut.cc
  PUBLIC_LIBRARIES refl
)

if(benchmark_FOUND)
  qv_executable(TARGET benchmark
    SOURCES ${_source_dir}/benchmark.cc
    PUBLIC_LIBRARIES refl benchmark::benchmark
  )
endif()

