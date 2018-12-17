cmake_minimum_required(VERSION 3.5)

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

add_library(refl ${_sources})
target_link_libraries(refl
  PUBLIC Boost::boost
         QV::qv
         ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES}
         ${ARMADILLO_LIBRARIES}
)
target_include_directories(refl
  PUBLIC ${ARMADILLO_INCLUDE_DIRS}
  PRIVATE ${CMAKE_SOURCE_DIR}/include
)

add_executable(qvrefl ${_source_dir}/qvrefl.cc)
target_link_libraries(qvrefl PUBLIC refl)
target_include_directories(qvrefl PRIVATE ${CMAKE_SOURCE_DIR}/include)

add_executable(cexch ${_source_dir}/cexch.cc)
target_link_libraries(cexch PUBLIC refl)
target_include_directories(cexch PRIVATE ${CMAKE_SOURCE_DIR}/include)

add_executable(cmut ${_source_dir}/cmut.cc)
target_link_libraries(cmut PUBLIC refl)
target_include_directories(cmut PRIVATE ${CMAKE_SOURCE_DIR}/include)

if(benchmark_FOUND)
  add_executable(benchmark ${_source_dir}/benchmark.cc)
  target_link_libraries(benchmark PUBLIC refl benchmark::benchmark)
  target_include_directories(benchmark PRIVATE ${CMAKE_SOURCE_DIR}/include)
endif()

