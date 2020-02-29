add_executable(CalibrateOsm CalibrateOsm.cc)
target_link_libraries(CalibrateOsm ${RoutingKit_LIBRARY})

add_executable(ComputeDijkstraRanks ComputeDijkstraRanks.cc)
target_compile_definitions(ComputeDijkstraRanks PRIVATE CSV_IO_NO_THREAD)

add_executable(ComputeDistances ComputeDistances.cc)
target_compile_definitions(ComputeDistances PRIVATE CSV_IO_NO_THREAD)

add_executable(ComputeUnionBoundary ComputeUnionBoundary.cc)
target_compile_definitions(ComputeUnionBoundary PRIVATE CSV_IO_NO_THREAD)
target_link_libraries(ComputeUnionBoundary CGAL::CGAL ${PROJ_LIBRARY})

add_executable(ConvertGraph ConvertGraph.cc)
target_compile_definitions(ConvertGraph PRIVATE CSV_IO_NO_THREAD)
target_link_libraries(ConvertGraph ${PROJ_LIBRARY} ${RoutingKit_LIBRARY})

add_executable(DrawNetwork DrawNetwork.cc)
target_compile_definitions(DrawNetwork PRIVATE CSV_IO_NO_THREAD)
target_link_libraries(DrawNetwork ${cairo_LIBRARY})

add_executable(GenerateMatchingODPairs GenerateMatchingODPairs.cc)
target_link_libraries(GenerateMatchingODPairs ${RoutingKit_LIBRARY})

add_executable(GenerateODPairs GenerateODPairs.cc)
target_compile_definitions(GenerateODPairs PRIVATE CSV_IO_NO_THREAD)
if(OpenMP_FOUND)
  target_link_libraries(GenerateODPairs OpenMP::OpenMP_CXX)
endif()

add_executable(MeanDistancePerRank MeanDistancePerRank.cc)
target_compile_definitions(MeanDistancePerRank PRIVATE CSV_IO_NO_THREAD)
target_link_libraries(MeanDistancePerRank ${PROJ_LIBRARY})
if(OpenMP_FOUND)
  target_link_libraries(MeanDistancePerRank OpenMP::OpenMP_CXX)
endif()

add_executable(ParseMobiToppFile ParseMobiToppFile.cc)
target_compile_definitions(ParseMobiToppFile PRIVATE CSV_IO_NO_THREAD)

add_executable(ParseMtxFile ParseMtxFile.cc)
target_compile_definitions(ParseMtxFile PRIVATE CSV_IO_NO_THREAD)

add_executable(RadiationModel RadiationModel.cc)
target_compile_definitions(RadiationModel PRIVATE CSV_IO_NO_THREAD)
target_link_libraries(RadiationModel ${PROJ_LIBRARY})

add_executable(VisumExtractPoi VisumExtractPoi.cc)
target_compile_definitions(VisumExtractPoi PRIVATE CSV_IO_NO_THREAD)
target_link_libraries(VisumExtractPoi ${PROJ_LIBRARY})
