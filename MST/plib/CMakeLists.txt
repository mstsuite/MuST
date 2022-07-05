target_sources(mst_lib
              PUBLIC
              GroupCommModule.F90
              MPPModule.F90
              ParallelIOModule.F90
              StatusModule.F90
              )


target_include_directories(mst_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
