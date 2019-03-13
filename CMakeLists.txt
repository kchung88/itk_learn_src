cmake_minimum_required(VERSION 3.9.5)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
 
add_executable(rtk_vtk MACOSX_BUNDLE rtk_vtk.cxx)
target_link_libraries(rtk_vtk
  ${Glue}  ${VTK_LIBRARIES} ${ITK_LIBRARIES})


