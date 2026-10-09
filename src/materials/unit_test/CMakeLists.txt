add_subdirectory (twoPoint_SLS_Relaxation/)
add_subdirectory (utPeridigm_ElasticPlastic/)

# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string (COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)
if (NOT ${buildInSource})
   set(UT_DIRS ${UT_DIRS}  -d src/materials/unit_test/twoPoint_SLS_Relaxation)
   set(UT_DIRS ${UT_DIRS}  -d src/materials/unit_test/utPeridigm_ElasticPlastic)
   execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/scripts/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${UT_DIRS} )
endif (NOT ${buildInSource})

add_executable(utPeridigm_ElasticMaterial ./utPeridigm_ElasticMaterial.cpp
#is the ) correct here?
)
target_link_libraries(utPeridigm_ElasticMaterial
  ${Peridigm_LIBRARY}
  ${PdMaterialUtilitiesLib}
  PdField
  QuickGrid
  ${REQUIRED_LIBS}
  ${Trilinos_LIBRARIES}
)
add_test (utPeridigm_ElasticMaterial python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_ElasticMaterial)


add_executable(utPeridigm_MultiphysicsElasticMaterial ./utPeridigm_MultiphysicsElasticMaterial.cpp)
target_link_libraries(utPeridigm_MultiphysicsElasticMaterial
  ${Peridigm_LIBRARY}
  ${PdMaterialUtilitiesLib}
  PdField
  QuickGrid
  ${REQUIRED_LIBS}
  ${Trilinos_LIBRARIES}
)
add_test (utPeridigm_MultiphysicsElasticMaterial python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_MultiphysicsElasticMaterial)

add_executable(utJH2MaterialPoint ./utJH2MaterialPoint.cpp)
target_link_libraries(utJH2MaterialPoint
  ${PdMaterialUtilitiesLib}
  ${REQUIRED_LIBS}
  ${Trilinos_LIBRARIES}
)
add_test(utJH2MaterialPoint python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utJH2MaterialPoint)

add_executable(utPeridigm_JH2CorrespondenceMaterial
  ./utPeridigm_JH2CorrespondenceMaterial.cpp)
target_link_libraries(utPeridigm_JH2CorrespondenceMaterial
  ${Peridigm_LIBRARY}
  PdField
  ${PDNEIGH_LIBS}
  QuickGrid
  ${UTILITIES_LIBS}
  ${PdMaterialUtilitiesLib}
  ${REQUIRED_LIBS}
  ${Trilinos_LIBRARIES}
)
add_test(utPeridigm_JH2CorrespondenceMaterial
  python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py
  ./utPeridigm_JH2CorrespondenceMaterial)
