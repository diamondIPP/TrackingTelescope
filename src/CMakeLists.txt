INCLUDE_DIRECTORIES( . ${PROJECT_SOURCE_DIR}/interface )

# source files
SET(LIB_SOURCE_FILES
        PLTTesterReader.cc
        PLTTelescope.cc
        PLTAlignment.cc
        PLTHistReader.cc
        PLTAnalysis.cc
        DoAlignment.cc
        PSIBinaryFileReader.cc
        GetNames.cc
        PLTTrack.cc
        PLTTextFileReader.cc
        PLTPlane.cc
        PSIGainInterpolator.cc
        PLTTracking.cc
        PLTHardwareMap.cc
        PLTHit.cc
        TestPlaneEfficiencySilicon.cc
        PLTU.cc
        RootItems.cc
        FileWriterTracking.cc
        PLTGainCal.cc
        PLTBinaryFileReader.cc
        PLTDACs.cc
        PSIFileReader.cc
        PLTCluster.cc
        PSIRootFileReader.cc
        PLTTesterEvent.cc
        Utils.cc)

# header files
SET(LIB_HEADER_FILES
        ${PROJECT_SOURCE_DIR}/interface/PSIRootFileReader.h)

# root dictionary
set(EVAL_DICT ${CMAKE_CURRENT_BINARY_DIR}/evalDict.cxx)

SET(MY_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/interface ${ROOT_INCLUDE_DIR} )

ROOT_GENERATE_DICTIONARY_OLD("${LIB_HEADER_FILES}" "${MY_INCLUDE_DIRECTORIES}" "${EVAL_DICT}")
#root_generate_dictionary(${EVAL_DICT} LIB_HEADER_FILES LINKDEF inc/LinkDef.h)
ADD_LIBRARY(${PROJECT_NAME} SHARED ${LIB_SOURCE_FILES} ${EVAL_DICT})
# Link necessary libraries:
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${ROOT_LIBRARIES})

INSTALL(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)