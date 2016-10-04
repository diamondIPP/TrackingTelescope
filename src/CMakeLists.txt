# source files
SET(LIB_SOURCE_FILES
        PLTAlignment.cc
        FileWriterTracking.cc)
# header files
SET(LIB_HEADER_FILES ${PROJECT_SOURCE_DIR}/include/PSIRootFileReader.h )
include_directories(../interface)
# root dictionary
#set(EVAL_DICT ${CMAKE_CURRENT_BINARY_DIR}/evalDict.cxx)

SET(MY_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} )

#ROOT_GENERATE_DICTIONARY_OLD("${LIB_HEADER_FILES}" "${MY_INCLUDE_DIRECTORIES}" "${EVAL_DICT}")
#root_generate_dictionary(${EVAL_DICT} LIB_HEADER_FILES LINKDEF inc/LinkDef.h)
#ADD_LIBRARY(${PROJECT_NAME} SHARED ${LIB_SOURCE_FILES} ${EVAL_DICT})
# Link necessary libraries:
#TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${ROOT_LIBRARIES})