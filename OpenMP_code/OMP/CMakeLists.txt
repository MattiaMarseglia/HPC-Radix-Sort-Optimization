cmake_minimum_required(VERSION 3.9)
project (Contest-OMP
    DESCRIPTION "RadixSort implementation with OpenMP"
	VERSION 0.1
	LANGUAGES C
)

SET (INCLUDE_LIBRARY_C RadixSort.c CACHE STRING "Library Name C")
SET (INCLUDE_LIBRARY_H RadixSort.h CACHE STRING "Library Name H")


#create Executable. The OX rappresents the opt option.
add_executable(${PROJECT_NAME}_O1 src/main.c src/${INCLUDE_LIBRARY_C} include/${INCLUDE_LIBRARY_H})
add_executable(${PROJECT_NAME}_O2 src/main.c src/${INCLUDE_LIBRARY_C} include/${INCLUDE_LIBRARY_H})
add_executable(${PROJECT_NAME}_O3 src/main.c src/${INCLUDE_LIBRARY_C} include/${INCLUDE_LIBRARY_H})
add_executable(${PROJECT_NAME}_seq_O0 src/main.c src/${INCLUDE_LIBRARY_C} include/${INCLUDE_LIBRARY_H})
add_executable(${PROJECT_NAME}_seq_O1 src/main.c src/${INCLUDE_LIBRARY_C} include/${INCLUDE_LIBRARY_H})
add_executable(${PROJECT_NAME}_seq_O2 src/main.c src/${INCLUDE_LIBRARY_C} include/${INCLUDE_LIBRARY_H})
add_executable(${PROJECT_NAME}_seq_O3 src/main.c src/${INCLUDE_LIBRARY_C} include/${INCLUDE_LIBRARY_H})



target_include_directories(${PROJECT_NAME}_O1 PRIVATE include)
target_include_directories(${PROJECT_NAME}_O2 PRIVATE include)
target_include_directories(${PROJECT_NAME}_O3 PRIVATE include)
target_include_directories(${PROJECT_NAME}_seq_O0 PRIVATE include)
target_include_directories(${PROJECT_NAME}_seq_O1 PRIVATE include)
target_include_directories(${PROJECT_NAME}_seq_O2 PRIVATE include)
target_include_directories(${PROJECT_NAME}_seq_O3 PRIVATE include)


#link OpenMP library for parallel program

find_package(OpenMP)
if (OpenMP_C_FOUND)
    target_link_libraries(${PROJECT_NAME}_O1 PRIVATE OpenMP::OpenMP_C)
    target_link_libraries(${PROJECT_NAME}_O2 PRIVATE OpenMP::OpenMP_C)
    target_link_libraries(${PROJECT_NAME}_O3 PRIVATE OpenMP::OpenMP_C)
endif()


#Set compile opt option.
target_compile_options(${PROJECT_NAME}_O1 PRIVATE -O1)
target_compile_options(${PROJECT_NAME}_O2 PRIVATE -O2)
target_compile_options(${PROJECT_NAME}_O3 PRIVATE -O3)
target_compile_options(${PROJECT_NAME}_seq_O0 PRIVATE -O0)
target_compile_options(${PROJECT_NAME}_seq_O1 PRIVATE -O1)
target_compile_options(${PROJECT_NAME}_seq_O2 PRIVATE -O2)
target_compile_options(${PROJECT_NAME}_seq_O3 PRIVATE -O3)



#set compiler features required when compiling a given target!
#Ci dice cosa ci serve per compilare codesto target
target_compile_features(${PROJECT_NAME}_O1 PRIVATE c_std_11)
target_compile_features(${PROJECT_NAME}_O2 PRIVATE c_std_11)
target_compile_features(${PROJECT_NAME}_O3 PRIVATE c_std_11)
target_compile_features(${PROJECT_NAME}_seq_O0 PRIVATE c_std_11)
target_compile_features(${PROJECT_NAME}_seq_O1 PRIVATE c_std_11)
target_compile_features(${PROJECT_NAME}_seq_O2 PRIVATE c_std_11)
target_compile_features(${PROJECT_NAME}_seq_O3 PRIVATE c_std_11)

#Changing TYPE_OF_ELEMENTS in all file
SET (TYPE_OF_ELEMENTS int CACHE STRING "Type of elements in the program")
target_compile_definitions(${PROJECT_NAME}_O1 PRIVATE -D TYPE_OF_ELEMENTS=${TYPE_OF_ELEMENTS})
target_compile_definitions(${PROJECT_NAME}_O2 PRIVATE -D TYPE_OF_ELEMENTS=${TYPE_OF_ELEMENTS})
target_compile_definitions(${PROJECT_NAME}_O3 PRIVATE -D TYPE_OF_ELEMENTS=${TYPE_OF_ELEMENTS})
target_compile_definitions(${PROJECT_NAME}_seq_O0 PRIVATE -D TYPE_OF_ELEMENTS=${TYPE_OF_ELEMENTS})
target_compile_definitions(${PROJECT_NAME}_seq_O1 PRIVATE -D TYPE_OF_ELEMENTS=${TYPE_OF_ELEMENTS})
target_compile_definitions(${PROJECT_NAME}_seq_O2 PRIVATE -D TYPE_OF_ELEMENTS=${TYPE_OF_ELEMENTS})
target_compile_definitions(${PROJECT_NAME}_seq_O3 PRIVATE -D TYPE_OF_ELEMENTS=${TYPE_OF_ELEMENTS})


# ---------------------------------------- TEST -----------------------------------------
#add_executable(${PROJECT_NAME}_test test/test.c )
#target_include_directories(${PROJECT_NAME}_test PRIVATE include)
#target_compile_definitions(${PROJECT_NAME}_test PRIVATE -D TYPE_OF_ELEMENTS=${TYPE_OF_ELEMENTS})

#add_test (NAME ${PROJECT_NAME}_test COMMAND ${PROJECT_NAME}_test)


# ------------------------------------ DOCUMENTATION ------------------------------------
#It's generates docs for our project
find_package(Doxygen
            REQUIRED dot)
#All file needs the docs 
if(DOXYGEN_FOUND)
	set(DOXYGEN_GENERATE_XML YES)
	doxygen_add_docs(
		docs
		./src/${INCLUDE_LIBRARY_C}
		)
endif()

# ---------------------------------- COMMAND  ----------------------------------
#With this command we execute the bash file
add_custom_target (
    first_command
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/measurements/analysis.bash ${CMAKE_CURRENT_BINARY_DIR} ${PROJECT_NAME}
    COMMENT "FIRST COMMAND EXECUTION [${CMAKE_CURRENT_SOURCE_DIR}/measurements/analisis.bash ${CMAKE_CURRENT_BINARY_DIR}] ${PROJECT_NAME}"
    DEPENDS	
        ${PROJECT_NAME}_O1
        ${PROJECT_NAME}_O2
        ${PROJECT_NAME}_O3
        ${PROJECT_NAME}_seq_O0
        ${PROJECT_NAME}_seq_O1
        ${PROJECT_NAME}_seq_O2
        ${PROJECT_NAME}_seq_O3
)

# --------------------------------- OTHER COMMAND  ---------------------------------
#Python3 for extract the data from analisis, PipEnv because we like it!
find_package(Python3 COMPONENTS Interpreter REQUIRED)
find_program(PipEnv pipenv REQUIRED)

add_custom_target(
	install_pipenv
	COMMAND PIPENV_PIPFILE=${CMAKE_CURRENT_SOURCE_DIR}/Pipfile ${PipEnv} install --ignore-pipfile
    COMMENT "Installing pipenv dependencies [${PipEnv} install]"
    DEPENDS	${CMAKE_CURRENT_SOURCE_DIR}/Pipfile.lock
)

add_custom_target(
    second_command
    COMMAND PIPENV_PIPFILE=${CMAKE_CURRENT_SOURCE_DIR}/Pipfile ${PipEnv} run python ${CMAKE_CURRENT_SOURCE_DIR}/measurements/extract.py
    COMMENT "Generating version header [${PipEnv} run python ${CMAKE_CURRENT_SOURCE_DIR}/measures/extract.py]"
    DEPENDS	
		install_pipenv
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/extract.py
)


