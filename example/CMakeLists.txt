include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

find_program(GNUPLOT_FOUND gnuplot)

add_executable(threshold_scan threshold_scan.cpp)
target_link_libraries(threshold_scan medipix)

add_executable(frequency_exposure frequency_pattern.cpp)
target_link_libraries(frequency_exposure medipix)

add_executable(pileup_scan pileup_scan.cpp)
target_link_libraries(pileup_scan medipix)

if (GNUPLOT_FOUND)
    set(GNUPLOT_INPUT_FILE_SPM \"input_file_spm='${CMAKE_CURRENT_BINARY_DIR}/threshold_scan_spm.txt'\")
    set(GNUPLOT_INPUT_FILE_CSM \"input_file_csm='${CMAKE_CURRENT_BINARY_DIR}/threshold_scan_csm.txt'\")

    add_custom_target(plot_threshold_scan
            COMMAND ${CMAKE_CURRENT_BINARY_DIR}/threshold_scan
            COMMAND gnuplot -e ${GNUPLOT_INPUT_FILE_SPM} -e ${GNUPLOT_INPUT_FILE_CSM} ${CMAKE_CURRENT_SOURCE_DIR}/gnuplot/plot_threshold_scan.gp
            DEPENDS threshold_scan
            BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/threshold_scan_spm.txt ${CMAKE_CURRENT_BINARY_DIR}/threshold_scan_csm.txt ${CMAKE_CURRENT_BINARY_DIR}/threshold_scan.pdf)

    set(GNUPLOT_IMAGE_INPUT_FILE \"input_file='${CMAKE_CURRENT_BINARY_DIR}/frequency_pattern_spm.raw'\")
    add_custom_target(plot_frequency_pattern
            COMMAND ${CMAKE_CURRENT_BINARY_DIR}/frequency_exposure
            COMMAND gnuplot -e ${GNUPLOT_IMAGE_INPUT_FILE} ${CMAKE_CURRENT_SOURCE_DIR}/gnuplot/plot_frequency_pattern.gp
            DEPENDS frequency_exposure
            BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/frequency_pattern_spm.raw ${CMAKE_CURRENT_BINARY_DIR}/frequency_pattern.pdf)


endif (GNUPLOT_FOUND)
