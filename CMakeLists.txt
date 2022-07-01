
cmake_minimum_required(VERSION 3.14)
project(dcc012-lab4-heap)

option(CLASSROOM "Enable classroom build" OFF)
option(CPU_TIME  "Enable cputime build" ON)

set(INPUT_DIR "${PROJECT_SOURCE_DIR}/tests/input/")
set(EXPECTED_OUTPUT_DIR "${PROJECT_SOURCE_DIR}/tests/expected_output/")

configure_file(${PROJECT_SOURCE_DIR}/siga/config.h.in ${PROJECT_SOURCE_DIR}/siga/include/config.h) 
include_directories(${PROJECT_SOURCE_DIR}/siga/include)


enable_testing()
add_subdirectory(siga)
add_subdirectory(tests)

# teste 1
add_test(
       NAME Setup
    COMMAND cmake -E remove -f dataset.bin && cmake -E touch dataset.bin
)

# teste 2
add_test(
  NAME   TestInsert
  COMMAND $<TARGET_FILE:test1>
  )

#   teste 3
add_test(
   NAME  TestFindMax
   COMMAND $<TARGET_FILE:test2> 
   )

#    teste 4
add_test(
    NAME  TestRemoveMax
    COMMAND $<TARGET_FILE:test3> 
)

# teste 5
add_test(
    NAME  TestMerge
    COMMAND $<TARGET_FILE:test4> 
)

# teste 6
add_test(
    NAME  TestFindRemoveMaxSiga
    COMMAND $<TARGET_FILE:test5> 
)

# teste 7
add_test(
    NAME  TestHeapMergeSiga
    COMMAND $<TARGET_FILE:test6> 
)


# add_test(
#     NAME  MergeSortCompare
#     COMMAND cmake -E compare_files test2.csv ${PROJECT_SOURCE_DIR}/tests/expected_output/test2.csv
#     )

# add_test(
#   NAME  TestQuickSort1000
#   COMMAND COMMAND $<TARGET_FILE:test3> 
# )

# add_test(
#     NAME  QuickSortCompare
#     COMMAND cmake -E compare_files test2.csv ${PROJECT_SOURCE_DIR}/tests/expected_output/test3.csv
#     )

# add_test(
#   NAME  Importe2000
#   COMMAND COMMAND $<TARGET_FILE:test4> 
# )

# add_test(
#   NAME  OrdenaCursoENome
#   COMMAND COMMAND $<TARGET_FILE:test5> 
# )

# add_test(
#     NAME  OrdenaCursoENomeCompare
#     COMMAND cmake -E compare_files test5.csv ${PROJECT_SOURCE_DIR}/tests/expected_output/test5.csv
#     )

# add_test(
#       NAME  Compara
#       COMMAND COMMAND $<TARGET_FILE:test6> 
# )
    