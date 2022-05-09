include_directories(include)

add_executable(dft_test src/test/dft_test.cpp)

add_executable(lemma1 src/experiment/lemma1.cpp src/experiment/util.cpp)
add_executable(lemma2 src/experiment/lemma2.cpp src/experiment/util.cpp)
add_executable(theorem1 src/experiment/theorem1.cpp src/experiment/util.cpp)
add_executable(corollary src/experiment/corollary.cpp src/experiment/util.cpp)
