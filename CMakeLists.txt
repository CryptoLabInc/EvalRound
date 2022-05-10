include_directories(include)

add_executable(dft_test src/test/dft_test.cpp)
add_executable(crt_test src/test/crt_test.cpp)
add_executable(ntt_test src/test/ntt_test.cpp)
add_executable(conv_test src/test/conv_test.cpp)

add_executable(lemma1 src/experiment/lemma1.cpp)
add_executable(lemma2 src/experiment/lemma2.cpp)
add_executable(theorem1 src/experiment/theorem1.cpp)
add_executable(corollary src/experiment/corollary.cpp)
