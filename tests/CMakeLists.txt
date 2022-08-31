

# excutable
#add_subdirectory(tests)
link_libraries(splocs yaml-cpp)
add_executable(neko_test "neko_test.cpp" "../eigen.natvis")
add_executable(eigen_test "eigen_test.cpp" "../eigen.natvis")
add_executable(face_test "face_test.cpp" "../eigen.natvis")
add_executable(face_test2 "face_test2.cpp" "../eigen.natvis")
