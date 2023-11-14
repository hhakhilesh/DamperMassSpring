# DamperMassSpring
1. Establish a sym link between numpy installation and /usr/include/numpy. Check:
https://stackoverflow.com/questions/44888925/fatal-error-numpy-arrayobject-h-no-such-file-or-directory/44935933
2. Compile with
   g++ -o main main.cxx -I/usr/include/numpy -I/usr/include/python3.8 -lpython3.8
