#ifort ./src/*.f ./src/*.f90 ./src/data/*.f90  ./src/hi_intgl/*.f90 ./src/property/*.f90 -I./build/module/ 
#mv *.mod ./build/module/
cd build
make
cd ..
