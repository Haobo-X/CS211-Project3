module load mpich-3.2.1/gcc-4.8.5
cd build
cmake ..
make
cd ..
cd script
./submit.sh
