cmake_minimum_required(VERSION 2.6)
project(gen3d)

enable_language(Fortran)

#set (CMAKE_FORTRAN_COMPILER "/usr/bin/gfortran")
set (CMAKE_FORTRAN_COMPILER "ifort")

set(FOPT "-w -O2 -Wno-dev -cpp -fimplicit-none -shared-intel")

set(CMAKE_FORTRAN_FLAGS "${FOPT}")

##################################
# include directories
##################################
include_directories(
src 
/home/s.engkadac/mylibs/parmetis-4.0.3-install/include
/app/libraries/impi/5.0.1.035/include64
/home/s.engkadac/mylibs/metis-5.1.0-install/include
)


##################################
# lib directories
##################################
link_directories(
#/usr/lib/gcc/x86_64-linux-gnu/7
#/usr/lib/x86_64-linux-gnu
/usr/lib 
#/usr/x86_64-linux-gnu/lib 
#/opt/mpi/openmpi-icc/2.0.2.10/lib
/app/libraries/impi/5.0.1.035/lib64
/home/s.engkadac/mylibs/parmetis-4.0.3-install/lib
/home/s.engkadac/mylibs/metis-5.1.0-install/lib
/app/libraries/petsc/3.7.5/el6/AVX/intel-16.0/intel-5.1/lib
)

##################################
# source files
##################################


file(GLOB_RECURSE sources1 src/*.f90)


add_executable(partfort src/meshpartitionfort.f90 ${sources1})

target_link_libraries(partfort parmetis metis mpi mpiif)

install(TARGETS partfort RUNTIME DESTINATION /home/s.engkadac/projects/partition/bin)









