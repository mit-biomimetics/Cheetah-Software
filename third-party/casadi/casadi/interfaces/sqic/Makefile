all: libwsqic.a

libwsqic.a: wsqic.f90
	gfortran -c -O  -fPIC -J /home/jg/programs/sqic/mod -I /home/jg/programs/sqic/mod -o wsqic.o wsqic.f90
	ar rcu libwsqic.a wsqic.o

#prog: wsqic.cpp wsqic.hpp libwsqic.a
#	gcc wsqic.cpp -o prog -L. -lstdc++ -lwsqic /home/jg/programs/sqic-install/lib/libsqic.a /home/jg/programs/sqic-install/lib/libsnblas.a -lgfortran
