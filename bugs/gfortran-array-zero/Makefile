
GFORTRAN_HOME=/home/wdobler/opt/misc/gfortran

test: test.x run

test.x: test.f90
	${GFORTRAN_HOME}/bin/gfortran test.f90 -o test.x

run: test.x
	env LD_LIBRARY_PATH=${GFORTRAN_HOME}/lib64 test.x

clean:
	rm -f *.o *.x a.out

