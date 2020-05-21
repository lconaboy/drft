all:
	f2py -c cic.f95 -m cic

strict:
	f2py -c cic.f95 -m cic --f90flags='-fbounds-check'

.PHONY : clean
clean :
	rm -f cic.cpython*
