all:
	python -m numpy.f2py -c cic.f95 -m cic

.PHONY : clean
clean :
	rm cic.cpython*
