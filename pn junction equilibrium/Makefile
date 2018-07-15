



# Libraries required
PYT_LIB := $(shell pkg-config --cflags --libs python3)
LOC_LIB := $(shell pwd)
NPY_LIB := $(shell python -c "import numpy ; print(str(numpy.get_include()))")


# Create the wrappers and compile the library
compile:
	# eq library
	swig -python -debug-tmsearch eq.i ; rm eq.py ; sed -i 's/_eq/eq/g' eq_wrap.c
	gcc -fPIC -shared eq.c eq_wrap.c $(PYT_LIB) -o eq.so

# Clean the directory
clean:
	rm -r *.o ; rm -r *.out ; rm -rf __pycache__ ; rm *.so

# Test
test:
	python -c "import eq ; print(eq.equilibrium_Poisson_VonNeumann([1,2,3],[4,5,6])) ; print(eq.input_2x_PyList2D__output_2x_PyList2D([[1,2,3],[4,5,6]],[[1,2,3],[4,5,6]]))"
