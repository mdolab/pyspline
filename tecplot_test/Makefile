default: intel

intel:
	ifort -c -O2 -fPIC tec_test.f90
	ifort -lstdc++  tec_test.o tecio.a -o tec_test
clean:
	-rm tec_test
	-rm tec_test.o