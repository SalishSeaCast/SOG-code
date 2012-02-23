# Makefile for SOG unit test suite

F90=g95
LD=g95

tests: fruit_driver
	./fruit_driver

fruit_driver:	fruit_util.f90 fruit.f90 test_*.f90 fruit_driver.f90
	$(F90) -c fruit_util.f90 fruit.f90
	$(F90) -c test_*.f90 fruit_driver.f90 -I../..
	$(LD) -o fruit_driver *.o ../../libSOG.a

clean:
	rm test*.o fruit*.o
	rm test*.mod fruit*.mod
	rm fruit_driver
