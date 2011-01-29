# Makefile for SOG unit test suite

tests: fruit_driver
	./fruit_driver

fruit_driver:	fruit_util.f90 fruit.f90 test_*.f90 fruit_driver.f90
	g95 -c fruit_util.f90 fruit.f90
	g95 -c test_*.f90 fruit_driver.f90 -I../..
	g95 -o fruit_driver *.o ../../libSOG.a

clean:
	rm test*.o fruit*.o
	rm test*.mod fruit*.mod
	rm fruit_driver
