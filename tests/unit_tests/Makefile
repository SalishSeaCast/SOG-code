# Makefile for SOG unit test suite

test: fruit fruit_driver
	./fruit_driver

fruit_driver:	test_*.f90 fruit_driver.f90
	g95 -c test_*.f90 fruit_driver.f90 -I../..
	g95 -o fruit_driver *.o ../../libSOG.a

fruit:	fruit_util.f90 fruit.f90
	g95 -c fruit_util.f90 fruit.f90

clean:
	rm test_*.o test_*.mod
	rm fruit_driver.o fruit_driver
