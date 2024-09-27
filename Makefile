all: test_gps


test_gps: test_gps.cc
	g++ -I/home/radar/.syowa/include  test_gps.cc -L/home/radar/.syowa/lib -o test_gps -luhd
