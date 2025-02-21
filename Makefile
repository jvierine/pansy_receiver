all: test_gps pansy_uhd_rx

test_gps: test_gps.cc
	g++ -I/home/radar/.syowa/include  test_gps.cc -L/home/radar/.syowa/lib -o test_gps -luhd

#`pkg-config --cflags hdf5 digital_rf`
#`pkg-config --libs hdf5 digital_rf`
pansy_uhd_rx:
	g++ `pkg-config --cflags hdf5 digital_rf` -I/home/radar/.syowa/include -o pansy_uhd_rx pansy_uhd_rx.cpp -pthread  -L/home/radar/.syowa/lib -lboost_program_options -lboost_system -lboost_thread -lboost_date_time -lboost_regex -lboost_serialization -ldigital_rf `pkg-config --libs hdf5 digital_rf`
