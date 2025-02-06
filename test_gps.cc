#include <stdio.h>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <boost/algorithm/string.hpp> //for split
#include <boost/program_options.hpp>
#include <iostream>
#include <sstream>
#include <uhd/usrp_clock/multi_usrp_clock.hpp>
#include <uhd/types/device_addr.hpp>


int main(int argc, char **argv)
{
  uhd::device_addr_t dev;
  dev["addr"] = "192.168.10.3";
  uhd::usrp_clock::multi_usrp_clock::sptr clock = uhd::usrp_clock::multi_usrp_clock::make(dev);
  /*  device_addr_t dev;
  dev["addr"] = 192.168.10.3;
  //multi_usrp::sptr usrp = multi_usrp::make(usrp_args);
  uhd::usrp_clock::multi_usrp_clock::sptr clock = uhd::usrp_clock::multi_usrp_clock::make(dev);*/
  const auto sensor_value = clock->get_sensor("gps_detected");
  std::cout << sensor_value.to_pp_string() <<  std::endl;

  const auto sensor_value2 = clock->get_sensor("gps_locked");
  std::cout << sensor_value2.to_pp_string() <<  std::endl;

  const auto sensor_value3 = clock->get_sensor("gps_time");
  std::cout << sensor_value3.to_pp_string() <<  std::endl;  

  const auto sensor_value4 = clock->get_sensor("gps_gpgga");
  std::cout << sensor_value4.to_pp_string() <<  std::endl;  

  const auto sensor_value5 = clock->get_sensor("gps_gprmc");
  std::cout << sensor_value5.to_pp_string() <<  std::endl;  

  const auto sensor_value6 = clock->get_sensor("gps_servo");
  std::cout << sensor_value6.to_pp_string() <<  std::endl;  


}
