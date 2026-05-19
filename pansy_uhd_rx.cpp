//
// Resilient multi-channel UHD receiver for PANSY raw voltage recording.
//
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/thread.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <complex>
#include <thread>
#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <digital_rf.h>
#include <string>
#include <ctime>
#include <cmath>

namespace po = boost::program_options;

using namespace uhd::usrp;
using namespace std;

namespace {

void stop_rx_stream(const uhd::rx_streamer::sptr& rx_stream)
{
    if (!rx_stream) {
        return;
    }

    try {
        uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS);
        stream_cmd.stream_now = true;
        rx_stream->issue_stream_cmd(stream_cmd);
    } catch (const std::exception& e) {
        std::cerr << "Failed to stop RX stream cleanly: " << e.what() << std::endl;
    }
}

uint64_t time_spec_to_sample(const uhd::time_spec_t& ts, double rate)
{
    // Digital RF indexes samples as integer offsets from a global sample epoch.
    return static_cast<uint64_t>(std::llround(ts.get_real_secs() * rate));
}

uhd::time_spec_t next_stream_start_time(const multi_usrp::sptr& usrp, double delay)
{
    // Restart on a future integer second so a recovered channel keeps the same
    // time base as the other channels.
    const double now = usrp->get_time_now().get_real_secs();
    return uhd::time_spec_t(std::ceil(now + delay));
}

void assert_mboard_times_aligned(const multi_usrp::sptr& usrp)
{
    // Catch bad startups before recording. A sub-microsecond disagreement here
    // means interferometric phases will be wrong even if data is flowing.
    const size_t num_mboards = usrp->get_num_mboards();
    const double reference_time = usrp->get_time_last_pps(0).get_real_secs();
    bool aligned = true;

    std::cout << "Checking motherboard PPS time alignment." << std::endl;
    for (size_t mb = 0; mb < num_mboards; mb++) {
        const double mb_time = usrp->get_time_last_pps(mb).get_real_secs();
        const double delta = mb_time - reference_time;
        std::cout << boost::format(" * mboard %d last PPS: %.6f (delta %.6f s)")
                         % mb % mb_time % delta
                  << std::endl;
        if (std::abs(delta) > 0.5 / 1000000.0) {
            aligned = false;
        }
    }

    if (!aligned) {
        throw std::runtime_error(
            "USRP motherboard times are not aligned after PPS synchronization.");
    }
}

}

void get_usrp_time(multi_usrp::sptr usrp, size_t mboard, std::vector<int64_t>* times)
{
    (*times)[mboard] = usrp->get_time_now(mboard).get_full_secs();
}

void streaming_by_channel(size_t chan,double rate,std::string subdev,std::string outdir, multi_usrp::sptr usrp, uhd::time_spec_t time_last_pps)
{
    Digital_rf_write_object * data_object = NULL; /* main object created by init */
    uint64_t global_start_index; /* start sample (unix time * sample_rate) of first measurement - set below */
    int i, result;
    std::vector<size_t> channel_number;
    channel_number.push_back(chan);
    uint64_t sample_rate_numerator = static_cast<uint64_t>(std::llround(rate));
    uint64_t sample_rate_denominator = 1;
    uint64_t subdir_cadence = 3600;
    uint64_t millseconds_per_file = 1000; 
    int compression_level = 0; /* low level of compression */
    int checksum = 0; /* no checksum */
    int is_complex = 1; /* complex values */
    int is_continuous = 1; /* continuous data written */
    int num_subchannels = 1; /* only one subchannel */
    int marching_periods = 1; /* marching periods when writing */
    char uuid[100] = "Fake UUID - use a better one!";
    uint64_t vector_length = 363; /* one packet */

    // setup streaming
    double tstart=time_last_pps.get_real_secs()+2.0;
    uhd::time_spec_t ts_t0=uhd::time_spec_t(tstart);
    printf("Streaming start at %f\n",time_last_pps.get_real_secs()+2.0);

    // start recording at global_start_sample
    global_start_index = time_spec_to_sample(ts_t0, rate);
    printf("%lu",global_start_index);

    std::string ch_dir = outdir + "/ch" + std::string(3 - std::to_string(chan).length(), '0') + std::to_string(chan);

    std::cout << "Writing complex short to multiple files and subdirectores in " << ch_dir << std::endl;
    std::string mkdir_cmd = "mkdir -p "+ch_dir;
    std::cout << mkdir_cmd << std::endl;
    result = system(mkdir_cmd.c_str());

    data_object = digital_rf_create_write_hdf5((char *)ch_dir.c_str(),
					       H5T_NATIVE_SHORT,
					       subdir_cadence,
					       millseconds_per_file,
					       global_start_index,
					       sample_rate_numerator,
					       sample_rate_denominator,
					       uuid,
					       compression_level,
					       checksum,
					       is_complex,
					       num_subchannels,
					       is_continuous,
					       marching_periods);

    if (!data_object){
      printf("no data object created\n");
      exit(-1);
    }

    uint64_t prev_tl=0;
    size_t prev_num_rx_samps = 0;
    uhd::time_spec_t next_start_time = ts_t0;
    uint64_t restart_count = 0;
    const int max_empty_recvs = 10;
    const double restart_delay = 2.0;

    while (1)
    {
      uhd::rx_streamer::sptr rx_stream;

      try
      {
        // create a receive streamer for this thread's channel
        uhd::stream_args_t stream_args("sc16", "sc16"); // complex shorts
        stream_args.channels = channel_number;
        // Request one radar packet per recv when the transport supports it.
        stream_args.args["spp"] = std::to_string(vector_length);
        rx_stream = usrp->get_rx_stream(stream_args);

        // metadata
        uhd::rx_metadata_t md;

        // allocate buffer to receive with samples
        std::vector<std::complex<short>> buff(rx_stream->get_max_num_samps());
        std::vector<void*> buffs;
        buffs.push_back(&buff.front());

        // Always use a timed start, including after recovery. stream_now=true
        // would restart each channel at a slightly different device time.
        uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
        stream_cmd.stream_now = false;
        stream_cmd.time_spec = next_start_time;
        rx_stream->issue_stream_cmd(stream_cmd);

        std::cout << "Channel " << chan << " streaming start requested at "
                  << next_start_time.get_real_secs() << std::endl;

        // The first call after a timed stream start can block longer.
        double timeout = restart_delay + 1.0;
        int n_empty = 0;

        while (1)
        {
          size_t num_rx_samps = rx_stream->recv(buffs, buff.size(), md, timeout, true);

          if (num_rx_samps > 0 &&
              md.error_code == uhd::rx_metadata_t::ERROR_CODE_NONE) {
            n_empty = 0;
            uint64_t tl = time_spec_to_sample(md.time_spec, rate);

            if(prev_tl!=0)
            {
              const uint64_t expected_tl = prev_tl + prev_num_rx_samps;
              if (tl != expected_tl) {
                const int64_t gap = static_cast<int64_t>(tl) - static_cast<int64_t>(expected_tl);
                std::cerr << "ch " << chan << " timestamp discontinuity of "
                          << gap << " samples at " << md.time_spec.get_real_secs()
                          << std::endl;
              }
            }

            if (tl < global_start_index) {
              std::cerr << "ch " << chan << " received data before global start; skipping"
                        << std::endl;
              continue;
            }

            // Use the UHD timestamp for the Digital RF offset. This preserves
            // gaps/restarts as real sample-time gaps instead of compressing them.
            short *a = (short *)buff.data();
            result = digital_rf_write_hdf5(data_object, tl - global_start_index, a, num_rx_samps);
            if(result != 0) {
              stop_rx_stream(rx_stream);
              return;
            }
            prev_tl=tl;
            prev_num_rx_samps = num_rx_samps;
          }
          else
          {
            printf("ch %zu recv problem %d: samps=%zu error=%s\n",
                   chan,
                   n_empty,
                   num_rx_samps,
                   md.strerror().c_str());
            n_empty+=1;
            if(n_empty > max_empty_recvs)
            {
              throw std::runtime_error("channel stalled");
            }
          }

          // use a small timeout for subsequent packets
          timeout = 0.1;
        }
      }
      catch (const std::exception& e)
      {
        stop_rx_stream(rx_stream);
        restart_count += 1;
        // Recreate the stream and schedule the next attempt in the near future.
        next_start_time = next_stream_start_time(usrp, restart_delay);
        std::cerr << "Channel " << chan << " stalled after " << restart_count
                  << " restart attempts: " << e.what()
                  << ". Recreating stream at "
                  << next_start_time.get_real_secs() << "." << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
      }
    }
}

int UHD_SAFE_MAIN(int argc, char* argv[])
{
    // variables to be set by po
    std::string usrp_args;
    std::string clock_args;
    std::string outdir;
    std::string killdir;
    std::string wire;
    std::string subdev;
    uint32_t max_interval, num_tests;
    double seconds_in_future;
    size_t total_num_samps;
    double rate;
    std::string channel_list;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help", "help message")
        ("usrp_args", po::value<std::string>(&usrp_args)->default_value("addr0=192.168.11.2,addr1=192.168.12.2,addr2=192.168.13.2,addr3=192.168.14.2,recv_buff_size=500000000"),"ettus device args")
        ("clock_args",po::value<std::string>(&clock_args)->default_value("addr=192.168.10.3"),"octoclock address args")
        ("outdir", po::value<std::string>(&outdir)->default_value("/media/archive"), "output directory")
        ("killdir", po::value<std::string>(&killdir)->default_value("/home/sdr/chirpsounder2"), "kill directory")
        ("wire", po::value<std::string>(&wire)->default_value(""), "the over the wire type, sc16, sc8, etc")
        ("subdev", po::value<std::string>(&subdev)->default_value("A:A A:"), "subdevice")
        ("rate", po::value<double>(&rate)->default_value(1000000), "rate of incoming samples")
        ("dilv", "specify to disable inner-loop verbose")
        ("channels", po::value<std::string>(&channel_list)->default_value("0,1,2,3,4,5,6,7"), "which channel(s) to use (specify \"0\", \"1\", \"0,1\", etc)")
        ("max-interval", po::value<uint32_t>(&max_interval)->default_value(10000), "Maximum interval between comparisons (in ms)")
        ("num-tests", po::value<uint32_t>(&num_tests)->default_value(2), "Number of times to compare device times")
    ;
    
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    bool verbose = vm.count("dilv") == 0;

    
    // Create a Multi-USRP device
    std::cout << boost::format("\nCreating the USRP device with: %s") % usrp_args
              << std::endl;
    multi_usrp::sptr usrp = multi_usrp::make(usrp_args);
    // use external pps and 10 mhz ref
    usrp->set_clock_source("external");
    usrp->set_time_source("external");
    usrp->set_rx_subdev_spec(uhd::usrp::subdev_spec_t("A:A A:B"));
    usrp->set_rx_rate(rate);
    for (int i=0; i<8; i++){
      usrp->set_rx_freq(-47e6,i);
    }

    std::this_thread::sleep_for(std::chrono::seconds(2));

    std::cout << std::endl << "Checking USRP devices for lock." << std::endl;
    bool all_locked = true;
    for (size_t ch = 0; ch < usrp->get_num_mboards(); ch++) {
        std::string ref_locked = usrp->get_mboard_sensor("ref_locked", ch).value;
        std::cout << boost::format(" * %d: %s") % ch % ref_locked << std::endl;

        if (ref_locked != "true")
            all_locked = false;
    }
    if (not all_locked)
    {
        std::cout << std::endl << "WARNING: One or more devices not locked!" << std::endl;
        exit(0);
    }
    // detect which channels to use
    std::vector<std::string> channel_strings;
    std::vector<size_t> channel_nums;
    boost::split(channel_strings, channel_list, boost::is_any_of("\"',"));
    for (size_t ch = 0; ch < channel_strings.size(); ch++) {
        size_t chan = std::stoi(channel_strings[ch]);
	      std::cout << chan << std::endl;
        if (chan >= usrp->get_rx_num_channels()) 
        {
            throw std::runtime_error("Invalid channel(s) specified.");
        } else
            channel_nums.push_back(std::stoi(channel_strings[ch]));
    }
    std::cout << "done channels" << std::endl;

    // Get GPS time to initially set USRP devices
    std::cout << std::endl
              << "Querying Clock for time and setting USRP times..." << std::endl
              << std::endl;

    std::time_t now = std::time(nullptr);
    double timestamp = static_cast<double>(now);
    usrp->set_time_unknown_pps(uhd::time_spec_t(timestamp+2));

    // Wait for it to apply
    // The wait is 2 seconds because N-Series has a known issue where
    // the time at the last PPS does not properly update at the PPS edge
    // when the time is actually set.
    std::this_thread::sleep_for(std::chrono::seconds(2));

    assert_mboard_times_aligned(usrp);

    uhd::time_spec_t time_last_pps = usrp->get_time_last_pps();
    printf("USRP time now %1.4f USRP last pps %1.4f\n",usrp->get_time_now().get_real_secs(),time_last_pps.get_real_secs());

    // Threading for each channel
    std::vector<std::thread> threads;
    for (size_t ch = 0; ch < channel_strings.size(); ch++) {
        threads.push_back(std::thread(streaming_by_channel, std::stoi(channel_strings[ch]), rate, subdev, outdir, usrp, time_last_pps));
    }  
    
    // Join threads
    for(auto& thread : threads){
        thread.join();
    }

    return EXIT_SUCCESS;
}
