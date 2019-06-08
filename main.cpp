#include "config.hpp"
#include "include/common.hpp"
#include "include/pi.hpp"
#include <boost/mpi.hpp>
#include <boost/program_options.hpp>
#include <iomanip>
#include <iostream>
#include <string>

namespace mpi = boost::mpi;
namespace po = boost::program_options;

constexpr int ROOT = 0;

namespace mpi_pi {

int main(int argc, char **argv) {
  mpi::environment env(argc, argv, false);
  mpi::communicator world;

  config::Configuration config;
  po::options_description description("options");
  if (world.rank() == ROOT) {
    description.add_options()("help,h", "print help message");
    description.add_options()(
        "method,m", po::value<std::string>()->default_value("area_integral"),
        "method to calculate pi");
    description.add_options()("terms,t",
                              po::value<std::size_t>()->default_value(1000),
                              "terms number to calculate");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);

    if (vm.count("help"))
      config.show_help = true;
    if (vm.count("method"))
      config.method = vm["method"].as<std::string>();
    if (vm.count("terms"))
      config.terms = vm["terms"].as<std::size_t>();
  }

  mpi::broadcast(world, config, ROOT);

  if (config.show_help) {
    if (world.rank() == ROOT) {
      std::cout << description << std::endl;
    }
    return 0; // all processes exit
  }

  common::number result;
  if (config.method == "area_integral") {
    result = area_integral::pi(world, ROOT, config.terms);
  } else if (config.method == "power_series") {
    result = power_series::pi(world, ROOT, config.terms);
  } else if (config.method == "improved_power_series") {
    result = improved_power_series::pi(world, ROOT, config.terms);
  } else if (config.method == "monte_carlo") {
    result = monte_carlo::pi(world, ROOT, config.terms);
  } else if (config.method == "monte_carlo_integral") {
    result = monte_carlo_integral::pi(world, ROOT, config.terms);
  } else if (config.method == "random_integral") {
    result = random_integral::pi(world, ROOT, config.terms);
  } else if (config.method == "borwein1987") {
    result = borwein1987::pi(world, ROOT, config.terms);
  } else if (config.method == "yasumasa2002") {
    result = yasumasa2002::pi(world, ROOT, config.terms);
  } else {
    std::cout << "invalid method \"" << config.method << "\"" << std::endl;
  }
  if (world.rank() == ROOT) {
    std::cout << std::setprecision(
                     std::numeric_limits<common::number>::max_digits10)
              << result << std::endl;
  }
}

} // namespace mpi_pi

int main(int argc, char **argv) { return mpi_pi::main(argc, argv); }