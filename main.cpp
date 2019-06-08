#include "config.hpp"
#include <boost/mpi.hpp>
#include <boost/program_options.hpp>
#include <iostream>

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

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);

    if (vm.count("help"))
      config.show_help = true;
  }

  if (config.show_help) {
    if (world.rank() == ROOT) {
      std::cout << description << std::endl;
    }
    return 0; // all processes exit
  }
}

} // namespace mpi_pi

int main(int argc, char **argv) { return mpi_pi::main(argc, argv); }