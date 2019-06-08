#ifndef MPI_PI_CONFIG_INCLUDED
#define MPI_PI_CONFIG_INCLUDED

#include <boost/serialization/access.hpp>
#include <boost/serialization/nvp.hpp>
#include <cstddef>
#include <string>

namespace mpi_pi::config {

struct Configuration {
public:
  bool show_help;
  std::string method;
  std::size_t terms;

  Configuration() { show_help = false; }
  
private:
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar &BOOST_SERIALIZATION_NVP(show_help);
    ar &BOOST_SERIALIZATION_NVP(method);
    ar &BOOST_SERIALIZATION_NVP(terms);
  }
};

} // namespace mpi_pi::config

#endif
