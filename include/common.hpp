#ifndef MPI_PI_COMMON_INCLUDED
#define MPI_PI_COMMON_INCLUDED

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/mpi.hpp>

namespace mp = boost::multiprecision;

namespace mpi_pi::common {
    
constexpr unsigned PRECISION = 100;
//using number = mp::number<mp::cpp_dec_float<PRECISION>>;
using number = double;

struct Division {
  std::size_t first;
  std::size_t last;

  explicit Division(const boost::mpi::communicator &comm, std::size_t total) noexcept;

  std::size_t count() const noexcept;
};

} // namespace mpi_pi::common

#endif