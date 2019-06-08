#ifndef MPI_PI_COMMON_INCLUDED
#define MPI_PI_COMMON_INCLUDED

#include <boost/mpi.hpp>

namespace mpi_pi::common {

struct Division {
  std::size_t first;
  std::size_t last;

  explicit Division(const boost::mpi::communicator &comm, std::size_t total) noexcept;

  std::size_t count() const noexcept;
};

} // namespace mpi_pi::common

#endif