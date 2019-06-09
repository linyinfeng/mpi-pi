#include "../include/common.hpp"
#include <boost/mpi.hpp>
#include <cstddef>
#include <algorithm>
#include <cmath>

namespace mpi_pi::common {

Division::Division(const boost::mpi::communicator &comm, std::size_t total) noexcept {
  const auto rank = static_cast<std::size_t>(comm.rank());
  const auto size = static_cast<std::size_t>(comm.size());

  const auto basic_local_count = static_cast<std::size_t>(std::ceil(static_cast<double>(total) / size));
  this->first = basic_local_count * rank;
  this->last = std::min(total, basic_local_count * (rank + 1));
}

std::size_t Division::count() const noexcept { return this->last - this->first; }

} // namespace mpi_pi::common
