#ifndef MPI_PI_PI_INCLUDED
#define MPI_PI_PI_INCLUDED

#include "../include/common.hpp"

namespace mpi_pi {

namespace area_integral {
extern common::number pi(const boost::mpi::communicator &comm, int root,
                         std::size_t terms);
} // namespace area_integral

namespace power_series {
extern common::number pi(const boost::mpi::communicator &comm, int root,
                         std::size_t terms);
} // namespace power_series

namespace improved_power_series {
extern common::number pi(const boost::mpi::communicator &comm, int root,
                         std::size_t terms);
} // namespace imporoved_power_series

namespace monte_carlo {
extern common::number pi(const boost::mpi::communicator &comm, int root,
                         std::size_t terms);
}

namespace monte_carlo_integral {
extern common::number pi(const boost::mpi::communicator &comm, int root,
                         std::size_t terms);
}

namespace random_integral {
extern common::number pi(const boost::mpi::communicator &comm, int root,
                         std::size_t terms);
}

namespace borwein1987 {
extern common::number pi(const boost::mpi::communicator &comm, int root,
                         std::size_t terms);
}

namespace yasumasa2002 {
extern common::number pi(const boost::mpi::communicator &comm, int root,
                         std::size_t terms);
}

} // namespace mpi_pi

#endif