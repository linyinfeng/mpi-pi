#ifndef MPI_RANDOM_RANDOM_HPP
#define MPI_RANDOM_RANDOM_HPP

#include <boost/mpi.hpp>
#include <cstdint>
#include <memory>
#include <random>

namespace mpi_pi::random {

template <typename UIntType, UIntType m> class Engine;

using MinimunStandardEngine = Engine<std::uint_fast32_t, 2147483647>;

// random engine for parallel random number generating
// the engine must be construct at the same time on all mpi process
template <typename UIntType, UIntType m> class Engine {
public:
  using result_type = UIntType;

  Engine(const boost::mpi::communicator &comm, int root, UIntType seed) {
    const UIntType origin_a = 48271;
    const UIntType origin_c = 0;

    const int size = comm.size();

    if (comm.rank() == root) {
      std::vector<UIntType> xs;
      xs.reserve(size);
      xs.push_back(seed);
      this->a = 1;
      this->c = 0;

      for (int i = 0; i < size; ++i) {
        if (i != 0)
          xs.push_back(run(origin_a, origin_c, xs[i - 1]));
        this->c += this->a;
        this->c %= m;
        this->a *= origin_a;
        this->a %= m;
      }
      this->c *= origin_c;

      boost::mpi::scatter(comm, xs, this->current, root);
    } else {
      boost::mpi::scatter(comm, this->current, root);
    }
    boost::mpi::broadcast(comm, this->a, root);
    boost::mpi::broadcast(comm, this->c, root);
  }

  Engine(const boost::mpi::communicator &comm, int root)
      : Engine(comm, root, new_seed_in_root(comm, root)) {}

  static constexpr result_type min() { return 1; }

  static constexpr result_type max() { return m - 1; }

  result_type operator()() {
    this->current = run(this->a, this->c, this->current);
    return this->current;
  }

private:
  static UIntType new_seed_in_root(const boost::mpi::communicator &comm,
                                   int root) {
    if (comm.rank() == root) {
      std::random_device rd;
      std::uniform_int_distribution<UIntType> dist(min(), max());
      return dist(rd);
    } else {
      return 0;
    }
  }

  UIntType run(UIntType pa, UIntType pc, UIntType x) {
    return (pa * x + pc) % m;
  }

  UIntType a;
  UIntType c;

  UIntType current;
};

} // namespace mpi_pi::random

#endif