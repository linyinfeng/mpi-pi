#ifndef MPI_PI_PI_INCLUDED
#define MPI_PI_PI_INCLUDED

#include "../include/common.hpp"
#include "../include/random.hpp"
#include <boost/multiprecision/number.hpp>
#include <boost/mpi.hpp>

namespace mpi_pi {

namespace area_integral {

template <typename Number>
Number pi(const boost::mpi::communicator &comm, int root, std::size_t terms) {
  common::Division division(comm, terms);
  Number nterms = terms;
  Number delta = 1 / nterms;
  Number local_sum = 0;
  Number fxdx = 0;
  for (std::size_t i = division.first; i != division.last; ++i) {
    Number ni = i;
    Number x = ni / nterms;
    Number x1 = (ni + 1) / nterms;
    Number fxdx1 = delta / (x1 * x1 + 1);
    if (i == division.first) {
      fxdx = delta / (x * x + 1);
    }
    local_sum += (fxdx + fxdx1) / 2;
    fxdx = fxdx1;
  }
  local_sum *= 4;
  Number result = 0;
  boost::mpi::reduce(comm, local_sum, result, std::plus<Number>(), root);
  return result;
}

} // namespace area_integral

namespace power_series {

template <typename Number>
Number pi(const boost::mpi::communicator &comm, int root, std::size_t terms) {
  common::Division division(comm, terms);
  Number local_sum = 0;
  for (std::size_t i = division.first; i != division.last; ++i) {
    Number ni = i;
    Number term = 1 / (2 * ni + 1);
    if (i % 2 == 1)
      term = -term;
    local_sum += term;
  }
  local_sum *= 4;
  Number result = 0;
  boost::mpi::reduce(comm, local_sum, result, std::plus<Number>(), root);
  return result;
}

} // namespace power_series

namespace improved_power_series {

template <typename Number>
Number pi(const boost::mpi::communicator &comm, int root, std::size_t terms) {
  using boost::multiprecision::pow;
  using std::pow;

  common::Division division(comm, terms);
  Number local_sum = 0;
  for (std::size_t i = division.first; i != division.last; ++i) {
    Number ni = i;
    Number x = ni * 2 + 1;
    Number term5 = 4 / (x * pow(5, x));
    Number term239 = 1 / (x * pow(239, x));
    if (i % 2 == 1) {
      term5 = -term5;
      term239 = -term239;
    }
    local_sum += term5 - term239;
  }
  local_sum *= 4;
  Number result = 0;
  boost::mpi::reduce(comm, local_sum, result, std::plus<Number>(), root);
  return result;
}

} // namespace improved_power_series

namespace monte_carlo {

template <typename Number>
Number pi(const boost::mpi::communicator &comm, int root, std::size_t terms) {
  common::Division division(comm, terms);
  random::MinimunStandardEngine random_engine(comm, root);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  std::size_t local_hit = 0;
  for (std::size_t i = division.first; i != division.last; ++i) {
    auto x = dist(random_engine);
    auto y = dist(random_engine);
    if ((x * x + y * y) <= 1)
      local_hit += 1;
  }
  std::size_t total_hit = 0;
  boost::mpi::reduce(comm, local_hit, total_hit, std::plus<std::size_t>(), root);
  Number result = 0;
  if (comm.rank() == root) {
    result = Number(total_hit) * 4 / terms;
  }
  return result;
}

} // namespace monte_carlo

namespace monte_carlo_integral {

template <typename Number>
Number pi(const boost::mpi::communicator &comm, int root, std::size_t terms) {
  common::Division division(comm, terms);
  random::MinimunStandardEngine random_engine(comm, root);
  std::uniform_real_distribution<double> dist(0, 1.0);
  std::size_t local_hit = 0;
  for (std::size_t i = division.first; i != division.last; ++i) {
    auto x = dist(random_engine);
    auto y = dist(random_engine);
    double fx = 1.0 / (1.0 + x * x);
    if (y <= fx)
      local_hit += 1;
  }
  std::size_t total_hit = 0;
  boost::mpi::reduce(comm, local_hit, total_hit, std::plus<std::size_t>(), root);
  Number result = 0;
  if (comm.rank() == root) {
    result = Number(total_hit) * 4 / terms;
  }
  return result;
}

} // namespace monte_carlo_integral

namespace random_integral {

template <typename Number>
Number pi(const boost::mpi::communicator &comm, int root, std::size_t terms) {
  common::Division division(comm, terms);
  random::MinimunStandardEngine random_engine(comm, root);
  std::uniform_real_distribution<double> dist(0, 1.0);
  Number nterms = terms;
  Number delta = 1 / nterms;
  Number local_sum = 0;
  for (std::size_t i = division.first; i != division.last; ++i) {
    Number x = dist(random_engine);
    Number x1 = x + delta;
    Number fxdx = delta / (x * x + 1);
    Number fxdx1 = delta / (x1 * x1 + 1);
    local_sum += (fxdx + fxdx1) / 2;
  }
  local_sum *= 4;
  Number result = 0;
  boost::mpi::reduce(comm, local_sum, result, std::plus<Number>(), root);
  return result;
}

} // namespace random_integral

namespace borwein1987 {

template <typename Number>
Number pi(const boost::mpi::communicator &comm, int root, std::size_t terms) {
  using boost::multiprecision::pow;
  using boost::multiprecision::sqrt;
  using std::pow;
  using std::sqrt;

  if (comm.rank() == root) {
    std::cerr << "borwein1987 is an iterative method, will not run parallelized"
              << std::endl;
  }
  Number x = sqrt(Number(2.0));
  Number p = 2 + sqrt(Number(2.0));
  Number y = pow(Number(2.0), Number(1.0 / 4.0));
  for (std::size_t i = 0; i != terms; ++i) {
    x = 0.5 * (sqrt(x) + 1 / sqrt(x));
    p *= (x + 1) / (y + 1);
    y = (y * sqrt(x) + 1 / sqrt(x)) / (y + 1);
  }
  return p;
}

} // namespace borwein1987

namespace yasumasa2002 {

template <typename Number>
Number pi(const boost::mpi::communicator &comm, int root, std::size_t terms) {
  using boost::multiprecision::pow;
  using std::pow;

  common::Division division(comm, terms);
  Number local_sum = 0;
  for (std::size_t i = division.first; i != division.last; ++i) {
    Number ni = i;
    Number x = ni * 2 + 1;
    Number term49 = 12 / (x * pow(49, x));
    Number term57 = 32 / (x * pow(57, x));
    Number term239 = -5 / (x * pow(239, x));
    Number term110443 = 12 / (x * pow(110443, x));
    if (i % 2 == 1) {
      term49 = -term49;
      term57 = -term57;
      term239 = -term239;
      term110443 = -term110443;
    }
    local_sum += term49 + term57 + term239 + term110443;
  }
  local_sum *= 4;
  Number result = 0;
  boost::mpi::reduce(comm, local_sum, result, std::plus<Number>(), root);
  return result;
}

} // namespace yasumasa2002

} // namespace mpi_pi

#endif