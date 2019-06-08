#include "../include/pi.hpp"
#include "../include/common.hpp"
#include "../include/random.hpp"
#include <boost/mpi.hpp>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>

namespace mpi = boost::mpi;
namespace mp = boost::multiprecision;

namespace mpi_pi {

namespace area_integral {

common::number pi(const mpi::communicator &comm, int root, std::size_t terms) {
  common::Division division(comm, terms);
  common::number nterms = terms;
  common::number delta = 1 / nterms;
  common::number local_sum = 0;
  common::number fxdx = 0;
  for (std::size_t i = division.first; i != division.last; ++i) {
    common::number ni = i;
    common::number x = ni / nterms;
    common::number x1 = (ni + 1) / nterms;
    common::number fxdx1 = delta / (x1 * x1 + 1);
    if (i == division.first) {
      fxdx = delta / (x * x + 1);
    }
    local_sum += (fxdx + fxdx1) / 2;
    fxdx = fxdx1;
  }
  local_sum *= 4;
  common::number result = 0;
  mpi::reduce(comm, local_sum, result, std::plus<common::number>(), root);
  return result;
}

} // namespace area_integral

namespace power_series {

common::number pi(const mpi::communicator &comm, int root, std::size_t terms) {
  common::Division division(comm, terms);
  common::number local_sum = 0;
  for (std::size_t i = division.first; i != division.last; ++i) {
    common::number ni = i;
    common::number term = 1 / (2 * ni + 1);
    if (i % 2 == 1)
      term = -term;
    local_sum += term;
  }
  local_sum *= 4;
  common::number result = 0;
  mpi::reduce(comm, local_sum, result, std::plus<common::number>(), root);
  return result;
}

} // namespace power_series

namespace improved_power_series {

common::number pi(const mpi::communicator &comm, int root, std::size_t terms) {
  using mp::pow;
  using std::pow;

  common::Division division(comm, terms);
  common::number local_sum = 0;
  for (std::size_t i = division.first; i != division.last; ++i) {
    common::number ni = i;
    common::number x = ni * 2 + 1;
    common::number term5 = 4 / (x * pow(5, x));
    common::number term239 = 1 / (x * pow(239, x));
    if (i % 2 == 1) {
      term5 = -term5;
      term239 = -term239;
    }
    local_sum += term5 - term239;
  }
  local_sum *= 4;
  common::number result = 0;
  mpi::reduce(comm, local_sum, result, std::plus<common::number>(), root);
  return result;
}

} // namespace improved_power_series

namespace monte_carlo {

common::number pi(const mpi::communicator &comm, int root, std::size_t terms) {
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
  mpi::reduce(comm, local_hit, total_hit, std::plus<std::size_t>(), root);
  common::number result = 0;
  if (comm.rank() == root) {
    result = common::number(total_hit) * 4 / terms;
  }
  return result;
}

} // namespace monte_carlo

namespace monte_carlo_integral {

common::number pi(const mpi::communicator &comm, int root, std::size_t terms) {
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
  mpi::reduce(comm, local_hit, total_hit, std::plus<std::size_t>(), root);
  common::number result = 0;
  if (comm.rank() == root) {
    result = common::number(total_hit) * 4 / terms;
  }
  return result;
}

} // namespace monte_carlo_integral

namespace random_integral {

common::number pi(const mpi::communicator &comm, int root, std::size_t terms) {
  common::Division division(comm, terms);
  random::MinimunStandardEngine random_engine(comm, root);
  std::uniform_real_distribution<double> dist(0, 1.0);
  common::number nterms = terms;
  common::number delta = 1 / nterms;
  common::number local_sum = 0;
  for (std::size_t i = division.first; i != division.last; ++i) {
    common::number x = dist(random_engine);
    common::number x1 = x + delta;
    common::number fxdx = delta / (x * x + 1);
    common::number fxdx1 = delta / (x1 * x1 + 1);
    local_sum += (fxdx + fxdx1) / 2;
  }
  local_sum *= 4;
  common::number result = 0;
  mpi::reduce(comm, local_sum, result, std::plus<common::number>(), root);
  return result;
}

} // namespace random_integral

namespace borwein1987 {

common::number pi(const boost::mpi::communicator &comm, int root,
                  std::size_t terms) {
  using mp::pow;
  using mp::sqrt;
  using std::pow;
  using std::sqrt;

  if (comm.rank() == root) {
    std::cerr << "borwein1987 is an iterative method, will not run parallelized"
              << std::endl;
  }
  common::number x = sqrt(common::number(2.0));
  common::number p = 2 + sqrt(common::number(2.0));
  common::number y = pow(common::number(2.0), common::number(1.0 / 4.0));
  for (std::size_t i = 0; i != terms; ++i) {
    x = 0.5 * (sqrt(x) + 1 / sqrt(x));
    p *= (x + 1) / (y + 1);
    y = (y * sqrt(x) + 1 / sqrt(x)) / (y + 1);
  }
  return p;
}

} // namespace borwein1987

namespace yasumasa2002 {

common::number pi(const mpi::communicator &comm, int root, std::size_t terms) {
  using mp::pow;
  using std::pow;

  common::Division division(comm, terms);
  common::number local_sum = 0;
  for (std::size_t i = division.first; i != division.last; ++i) {
    common::number ni = i;
    common::number x = ni * 2 + 1;
    common::number term49 = 12 / (x * pow(49, x));
    common::number term57 = 32 / (x * pow(57, x));
    common::number term239 = -5 / (x * pow(239, x));
    common::number term110443 = 12 / (x * pow(110443, x));
    if (i % 2 == 1) {
      term49 = -term49;
      term57 = -term57;
      term239 = -term239;
      term110443 = -term110443;
    }
    local_sum += term49 + term57 + term239 + term110443;
  }
  local_sum *= 4;
  common::number result = 0;
  mpi::reduce(comm, local_sum, result, std::plus<common::number>(), root);
  return result;
}

} // namespace yasumasa2002

} // namespace mpi_pi
