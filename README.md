# mpi-pi

Implementations of several pi calculation algorithms in MPI and Boost.Multiprecision.

## Algorithms

- `area_integral`: Area integral to calculate arctan(1) - arctan(0)
- `power_series`: Power series of arctan(1)
- `improved_power_series`: pi / 4 = 4arctan(1 / 5) - arctan(1 / 239)
- `monte_carlo`: Monte carlo method with circle in square
- `monte_carlo_integral`: Monte carlo method with area integral
- `random_integral`: Integral with random x value
- `borwein1987`: *Pi and the AGM: A Study in Analytic Number Theory and Computational Complexity* by Jonathan M. Borwein, Peter B. Borwein
- `yasumasa2002`: 2002 world record by 金田康正
- `chudnovsky`: Chudnovsky algorithm
- `bbp`: Bailey–Borwein–Plouffe formula

## Usage

```plaintext
options:
  -h [ --help ]                        print help message
  -m [ --method ] arg (=area_integral) method to calculate pi
  -t [ --terms ] arg (=1000)           terms number to calculate
```
