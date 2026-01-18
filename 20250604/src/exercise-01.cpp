#include <iostream>

#include "ADR1D.hpp"

// Main function.
int
main(int /*argc*/, char * /*argv*/[])
{
  constexpr unsigned int dim = ADR1D::dim;

  const unsigned int N_el = 40;
  const unsigned int r    = 1;
  const auto         mu   = [](const Point<dim>           &/*p*/) { return 1.0; };
  const auto         f    = [](const Point<dim> &p) {
    if (p[0] <= 0.125 || p[0] > 0.25)
      return 0.0;
    else
      return -1.0;
  };
  const auto         b    = [](const Point<dim> & /*p*/) { return 0.0; };
  const auto         sigma= [](const Point<dim> & /*p*/) { return 0.0; };

  ADR1D problem(N_el, r, mu, b, sigma, f);

  problem.setup();
  problem.assemble();
  problem.solve();
  problem.output();

  return 0;
}
