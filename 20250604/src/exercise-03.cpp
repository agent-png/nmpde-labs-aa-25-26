#include <iostream>

#include "ADR2D.hpp"

// Main function.
int
main(int /*argc*/, char * /*argv*/[])
{
  constexpr unsigned int dim = ADR2D::dim;

  const std::string  mesh_file_name = "../mesh/mesh-square-20.msh";
  const unsigned int r              = 1;

  const auto mu = [](const Point<dim> & /*p*/) { return 1.0; };
  auto b = [](const Point<dim> & /*p*/){
    Tensor<1, dim> v;
    v[0] = 0.0;
    v[1] = 0.0;
    return v;
  };

  // Reaction
  auto sigma = [](const Point<dim> & /*p*/) { return 1.0; };
  const auto f  = [](const Point<dim>  &p) { return (20.0 * M_PI * M_PI + 1.0) * std::sin(2.0 * M_PI * p[0]) *
           std::sin(4.0 * M_PI * p[1]); };
  const auto gamma  = [](const Point<dim> & /*p*/) { return 0.0; };

  ADR2D problem(mesh_file_name, r, mu, b, sigma, f, gamma);

  problem.setup();
  problem.assemble();
  problem.solve();
  problem.output();

  return 0;
}