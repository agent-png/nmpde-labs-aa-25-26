#include <iostream>
#include <deal.II/base/convergence_table.h>

#include "ADR1D.hpp"
constexpr unsigned int dim = ADR1D::dim;

// Exact solution.
class ExactSolution : public Function<dim>
{
public:
  // Constructor.
  ExactSolution()
  {}

  // Evaluation.
  virtual double
  value(const Point<dim> &p,
        const unsigned int /*component*/ = 0) const override
  {
    return p[0] - p[0]*p[0];
  }

  // Gradient evaluation.
  virtual Tensor<1, dim>
  gradient(const Point<dim> &p,
           const unsigned int /*component*/ = 0) const override
  {
    Tensor<1, dim> result;
    result[0] = 1 - 2*p[0];
    return result;
  }
};

// Main function.
int
main(int /*argc*/, char * /*argv*/[])
{
  ConvergenceTable table;

  const std::vector<unsigned int> N_el_values = {5, 10, 20, 40};
  const unsigned int              r           = 1;

  const auto         mu   = [](const Point<dim> &/*p*/) { return 1.0; };
  const auto         b    = [](const Point<dim> &p) { return -p[0]; };
  const auto         sigma= [](const Point<dim> & /*p*/) { return 1.0; };
  const auto         f    = [](const Point<dim> &p) { return p[0]*p[0] + 2.0; };

  const ExactSolution exact_solution;

  std::ofstream convergence_file("convergence.csv");
  convergence_file << "h,eL2,eH1" << std::endl;

  for (const auto &N_el : N_el_values)
    {
      ADR1D problem(N_el, r, mu, b, sigma, f);

      problem.setup();
      problem.assemble();
      problem.solve();
      problem.output();

      const double h = 1.0 / N_el;

      const double error_L2 =
        problem.compute_error(VectorTools::L2_norm, exact_solution);
      const double error_H1 =
        problem.compute_error(VectorTools::H1_norm, exact_solution);

      table.add_value("h", h);
      table.add_value("L2", error_L2);
      table.add_value("H1", error_H1);

      convergence_file << h << "," << error_L2 << "," << error_H1 << std::endl;
    }

  table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);

  table.set_scientific("L2", true);
  table.set_scientific("H1", true);

  table.write_text(std::cout);


  return 0;
}
