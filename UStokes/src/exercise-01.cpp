#include "UStokes.hpp"

int main(int argc, char *argv[]){
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  constexpr unsigned int dim = UStokes::dim;

  const std::string  mesh_file_name  = "../mesh/mesh-pipe.msh";
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;
  const auto mu = [](const Point<dim> & /*p*/) { return 1.0; };
  const auto f  = [](const Point<dim>  &/*p*/, const double  &/*t*/) {
    return 0.0;
  };
  const double alpha = 100.0;

  UStokes problem(mesh_file_name,
                  degree_velocity,
                  degree_pressure,
                  1.0,
                  1.0,
                  0.0025,
                  alpha,
                  mu,
                  f);

  problem.run();

  return 0;
}