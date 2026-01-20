#include "UStokes.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string  mesh_file_name  = "../mesh/mesh-pipe.msh";
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;

  UStokes problem(mesh_file_name, degree_velocity, degree_pressure);

  // Set time-stepping parameters
  problem.time_step = 0.01;
  problem.final_time = 1.0;

  problem.setup();
  problem.solve();

  return 0;
}