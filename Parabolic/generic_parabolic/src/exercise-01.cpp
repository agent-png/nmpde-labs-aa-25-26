#include "Parabolic.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  constexpr unsigned int dim = Parabolic::dim;

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const auto mu = [](const Point<dim> & /*p*/) { return 0.1; };
  const auto b=[](const Point<dim> & /*p*/){
    Tensor<1,dim> v;
    v[0]=0.0;
    v[1]=0.0;
    v[2]=0.0;
    return v;
  };
  const auto sigma = [](const Point<dim> & /*p*/) { return 0.0; };
  const auto f  = [](const Point<dim>  &/*p*/, const double  &/*t*/) {
    return 0.0;
  };

  //Initial condition in .hpp
  
  Parabolic problem(/*mesh_filename = */ "../mesh/mesh-cube-10.msh",
               /* degree = */ 1,
               /* T = */ 1.0,
               /* theta = */ 0.0,
               /* delta_t = */ 0.0025,
               mu,
               b,
               sigma,
               f);

  problem.run();

  return 0;
}