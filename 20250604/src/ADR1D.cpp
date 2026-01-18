#include "ADR1D.hpp"

void
ADR1D::setup()
{
  std::cout << "===============================================" << std::endl;

  // Create the mesh.
  {
    std::cout << "Initializing the mesh" << std::endl;
    GridGenerator::subdivided_hyper_cube(mesh, N_el, 0.0, 1.0, true);
    std::cout << "  Number of elements = " << mesh.n_active_cells()
              << std::endl;

    // Write the mesh to file.
    //
    // Since we generate the mesh internally, we also write it to file for
    // possible inspection by the user. This would not be necessary if we read
    // the mesh from file, as we will do later on.
    const std::string mesh_file_name = "mesh-" + std::to_string(N_el) + ".vtk";
    GridOut           grid_out;
    std::ofstream     grid_out_file(mesh_file_name);
    grid_out.write_vtk(mesh, grid_out_file);
    std::cout << "  Mesh saved to " << mesh_file_name << std::endl;
  }

  std::cout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    std::cout << "Initializing the finite element space" << std::endl;

    fe = std::make_unique<FE_SimplexP<dim>>(r);

    std::cout << "  Degree                     = " << fe->degree << std::endl;
    std::cout << "  DoFs per cell              = " << fe->dofs_per_cell
              << std::endl;

    // Construct the quadrature formula of the appopriate degree of exactness.
    // This formula integrates exactly the mass matrix terms (i.e. products of
    // basis functions).
    quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);

    std::cout << "  Quadrature points per cell = " << quadrature->size()
              << std::endl;
  }

  std::cout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    std::cout << "Initializing the DoF handler" << std::endl;

    // Initialize the DoF handler with the mesh we constructed.
    dof_handler.reinit(mesh);

    // Distribute the degrees of freedom.
    dof_handler.distribute_dofs(*fe);

    std::cout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  std::cout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    std::cout << "Initializing the linear system" << std::endl;

    std::cout << "  Initializing the sparsity pattern" << std::endl;
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    // Then, we use the sparsity pattern to initialize the system matrix
    std::cout << "  Initializing the system matrix" << std::endl;
    system_matrix.reinit(sparsity_pattern);

    // Finally, we initialize the right-hand side and solution vectors.
    std::cout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(dof_handler.n_dofs());
    std::cout << "  Initializing the solution vector" << std::endl;
    solution.reinit(dof_handler.n_dofs());
  }
}

void
ADR1D::assemble()
{
  std::cout << "===============================================" << std::endl;

  std::cout << "  Assembling the linear system" << std::endl;

  // Number of local DoFs for each element.
  const unsigned int dofs_per_cell = fe->dofs_per_cell;

  // Number of quadrature points for each element.
  const unsigned int n_q = quadrature->size();

  // FEValues instance. This object allows to compute basis functions, their
  // derivatives, the reference-to-current element mapping and its
  // derivatives on all quadrature points of all elements.
  FEValues<dim> fe_values(
    *fe,
    *quadrature,
    // Here we specify what quantities we need FEValues to compute on
    // quadrature points. For our test, we need:
    // - the values of shape functions (update_values);
    // - the derivative of shape functions (update_gradients);
    // - the position of quadrature points (update_quadrature_points);
    // - the quadrature weights (update_JxW_values).
    update_values | update_gradients /*| update_hessians*/ | update_quadrature_points |
      update_JxW_values);

  // Local matrix and right-hand side vector. We will overwrite them for
  // each element within the loop.
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  // We will use this vector to store the global indices of the DoFs of the
  // current element within the loop.
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  // Reset the global matrix and vector, just in case.
  system_matrix = 0.0;
  system_rhs    = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      // Reinitialize the FEValues object on current element. This
      // precomputes all the quantities we requested for all quadrature nodes of
      // the current cell.
      fe_values.reinit(cell);

      // We reset the cell matrix and vector.
      cell_matrix = 0.0;
      cell_rhs    = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        {
          // Here we assemble the local contribution for current cell and
          // current quadrature point, filling the local matrix and vector.
          const double mu_loc = mu(fe_values.quadrature_point(q));
          const double b_loc = b(fe_values.quadrature_point(q));
          const double sigma_loc = sigma(fe_values.quadrature_point(q));
          const double f_loc  = f(fe_values.quadrature_point(q));

          // Here we iterate over *local* DoF indices.
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  // Diffusion term.
                  cell_matrix(i, j) += mu_loc *                     //
                                       fe_values.shape_grad(i, q) * //
                                       fe_values.shape_grad(j, q) * //
                                       fe_values.JxW(q);
                  // Transport term.
                  cell_matrix(i, j) += b_loc *                          //
                                       fe_values.shape_value(i, q)    * //
                                       fe_values.shape_grad(j, q)[0]  * //
                                       fe_values.JxW(q);
                  // Reaction term.
                  cell_matrix(i, j) += sigma_loc *
                                       fe_values.shape_value(i, q) * //
                                       fe_values.shape_value(j, q) * //
                                       fe_values.JxW(q);
                }
              // Forcing funciton f contribution.
              cell_rhs(i) += f_loc *                       //
                             fe_values.shape_value(i, q) * //
                             fe_values.JxW(q);
            }
        }

      // Neumann boundary conditions.
      /*
      gamma  = [](const Point<dim> &p) { return p[0]; //Replace accordingly!! };
      if (cell->at_boundary())
      {
        // We loop over its edges.
        for (unsigned int face_number = 0; face_number < cell->n_faces(); ++face_number)
        {
          // If current face lies on the boundary, and its boundary ID (or
          // tag) is that of one of the Neumann boundaries, we assemble the
          // boundary integral.
          if (cell->face(face_number)->at_boundary() &&
          //replace conditions on boundary_id as necessary!!!
            (cell->face(face_number)->boundary_id() == 0 ||
              cell->face(face_number)->boundary_id() == 1))
          {
            fe_values_boundary.reinit(cell, face_number);

            for (unsigned int q = 0; q < quadrature_boundary->size(); ++q)
            {
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                cell_rhs(i) +=
                  gamma(fe_values_boundary.quadrature_point(q)) * //
                  fe_values_boundary.shape_value(i, q) *      //
                  fe_values_boundary.JxW(q);
            }
          }
        }
      }
      */

      // Retrieve the global indices of the DoFs of current cell.
      cell->get_dof_indices(dof_indices);

      // Add the local matrix and vector into the corresponding
      // positions of the global matrix and vector.
      system_matrix.add(dof_indices, cell_matrix);
      system_rhs.add(dof_indices, cell_rhs);
    }

  // Boundary conditions (Dirichlet).
  {
    // We construct a map that stores, for each DoF corresponding to a Dirichlet
    // condition, the corresponding value. E.g., if the Dirichlet condition is
    // u_i = b_i, the map will contain the pair (i, b_i).
    std::map<types::global_dof_index, double> boundary_values;

    // Dirichlet function. Some ways of having it (zero, constant, custom)
    Functions::ZeroFunction<dim> zero_function;
    Functions::ConstantFunction<dim> one_function(/*constant = */ 1.0);
    class G_Function : public Function<dim>
    {
    public:
      G_Function() : Function<dim>(1){} // Constructor.

      virtual double value(const Point<dim> &p, const unsigned int component = 0) const override
      {
        return p[0]; // just an example, replace accordingly.
      }
    };
    G_Function g;

    // Then, we build a map that, for each boundary tag, stores a pointer to the
    // corresponding boundary function.
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;
    boundary_functions[0] = &zero_function;
    boundary_functions[1] = &zero_function;

    // interpolate_boundary_values fills the boundary_values map.
    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values);

    // Finally, we modify the linear system to apply the boundary conditions.
    // This replaces the equations for the boundary DoFs with the corresponding
    // u_i = 0 equations.
    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, solution, system_rhs, true);
  }
}

void
ADR1D::solve()
{
  std::cout << "===============================================" << std::endl;

  // Here we specify the maximum number of iterations of the iterative solver,
  // and its absolute and relative tolerances.
  ReductionControl solver_control(/* maxiter = */ 1000,
                                  /* tolerance = */ 1.0e-16,
                                  /* reduce = */ 1.0e-6);

  // Since the system matrix is symmetric and positive definite, we solve the
  // system using the conjugate gradient method.
  SolverGMRES<Vector<double>> solver(solver_control);

  std::cout << "  Solving the linear system" << std::endl;
  // We use the identity preconditioner for now.
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
  std::cout << "  " << solver_control.last_step() << " GMRES iterations"
            << std::endl;
}

void
ADR1D::output() const
{
  std::cout << "===============================================" << std::endl;

  // The DataOut class manages writing the results to a file.
  DataOut<dim> data_out;

  // It can write multiple variables (defined on the same mesh) to a single
  // file. Each of them can be added by calling add_data_vector, passing the
  // associated DoFHandler and a name.
  data_out.add_data_vector(dof_handler, solution, "solution");

  // Once all vectors have been inserted, call build_patches to finalize the
  // DataOut object, preparing it for writing to file.
  data_out.build_patches();

  // Then, use one of the many write_* methods to write the file in an
  // appropriate format.
  const std::string output_file_name =
    "output-" + std::to_string(N_el) + ".vtk";
  std::ofstream output_file(output_file_name);
  data_out.write_vtk(output_file);

  std::cout << "Output written to " << output_file_name << std::endl;

  std::cout << "===============================================" << std::endl;
}

double
ADR1D::compute_error(const VectorTools::NormType &norm_type,
                         const Function<dim>         &exact_solution) const
{
  // The error is an integral, and we approximate that integral using a
  // quadrature formula. To make sure we are accurate enough, we use a
  // quadrature formula with one node more than what we used in assembly.
  const QGaussSimplex<dim> quadrature_error(r + 2);

  // First we compute the norm on each element, and store it in a vector.
  Vector<double> error_per_cell(mesh.n_active_cells());
  VectorTools::integrate_difference(dof_handler,
                                    solution,
                                    exact_solution,
                                    error_per_cell,
                                    quadrature_error,
                                    norm_type);

  // Then, we add out all the cells.
  const double error =
    VectorTools::compute_global_error(mesh, error_per_cell, norm_type);

  return error;
}