function [u_x, u_y, p, newton_errors] = solve(mesh, force_f, ctx)

% problem sizes:
n_nodes = size(mesh.nodes, 1);
n_elems = size(mesh.elements, 1);
n_u_dof = n_nodes + n_elems;
n_p_dof = n_nodes;
n_dof = 2 * n_u_dof + n_p_dof;

n_max_newton_iterations = 20;

% build the degree of freedom mapping between the elements and the linear system:
dof_map = [mesh.elements, (n_nodes + 1:n_u_dof)'];


% get the elementary integrals:
[ints, basis] = functional.elementary_integrals_p1bulle();


% initialise the solution with a guess to start the newton iteration:
x = navierstokes2d.build_newton_initial_solution(mesh);
u_x = x(1:n_u_dof);
u_y = x(n_u_dof + (1:n_u_dof));
p = x(2 * n_u_dof + (1:n_nodes));

% newton iteration:
newton_errors = [];
for k = 1:n_max_newton_iterations
  mat = navierstokes2d.assemble_p1bullep1_stationary_matrix(mesh, dof_map, ints, basis, ...
							    u_x, u_y, ...
							    ctx);
  for j = 1:1
    rhs = navierstokes2d.assemble_p1bullep1_stationary_rhs(mesh, dof_map, ints, basis, ...
							   u_x, u_y, p, ...
							   force_f, ...
							   ctx);

    real_u_dof = setdiff(dof_map(:), mesh.boundary_nodes);
    real_p_dof = 1:(n_nodes - 1);
    rhs_dofs = rhs([real_u_dof', n_u_dof + real_u_dof' , 2*n_u_dof + real_p_dof]);
    rhs_norm_l2 = sqrt(sum(rhs_dofs.^2));
    printf('rhs l2 norm: %g\n', rhs_norm_l2);
    newton_errors(end + 1) = rhs_norm_l2;

    [mat_bc, rhs_bc, T_inv] = navierstokes2d.set_boundary_conditions(mesh, mat, rhs);
    delta_x = mat_bc \ rhs_bc;
    delta_x = T_inv * delta_x;
    
    printf('delta_x l2 norm: %g\n', sqrt(sum(delta_x.^2)));

    x = x - delta_x(1:(2 * n_u_dof + n_p_dof));

    u_x = x(1:n_u_dof);
    u_y = x(n_u_dof + (1:n_u_dof));
    p = x(2 * n_u_dof + (1:n_nodes));

    vis = 1;
    if vis
      figure(1); hold on; cla;
      gplot(mesh.adj, mesh.nodes);
      h = quiver(mesh.nodes(:, 1), mesh.nodes(:, 2), u_x(1:n_nodes), u_y(1:n_nodes)); set(h, 'color', 'r');
      figure(2); cla;
      trisurf(mesh.elements, mesh.nodes(:, 1), mesh.nodes(:, 2), p);
      figure(3); cla
      semilogy(newton_errors);
    end
  end

  if max(abs(delta_x)) < 1e-10
    break;
  end

end

u_x = u_x(1:n_nodes);
u_y = u_y(1:n_nodes);
