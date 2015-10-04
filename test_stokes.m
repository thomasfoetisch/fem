clear all;

l_x = 1;
l_y = 1;


n_el_x = 5;
n_el_y = 5;

% pde parameters:
mu = 1.0;
%force_f = @(x) [ones(size(x, 1), 1), zeros(size(x, 1), 1)];
force_f = @(x) [x(:, 2) > 0.0, zeros(size(x, 1), 1)];


% elementary integrals:
ints = elementary_integrals_p1bulle();

% build a square mesh:
mesh = build_square_mesh(l_x, l_y, n_el_x, n_el_y, 0);


% build the matrix:
n_nodes = size(mesh.nodes, 1);
n_elems = size(mesh.elements, 1);
n_u_dof = n_nodes + n_elems;
n_dof = n_u_dof * 2 + n_nodes;

mat = spalloc(n_dof, n_dof, 0);


% build the degree of freedom map to take the bubbles into account:
dof_map = [mesh.elements, ((n_nodes + 1):n_u_dof)'];

tic();
for el = 1:n_elems
  % elliptic term:
  elem_u = zeros(4, 2, 4, 2);
  for sol_dim = 1:2
    for test_dim = 1:2
      if sol_dim == test_dim
	for i = 1:4
	  for j = 1:4
	    for k = 1:2
	      for p = 1:2
		for q = 1:2
		  elem_u(i, sol_dim, j, test_dim) ...
		  = elem_u(i, sol_dim, j, test_dim) ...
		    + mu * mesh.jac(el) ...
		      * mesh.jmt(k, p, el) ...
		      * mesh.jmt(k, q, el) ...
		      * ints.dphidphi(p, i, q, j);
		end
	      end
	    end
	  end
	end
      end
    end
  end

  % pressure term:
  elem_p = zeros(3, 4, 2);
  for test_dim = 1:2 
    for i = 1:3
      for j = 1:4
	for p = 1:2
	  elem_p(i, j, test_dim) ...
	  = elem_p(i, j, test_dim) ...
	    - mesh.jac(el) ...
	      * mesh.jmt(test_dim, p, el) ...
	      * ints.phidphi(i, p, j);
	end
      end
    end
  end

  % incompressibility constraint:
  elem_div = zeros(4, 2, 3);
  for sol_dim = 1:2
    for j = 1:3
      for i = 1:4
	for p = 1:2
	  elem_div(i, sol_dim, j) ...
	  = elem_div(i, sol_dim, j) ...
	    + mesh.jac(el) ...
	      * mesh.jmt(sol_dim, p, el) ...
	      * ints.phidphi(j, p, i);
	end
      end
    end
  end


  % accumulate everything in the linear system:
  for i = 1:4
    for j = 1:4
      for sol_dim = 1:2
	for test_dim = 1:2
	  if test_dim == sol_dim
	    mat((test_dim - 1) * n_u_dof + dof_map(el, j),...
		(sol_dim - 1) * n_u_dof + dof_map(el, i)) ...
	    = mat((test_dim - 1) * n_u_dof + dof_map(el, j), ...
		  (sol_dim - 1) * n_u_dof + dof_map(el, i)) ...
	      + elem_u(i, sol_dim, j, test_dim);
	  end
	end
      end
    end
  end

  for i = 1:3
    for j = 1:4
      for test_dim = 1:2
	mat((test_dim - 1) * n_u_dof + dof_map(el, j),
	    2 * n_u_dof + dof_map(el, i)) ...
	= mat((test_dim - 1) * n_u_dof + dof_map(el, j),
	      2 * n_u_dof + dof_map(el, i)) ...
	  + elem_p(i, j, test_dim);
      end
    end
  end

  for i = 1:4
    for sol_dim = 1:2
      for j = 1:3
	  mat(2 * n_u_dof + dof_map(el, j),
	      (sol_dim - 1) * n_u_dof +  dof_map(el, i)) ...
	  = mat(2 * n_u_dof + dof_map(el, j),
		(sol_dim - 1) * n_u_dof +  dof_map(el, i)) ...
	  + elem_div(i, sol_dim, j);
      end
    end
  end
end

build_time = toc();
printf('per element build time: %f\n', build_time / n_elems);


% build the rhs:
rhs = spalloc(size(mat, 1), 1, 2 * n_u_dof);
f_dofs = force_f(mesh.nodes);
for el = 1:n_elems
  for sol_dim = 1:2
    for i = 1:4
      for k = 1:3 % no bubble for the force field
	rhs((sol_dim - 1) * n_u_dof + dof_map(el, i), 1) ... 
	= rhs((sol_dim - 1) * n_u_dof + dof_map(el, i), 1) ...
	  + mesh.jac(el) * f_dofs(dof_map(el, k), sol_dim) * ints.phiphi(i, k);
      end
    end
  end
end

% set the dirichlet boundary conditions for the velocity:
for sol_dim = 1:2
  for n = 1:length(mesh.boundary_nodes)
    mat((sol_dim - 1) * n_u_dof + mesh.boundary_nodes(n), :) = 0;
    mat((sol_dim - 1) * n_u_dof + mesh.boundary_nodes(n), (sol_dim - 1) * n_u_dof + mesh.boundary_nodes(n)) = 1;
    rhs((sol_dim - 1) * n_u_dof + mesh.boundary_nodes(n), 1) = 0;
  end
end

% set the dirichlet boundary conditions for the pressure:
mean_pressure = 0;
if mean_pressure
  mean_p = spalloc(1, n_dof + 1, n_nodes + 1);
  for el = 1:n_elems
    for i = 1:3
      mean_p(1, 2 * n_u_dof + dof_map(el, i)) = mean_p(1, 2 * n_u_dof + dof_map(el, i)) + mesh.jac(el) * ints.phi(i);
    end
  end
  mean_p(1, end) = 1;
  mat = [mat, mean_p(1:end-1)'; mean_p];
  rhs = [rhs; 0];
else
  mat(end, :) = 0;
  mat(end, end) = 1;
end

% solve the system and extract the different variables:
x = mat \ rhs;
u_x = x(1:n_nodes);
u_y = x((n_u_dof + 1):(n_u_dof + n_nodes));
p = x((2 * n_u_dof + 1): 2 * n_u_dof + n_nodes);


% visualisation of the result:
vis = 1;
if vis
  close all;
  gplot(mesh.adj, mesh.nodes, 'color', 0.5 * [1 1 1]); hold on;
  quiver(mesh.nodes(:, 1), mesh.nodes(:, 2), u_x, u_y, 'r');
  quiver(mesh.nodes(:, 1), mesh.nodes(:, 2), f_dofs(:, 1), f_dofs(:, 2), 'g');
  figure(); trisurf(mesh.elements, mesh.nodes(:, 1), mesh.nodes(:, 2), p);
end
