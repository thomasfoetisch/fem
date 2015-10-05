clear all;

% geometric dimensions of the domain:
l_x = 1;
l_y = 1;

% subdivisions of the domain:
n_el_x = 5;
n_el_y = 5;

% build a square mesh:
mesh = build_square_mesh(l_x, l_y, n_el_x, n_el_y, 0);

% pde parameters:
mu = 1.0;
force_f = @(x) [x(:, 2) > 0.0, zeros(size(x, 1), 1)];

[u_x, u_y, p] = stokes2d.solve(mesh, mu, force_f);

% visualisation of the result:
vis = 1;
if vis
  f_dof = force_f(mesh.nodes);
   
  close all;
  gplot(mesh.adj, mesh.nodes); hold on;
  quiver(mesh.nodes(:, 1), mesh.nodes(:, 2), u_x, u_y, 'r');
  quiver(mesh.nodes(:, 1), mesh.nodes(:, 2), f_dof(:, 1), f_dof(:, 2), 'g');
  figure(); trisurf(mesh.elements, mesh.nodes(:, 1), mesh.nodes(:, 2), p);
end
