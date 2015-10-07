clear all;

% geometrical parameters of the domain:
l_x = 1;
l_y = 1;


% subdivisions of the domain:
n_el_x = 10;
n_el_y = 10;


% pde parameters:
laminar_viscosity = 0.1;
smagorinsky_coefficient = 0.1;
smagorinsky_caracteristic_length = max(l_x, l_y);
force_f = @(x) [ones(size(x, 1), 1), zeros(size(x, 1), 1)];


% build the mesh:
mesh = geometry.build_square_mesh(l_x, l_y, n_el_x, n_el_y, 0);


% solve the system:
[u_x, u_y, p] = navierstokes2d.solve(mesh, force_f, ...
				     laminar_viscosity, ...
				     smagorinsky_coefficient, ...
				     smagorinsky_caracteristic_length);


% visualize the result:
vis = 0;
if vis
  force_nodes = force_f(mesh.nodes);

  figure(1); hold on; cla;
  gplot(mesh.adj, mesh.nodes);
  h = quiver(mesh.nodes(:, 1), mesh.nodes(:, 2), u_x, u_y); set(h, 'color', 'r');
  h = quiver(mesh.nodes(:, 1), mesh.nodes(:, 2), force_nodes(:, 1), force_nodes(:, 2)); set(h, 'color', 'g');
  figure(2); cla;
  trisurf(mesh.elements, mesh.nodes(:, 1), mesh.nodes(:, 2), p);
end
