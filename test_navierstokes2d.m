clear all;
more off;

% geometrical parameters of the domain:
l_x = 1;
l_y = 1;


% subdivisions of the domain:
n_el_x = 10;
n_el_y = 10;


% pde parameters:
ctx = struct();
ctx.rho = 1;
ctx.laminar_viscosity = 0.01;
ctx.smagorinsky_coefficient = 0.01;
ctx.smagorinsky_caracteristic_length = max(l_x, l_y);
force_f = @(x) 0.1 * [x(:, 2) > 0.25, zeros(size(x, 1), 1)];


% build the mesh:
mesh = geometry.build_square_mesh(l_x, l_y, n_el_x, n_el_y, 0);


% solve the system:
[u_x, u_y, p, newton_errors] = navierstokes2d.solve(mesh, force_f, ctx);


% visualize the result:
vis = 1;
if vis
  force_nodes = force_f(mesh.nodes);

  figure(1); hold on; cla;
  gplot(mesh.adj, mesh.nodes);
  h = quiver(mesh.nodes(:, 1), mesh.nodes(:, 2), u_x, u_y); set(h, 'color', 'r');
  h = quiver(mesh.nodes(:, 1), mesh.nodes(:, 2), force_nodes(:, 1), force_nodes(:, 2)); set(h, 'color', 'g');
  figure(2); cla;
  trisurf(mesh.elements, mesh.nodes(:, 1), mesh.nodes(:, 2), p);
end
