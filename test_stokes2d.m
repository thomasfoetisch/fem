clear all;
more off;
graphics_toolkit('gnuplot');

% geometric dimensions of the domain:
l_x = 1;
l_y = 1;

% subdivisions of the domain:
n_el_x = 11;
n_el_y = 11;

% build a square mesh:
mesh = geometry.build_square_mesh(l_x, l_y, n_el_x, n_el_y, 0*pi/6);


% define useful partition of the boundary nodes:
all_nodes_id = 1:length(mesh.boundary_nodes);
corner_nodes_id = find(ismember(mesh.boundary_nodes, mesh.corner_nodes));
face_nodes_id = setdiff(all_nodes_id, corner_nodes_id);


% define dirichlet boundary conditions:
boundary_conditions = struct();
boundary_conditions.u_n = {face_nodes_id, zeros(length(all_nodes_id), 1)}; %{[], zeros(length(corner_nodes_id), 1)};
boundary_conditions.u_t = {[], zeros(length(all_nodes_id), 1)}; %{[], zeros(length(corner_nodes_id), 1)};
boundary_conditions.u_x = {corner_nodes_id, zeros(length(corner_nodes_id), 1)}; %{[corner_nodes_id; face_nodes_id'], @(x) zeros(size(x, 1), 1)};
boundary_conditions.u_y = {corner_nodes_id, zeros(length(corner_nodes_id), 1)}; %{corner_nodes_id, @(x) zeros(size(x, 1), 1)};


% pde parameters:
mu = 1.0;
force_f = @(x) [x(:, 2) > 0.0, zeros(size(x, 1), 1)];


% solve the problem and remove the bubbles:
[u_x, u_y, p] = stokes2d.solve(mesh, mu, force_f, boundary_conditions);
u_x = u_x(1:size(mesh.nodes, 1));
u_y = u_y(1:size(mesh.nodes, 1));

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
