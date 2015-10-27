clear all;
more off;
graphics_toolkit('gnuplot');

% geometric dimensions of the domain:
l_x = 1;
l_y = 8;

% subdivisions of the domain:
n_el_x = 6;
n_el_y = 30;

% build a square mesh:
mesh = geometry.build_square_mesh(l_x, l_y, n_el_x, n_el_y, 0*pi/6);


% define useful partition of the boundary nodes:
all_nodes_id = 1:length(mesh.boundary_nodes);
corner_nodes_id = find(ismember(mesh.boundary_nodes, mesh.corner_nodes));
face_nodes_id = setdiff(all_nodes_id, corner_nodes_id)';
upper_face_nodes_id = setdiff(find(mesh.nodes(mesh.boundary_nodes, 2) >= 4.0), corner_nodes_id);
lower_face_nodes_id = setdiff(find(mesh.nodes(mesh.boundary_nodes, 2) <= -4.0), corner_nodes_id);
other_faces_nodes_id = setdiff(face_nodes_id, [upper_face_nodes_id; lower_face_nodes_id]);


% define dirichlet boundary conditions:
boundary_conditions = struct();
boundary_conditions.u_n = {lower_face_nodes_id, -ones(length(lower_face_nodes_id), 1)};
boundary_conditions.u_t = {lower_face_nodes_id, zeros(length(lower_face_nodes_id), 1)};
boundary_conditions.u_x = {[corner_nodes_id; other_faces_nodes_id], zeros(length([corner_nodes_id; other_faces_nodes_id]), 1)};
boundary_conditions.u_y = {[corner_nodes_id; other_faces_nodes_id], zeros(length([corner_nodes_id; other_faces_nodes_id]), 1)};


% pde parameters:
mu = 1.0;
force_f = @(x) 0.0 * [x(:, 2) > 0.0, zeros(size(x, 1), 1)];


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
