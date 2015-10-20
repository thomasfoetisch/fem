clear all;
more off;
graphics_toolkit('gnuplot');


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
force_f = @(x) 0 * [x(:, 2) > 0.25, zeros(size(x, 1), 1)];


% build the mesh:
mesh = geometry.build_square_mesh(l_x, l_y, n_el_x, n_el_y, 0);


% define useful partition of the boundary nodes:
all_nodes_id = 1:length(mesh.boundary_nodes);
corner_nodes_id = find(ismember(mesh.boundary_nodes, mesh.corner_nodes));
face_nodes_id = setdiff(all_nodes_id, corner_nodes_id)';
upper_face_nodes_id = setdiff(find(mesh.nodes(mesh.boundary_nodes, 2) >= 0.5), corner_nodes_id);
other_faces_nodes_id = setdiff(face_nodes_id, upper_face_nodes_id);


% define dirichlet boundary conditions:
boundary_conditions = struct();
boundary_conditions.u_x = {corner_nodes_id, [0, 0, 0, 0]};
boundary_conditions.u_y = {corner_nodes_id, [0, 0, 0, 0]};
boundary_conditions.u_n = {face_nodes_id, @(x) zeros(size(x, 1), 1)};
boundary_conditions.u_t = {[upper_face_nodes_id; ...
                            other_faces_nodes_id],
			   [0.1 * ones(length(upper_face_nodes_id), 1); ...
                            zeros(length(other_faces_nodes_id), 1)]};


% solve the system:
[u_x, u_y, p, newton_errors] = navierstokes2d.solve(mesh, force_f, ctx, boundary_conditions);


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
