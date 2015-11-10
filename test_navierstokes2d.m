% prologue
clear all;
close all;
more off;
graphics_toolkit('gnuplot');
debug_on_interrupt(1);


% geometrical parameters of the domain:
l_x = 1;
l_y = 1;


% subdivisions of the domain:
n_el_x = 80;
n_el_y = 80;


% pde parameters:
ctx = struct();
ctx.rho = 1;
ctx.laminar_viscosity = 0.001;
ctx.smagorinsky_coefficient = 0.001;
ctx.smagorinsky_caracteristic_length = max(l_x, l_y);
force_f = @(x) 0 * [x(:, 2) > 0.25, zeros(size(x, 1), 1)];


% build the mesh:
mesh = geometry.build_square_mesh(l_x, l_y, n_el_x, n_el_y, 0);


% define useful partition of the boundary nodes:
all_nodes_id = 1:length(mesh.boundary_nodes);
corner_nodes_id = find(ismember(mesh.boundary_nodes, mesh.corner_nodes));
face_nodes_id = setdiff(all_nodes_id, corner_nodes_id)';
upper_face_nodes_id = setdiff(find(mesh.nodes(mesh.boundary_nodes, 2) >= 0.5), corner_nodes_id);
lower_face_nodes_id = setdiff(find(mesh.nodes(mesh.boundary_nodes, 2) <= -0.5), corner_nodes_id);
other_faces_nodes_id = setdiff(face_nodes_id, upper_face_nodes_id);


% define dirichlet boundary conditions:
boundary_conditions = struct();
boundary_conditions.u_x = {corner_nodes_id, [0, 0, 0, 0]};
boundary_conditions.u_y = {corner_nodes_id, [0, 0, 0, 0]};
boundary_conditions.u_n = {face_nodes_id, @(x) zeros(size(x, 1), 1)};
boundary_conditions.u_t = {[upper_face_nodes_id; ...
                            other_faces_nodes_id],
			   [1 * ones(length(upper_face_nodes_id), 1); ...
                            zeros(length(other_faces_nodes_id), 1)]};


% solve the system:
[u_x, u_y, p, newton_errors] = navierstokes2d.solve(mesh, force_f, ctx, boundary_conditions);


% visualize the result:
vis = 1;
if vis
  % plot the mesh, velocity field, force field:
  force_nodes = force_f(mesh.nodes);

  figure(1); hold on; cla;
  gplot(mesh.adj, mesh.nodes);
  h = quiver(mesh.nodes(:, 1), mesh.nodes(:, 2), ...
	     u_x(1:size(mesh.nodes, 1)), ...
	     u_y(1:size(mesh.nodes))); set(h, 'color', 'r');
  h = quiver(mesh.nodes(:, 1), mesh.nodes(:, 2), ...
	     force_nodes(:, 1), ...
	     force_nodes(:, 2)); set(h, 'color', 'g');

  % plot the pressure field:
  figure(2); cla;
  trisurf(mesh.elements, mesh.nodes(:, 1), mesh.nodes(:, 2), p);


  % plot the isolevels of the velocity magnitude:
  u_mag = sqrt(u_x.^2 + u_y.^2)(1:size(mesh.nodes, 1));
  
  figure(3); cla;
  contour(permute(reshape(u_mag, n_el_x + 1, n_el_y + 1), [2 1]));


  % plot some particle traces and domain boundary:
  n_emitter = 10;
  emitters = [zeros(n_emitter, 1), linspace(-l_y / 2, l_y / 2, n_emitter + 2)(2:end-1)' ];
  traces = postprocessing.particle_trace(mesh, [u_x, u_y], emitters);

  x = [mesh.nodes(mesh.boundary_edges(:, 1), 1), ...
       mesh.nodes(mesh.boundary_edges(:, 2), 1), ...
       NaN(size(mesh.boundary_edges, 1), 1)]'(:);
  y = [mesh.nodes(mesh.boundary_edges(:, 1), 2), ...
       mesh.nodes(mesh.boundary_edges(:, 2), 2), ...
       NaN(size(mesh.boundary_edges, 1), 1)]'(:);

  figure(4); cla; hold on;
  plot(x, y, '-k');
  plot(traces(:, 1), traces(:, 2), '-r');
end
