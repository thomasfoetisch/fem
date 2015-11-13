
% prologue
clear all;
close all;
more off;
graphics_toolkit('fltk');
debug_on_interrupt(1);
debug_on_error(1);
pkg load odepkg


% geometric dimensions of the domain:
l_x = 1;
l_y = 1;

% subdivisions of the domain:
n_el_x = 128;
n_el_y = 128;

% build a square mesh:
mesh = geometry.build_square_mesh(l_x, l_y, n_el_x, n_el_y, 0*pi/6);


% define useful partition of the boundary nodes:
all_nodes_id = 1:length(mesh.boundary_nodes);
corner_nodes_id = find(ismember(mesh.boundary_nodes, mesh.corner_nodes));
upper_corner_node_id = corner_nodes_id(find(mesh.nodes(mesh.boundary_nodes(corner_nodes_id), 2) > 0.0));
lower_corner_node_id = corner_nodes_id(find(mesh.nodes(mesh.boundary_nodes(corner_nodes_id), 2) < 0.0));

face_nodes_id = setdiff(all_nodes_id, corner_nodes_id)';
upper_face_nodes_id = setdiff(find(mesh.nodes(mesh.boundary_nodes, 2) >= 0.5), corner_nodes_id);
lower_face_nodes_id = setdiff(find(mesh.nodes(mesh.boundary_nodes, 2) <= -0.5), corner_nodes_id);
other_faces_nodes_id = setdiff(face_nodes_id, [upper_face_nodes_id; lower_face_nodes_id]);


% define dirichlet boundary conditions:
poiseuille = 0;
if poiseuille
  mu = 1.0;
  boundary_conditions = struct();
  boundary_conditions.u_n = {lower_face_nodes_id, -ones(length(lower_face_nodes_id), 1)};
  boundary_conditions.u_t = {lower_face_nodes_id, zeros(length(lower_face_nodes_id), 1)};
  boundary_conditions.u_x = {[corner_nodes_id; other_faces_nodes_id], zeros(length([corner_nodes_id; other_faces_nodes_id]), 1)};
  boundary_conditions.u_y = {[corner_nodes_id; other_faces_nodes_id], zeros(length([corner_nodes_id; other_faces_nodes_id]), 1)};
  force_f = @(x) 0.0 * [x(:, 2) > 0.0, zeros(size(x, 1), 1)];
end

driven_cavity = 0;
if driven_cavity
  mu = 1.0;
  boundary_conditions = struct();
  boundary_conditions.u_x = {[upper_corner_node_id; lower_corner_node_id], [1, 1, 0, 0]};
  boundary_conditions.u_y = {corner_nodes_id, [0, 0, 0, 0]};
  boundary_conditions.u_n = {face_nodes_id, @(x) zeros(size(x, 1), 1)};
  boundary_conditions.u_t = {[upper_face_nodes_id; ...
                              other_faces_nodes_id],
			     [ones(length(upper_face_nodes_id), 1); ...
                              zeros(length(other_faces_nodes_id), 1)]};
  force_f = @(x) 0.0 * [x(:, 2) > 0.0, zeros(size(x, 1), 1)];
end

kim_and_moin = 1;
if kim_and_moin
  mu = 1.0;
  u_x_exact = @(x) - cos(pi * x(:, 1)) .* sin(pi * x(:, 2));
  u_y_exact = @(x)   sin(pi * x(:, 1)) .* cos(pi * x(:, 2));
  p_exact = @(x) zeros(size(x, 1), 1);

  boundary_conditions = struct();
  boundary_conditions.u_x = {all_nodes_id, u_x_exact};
  boundary_conditions.u_y = {all_nodes_id, u_y_exact};
  boundary_conditions.u_n = {};
  boundary_conditions.u_t = {};

  force_f = @(x) 2 * mu * pi^2 * [u_x_exact(x), u_y_exact(x)];
end


% solve the problem
[u_x, u_y, p] = stokes2d.solve(mesh, mu, force_f, boundary_conditions);

% compute the errors
velocity_x_error = sqrt(functional.error_l2_exact_p1bulle(mesh, u_x_exact, u_x));
velocity_y_error = sqrt(functional.error_l2_exact_p1bulle(mesh, u_y_exact, u_y));
velocity_error = sqrt(velocity_x_error^2 + velocity_y_error^2);
pressure_error = sqrt(functional.error_l2_exact_p1(mesh, p_exact, p));

% display errors
printf('Error on velocities: %f\n', velocity_error);
printf('Error on pressure:   %f\n', pressure_error);

% remove the bubbles:
u_x = u_x(1:size(mesh.nodes, 1));
u_y = u_y(1:size(mesh.nodes, 1));


% visualisation of the result:
vis = 0;
if vis
  f_dof = force_f(mesh.nodes) / mu / pi^2 / 2;
   
  figure(1); hold on;
  quiver(mesh.nodes(:, 1), mesh.nodes(:, 2), u_x/10, u_y/10, 0, 'r', 'Autoscale', 'off', 'LineWidth', 2);
  quiver(mesh.nodes(:, 1), mesh.nodes(:, 2), f_dof(:, 1)/10, f_dof(:, 2)/10, 0, 'g', 'Autoscale', 'off');
  gplot(mesh.adj, mesh.nodes);
  
  figure(2);
  trisurf(mesh.elements, mesh.nodes(:, 1), mesh.nodes(:, 2), p);


  % plot the isolevels of the velocity magnitude:
  u_mag = sqrt(u_x.^2 + u_y.^2)(1:size(mesh.nodes, 1));
  
  figure(4);
  contour(permute(reshape(u_mag, n_el_x + 1, n_el_y + 1), [2 1]));


  % plot some particle traces and domain boundary:
  if 0
    n_emitter = 10;
    emitters = [zeros(n_emitter, 1), linspace(-l_y / 2, l_y / 2, n_emitter + 2)(2:end-1)' ];
    traces = postprocessing.particle_trace(mesh, [u_x, u_y], emitters);

    x = [mesh.nodes(mesh.boundary_edges(:, 1), 1), ...
	 mesh.nodes(mesh.boundary_edges(:, 2), 1), ...
	 NaN(size(mesh.boundary_edges, 1), 1)]'(:);
    y = [mesh.nodes(mesh.boundary_edges(:, 1), 2), ...
	 mesh.nodes(mesh.boundary_edges(:, 2), 2), ...
	 NaN(size(mesh.boundary_edges, 1), 1)]'(:);

    figure(5); cla; hold on;
    plot(x, y, '-k');
    plot(traces(:, 1), traces(:, 2), '-r');
  end
end
