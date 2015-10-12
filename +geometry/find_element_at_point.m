function [element, coordinates, path_length] = find_element_at_point(mesh, edge_adj, baryc, point, hint)
  epsilon_tol = 1.e-12;
  init_hint = hint;
  path = [ hint ];
  
  found = 0;
  max_iterations = size(mesh.elements, 1);
  n_iterations = 0;
  current_element = hint;
  while !found && n_iterations < max_iterations
    point_coordinates = baryc(:, :, current_element) ...
		      * [point, 1]';

    valid_coordinates = point_coordinates >= - epsilon_tol && point_coordinates <= 1 + epsilon_tol;
    if valid_coordinates
      found = 1;
    else 
      [~, i] = min(point_coordinates);
      if (edge_adj(current_element, i))
	current_element = edge_adj(current_element, i);
	path(end + 1) = current_element;
      else
	error('The mesh is probably non convex.');
      end
    end
    n_iterations = n_iterations + 1;
  end

  if found
    element = current_element;
    coordinates = point_coordinates;
    path_length = n_iterations - 1;
  else
    figure();
    gplot(mesh.adj, mesh.nodes);hold on;
    plot(mesh.barycenters(path, 1), mesh.barycenters(path, 2), 'og-');
    plot(mesh.barycenters(init_hint, 1), mesh.barycenters(init_hint, 2), 'ob');
    plot(point(1), point(2), 'r');
    keyboard();
    error(sprintf('Impossible to find a corresponding element for point (%f, %f) with hint %i.\n', point(1), point(2), init_hint));
  end
