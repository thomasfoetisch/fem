function [element, coordinates] = find_element_at_point_naive(mesh, edge_adj, baryc, point)

  found = 0;
  for n = 1:size(mesh.elements, 1)
    point_coordinates = baryc(:, :, n) ...
		      * [point, 1]';

    valid_coordinates = point_coordinates >= 0 && point_coordinates <= 1;
    if valid_coordinates
       element = n;
       coordinates = point_coordinates;
       found = 1;
    else 
    end
  end

  if !found
    error(sprintf('Impossible to find a corresponding element for point (%f, %f).\n', point(1), point(2)));
  end
