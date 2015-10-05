function nodes = input_flow_surface_nodes(mesh, w_f)
  
  % build the boundary node list:
  %  edges = [n1, n2, edge_count, edge_id, element_id]
  edges = int16(zeros(3*size(mesh.elements, 1), 5));
  edges_number = 0;
  for el = 1:size(mesh.elements)
    for k = 1:3
      kk = mod(k, 3) + 1;
      edge = sort(mesh.elements(el, [k, kk]));
      loc = find(edges(:, 1) == edge(1) & edges(:, 2) == edge(2));
      if isempty(loc)
	edges(edges_number + 1, :) = [edge, 1, mesh.elements(el, mod(kk, 3) + 1), el];
	edges_number = edges_number + 1;
      else
	edges(loc, 3) = edges(loc, 3) + 1;
      end
    end
  end


  % keep only the edges where the count is one:
  boundary_edges = edges(find(edges(:, 3) == 1), :);
  

  % compute the normals to the selected edges:
  normals = zeros(size(boundary_edges, 1), 2);
  rot = [0, -1; 1,0];
  for ed = 1:size(boundary_edges, 1)
    normals(ed, :) = (mesh.nodes(boundary_edges(ed, 2), :) ...
		      - mesh.nodes(boundary_edges(ed, 1), :)) * rot;
    
    second_edge = mesh.nodes(boundary_edges(ed, 4), :) ...
		      - mesh.nodes(boundary_edges(ed, 1), :);
    if normals(ed, :) * second_edge' > 0
      normals(ed, :) = - normals(ed, :) / sqrt(sum(normals(ed, :).^2));
    else
      normals(ed, :) = normals(ed, :) / sqrt(sum(normals(ed, :).^2));
    end
  end

  % compute the coordinates of the mid points of the selected boundary edges:
  mid_edges = (mesh.nodes(boundary_edges(:, 1), :) 
	       + mesh.nodes(boundary_edges(:, 2), :)) / 2;


  % compute the R^2 scalar product between the normal and w_f:
  input_edges = logical(zeros(size(boundary_edges, 1), 1));
  for ed = 1:size(boundary_edges, 1)
    input_edges(ed) = w_f(mid_edges(ed, :)) ...
		      * normals(ed, :)' < 0;
  end

  nodes = unique(sort(boundary_edges(input_edges, [1 2])(:)));
end
