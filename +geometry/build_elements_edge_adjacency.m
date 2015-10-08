function adj = build_elements_edge_adjacency(mesh)
  edges = [mesh.elements(:, [1 2]); ...
           mesh.elements(:, [1 3]); ...
           mesh.elements(:, [2 3])];
  edges = sort(edges, 2);
  edges = [edges, repmat((1:size(mesh.elements, 1))', [3 1])];
  edges = sortrows(edges);

  adj = zeros(size(mesh.elements, 1), 3);
  edge = 1;
  while (edge < size(edges, 1))
    if edges(edge, 1:2) == edges(edge + 1, 1:2)
       edge_id_1 = find(mesh.elements(edges(edge,     3), :) == setdiff(mesh.elements(edges(edge,     3), :), edges(edge,     1:2)));
       edge_id_2 = find(mesh.elements(edges(edge + 1, 3), :) == setdiff(mesh.elements(edges(edge + 1, 3), :), edges(edge + 1, 1:2)));
       adj(edges(edge, 3), edge_id_1) = edges(edge + 1, 3);
       adj(edges(edge + 1, 3), edge_id_1) = edges(edge, 3);
       edge = edge + 2;
    else
      edge = edge + 1;
    end
  end
