mesh = geometry.build_square_mesh(1, 1, 10, 10, 0);

edge_adj = geometry.build_elements_edge_adjacency(mesh);

for el = 1:size(mesh.elements, 1)
  for edge = 1:3
    if edge_adj(el, edge)
      edge_nodes = mesh.elements(el, setdiff(1:3, edge));

      if length(setdiff(mesh.elements(edge_adj(el, edge), :), edge_nodes)) != 1
	 error(sprintf('edge %i of element %i is wrong', edge, el));
      end
    end
  end
end
