function g = build_clement_interpolant(mesh, field_p0)
  g = zeros(size(mesh.nodes, 1), 1);
  for node = 1:size(mesh.nodes, 1)
    elements = mesh.node_el_adj{node};
    g(node) = sum(mesh.jac(elements) .* field_p0(elements)) / sum(mesh.jac(elements));
  end
end
