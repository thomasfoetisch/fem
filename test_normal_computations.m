

mesh = geometry.build_square_mesh(1, 1, 5, 5, 0);

x = [mesh.nodes(mesh.boundary_edges(:, 1), 1), mesh.nodes(mesh.boundary_edges(:, 2), 1), NaN(size(mesh.boundary_edges, 1), 1)]'(:);
y = [mesh.nodes(mesh.boundary_edges(:, 1), 2), mesh.nodes(mesh.boundary_edges(:, 2), 2), NaN(size(mesh.boundary_edges, 1), 1)]'(:);


plot(x, y, '-o'); hold on;
quiver(mesh.boundary_barycenters(:, 1), mesh.boundary_barycenters(:, 2), mesh.tangents(:, 1), mesh.tangents(:, 2), 0.25, 'g');
quiver(mesh.boundary_barycenters(:, 1), mesh.boundary_barycenters(:, 2), mesh.normals(:, 1), mesh.normals(:, 2), 0.25, 'r');



normal_nodes = zeros(size(mesh.boundary_nodes, 1), 2);
for node = 1:size(mesh.boundary_nodes, 1)
  normal_nodes(node, :) = (mesh.normals(mesh.node_bel_adj(node, 1), :) + mesh.normals(mesh.node_bel_adj(node, 2), :)) / 2;
end


quiver(mesh.nodes(mesh.boundary_nodes, 1), mesh.nodes(mesh.boundary_nodes, 2), normal_nodes(:, 1), normal_nodes(:, 2), 'c')
