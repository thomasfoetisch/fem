

mesh = geometry.build_square_mesh(1, 1, 2, 2, pi/6);

x = [mesh.nodes(mesh.boundary_edges(:, 1), 1), mesh.nodes(mesh.boundary_edges(:, 2), 1), NaN(size(mesh.boundary_edges, 1), 1)]'(:);
y = [mesh.nodes(mesh.boundary_edges(:, 1), 2), mesh.nodes(mesh.boundary_edges(:, 2), 2), NaN(size(mesh.boundary_edges, 1), 1)]'(:);


plot(x, y, '-o'); hold on;
quiver(mesh.boundary_barycenters(:, 1), mesh.boundary_barycenters(:, 2), mesh.tangents(:, 1), mesh.tangents(:, 2), 0.25, 'g');
quiver(mesh.boundary_barycenters(:, 1), mesh.boundary_barycenters(:, 2), mesh.normals(:, 1), mesh.normals(:, 2), 0.25, 'r');

quiver(mesh.nodes(mesh.boundary_nodes, 1), mesh.nodes(mesh.boundary_nodes, 2), mesh.normals_nodes(:, 1), mesh.normals_nodes(:, 2), 'c');
quiver(mesh.nodes(mesh.boundary_nodes, 1), mesh.nodes(mesh.boundary_nodes, 2), mesh.tangents_nodes(:, 1), mesh.tangents_nodes(:, 2), 'k');

plot(mesh.nodes(mesh.boundary_nodes, 1), mesh.nodes(mesh.boundary_nodes, 2), 'ok');
