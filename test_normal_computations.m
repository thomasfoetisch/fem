

mesh = geometry.build_square_mesh(1, 1, 5, 5, 0);

x = [mesh.nodes(mesh.boundary_edges(:, 1), 1), mesh.nodes(mesh.boundary_edges(:, 2), 1), NaN(size(mesh.boundary_edges, 1), 1)]'(:);
y = [mesh.nodes(mesh.boundary_edges(:, 1), 2), mesh.nodes(mesh.boundary_edges(:, 2), 2), NaN(size(mesh.boundary_edges, 1), 1)]'(:);


plot(x, y, '-o'); hold on;
quiver(mesh.boundary_barycenters(:, 1), mesh.boundary_barycenters(:, 2), mesh.tangents(:, 1), mesh.tangents(:, 2), 0.25, 'g');
quiver(mesh.boundary_barycenters(:, 1), mesh.boundary_barycenters(:, 2), mesh.normals(:, 1), mesh.normals(:, 2), 0.25, 'r');

