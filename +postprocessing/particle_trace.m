function coords = particle_trace(mesh, field, emitters)
  edge_adj = geometry.build_elements_edge_adjacency(mesh);
  bc = geometry.build_barycentric_coordinates(mesh);
  
  vopts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4, 'InitialStep', 0.1, 'MaxStep', 0.1);

  coords = zeros(0, 2);
  for n = 1:size(emitters, 1)
    [time, trace] = ode23(@ode_adapter,
			  [0, 200],
			  emitters(n, :), vopts, mesh, edge_adj, bc, field);
    coords = [coords; [NaN, NaN]; trace];
  end
