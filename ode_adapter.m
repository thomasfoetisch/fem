function out = ode_adapter(t, x, mesh, edge_adj, baryc, field_p1)
  [values, ~, ~] = geometry.interpolate_on_mesh(mesh, edge_adj, baryc, field_p1, x');
  out = values;
