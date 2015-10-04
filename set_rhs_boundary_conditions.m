function rhs = set_rhs_boundary_conditions(mesh, rhs, dbc_f, dbc_nodes)

  for bcn = 1:size(dbc_nodes, 1)
    rhs(dbc_nodes(bcn)) = dbc_f(mesh.nodes(dbc_nodes(bcn), :));
  end
end
    
