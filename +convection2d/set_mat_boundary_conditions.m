function mat = set_mat_boundary_conditions(mesh, mat, dbc_f, dbc_nodes)

  [i, j, v] = find(mat);
  for bcn = 1:size(dbc_nodes, 1)
    idx = find(i == dbc_nodes(bcn));
    mat(i(idx), j(idx)) = 0;
    mat(dbc_nodes(bcn), ...
	dbc_nodes(bcn)) = 1;
  end
end
    
