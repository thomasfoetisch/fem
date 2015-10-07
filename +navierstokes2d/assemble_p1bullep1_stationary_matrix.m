function mat = assemble_p1bullep1_stationary_matrix(mesh, dof_map, ints, laminar_viscosity, smagorinsky_coefficient, smagorinsky_caracteristic_length)

  % problem sizes:
  n_nodes = size(mesh.nodes, 1);
  n_elems = size(mesh.elements, 1);
  n_u_dof = n_nodes + n_elems;
  n_p_dof = n_nodes;
  n_dof = 2 * n_u_dof + n_p_dof;

  mat = spalloc(n_dof, n_dof, 0);
  
  for el = 1:n_elems
    elem_u = zeros(4, 2, 4, 2);
    

    elem_p = zeros(4, 2, 3);


    elem_div = zeros(3, 4, 2);

  end
  
