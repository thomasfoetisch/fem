% parameters:
l_x = 1;

n_el_x = 10;


% build the nodes:
nodes = linspace(-l_x/2, l_x/2, n_el_x + 1);


% build the elements:
elements = int16(zeros(n_el_x, 2));

for el = 1:size(elements, 1)
    elements(el, 1) = el;
    elements(el, 2) = el + 1;
end


% add randomness
perm = randperm(size(nodes, 2));
%perm = 1:size(nodes, 2);
[~, inv_perm] = sort(perm);

nodes = nodes(perm);
elements = reshape(inv_perm(reshape(elements, 1, [])), size(elements));

for e = 1:size(elements, 1)
  %elements(e, :) = circshift(elements(e, :), [0, randi(3) - 1]);
  elements(e, :) = elements(e, randperm(2));
end
elements = int16(elements);


% build the adjacency matrix:
adj = int8(zeros(size(nodes, 1), size(nodes, 1)));
for el = 1:size(elements, 1)
  adj(elements(el, 1), elements(el, 2)) = 1;
  adj(elements(el, 2), elements(el, 1)) = 1;
end
for n = 1:size(nodes, 1)
    adj(n, n) = 1;
end


% build the jacobian matrix:
jmt = zeros(1, size(elements, 1));
jac = zeros(1, size(elements, 1));

for el = 1:size(elements, 1)
  m = nodes(elements(el, 2)) - nodes(elements(el, 1));
  jmt(el) = 1/m;
  jac(el) = abs(m);
end


% elementary integrals:
one = 1;
phi = [1/2, 1/2];
dphi = [-1, 1];
phiphi = [1/3, 1/6; 1/6, 1/3];

% assemble the matrix:
mat = zeros(size(nodes, 2));

for el = 1:size(elements, 1)
  for i = 1:2
    for j = 1:2
      contrib = jac(el) * jmt(el) * dphi(i) * jmt(el) * dphi(j);
      mat(elements(el, i), elements(el, j)) = mat(elements(el, i), elements(el, j)) + contrib;
    end
  end
end


% assemble the rhs:
rhs = zeros(size(nodes, 2), 1);

for el = 1:size(elements, 1)
  for j = 1:2
      contrib = jac(el) * (2) * phi(j);
      rhs(elements(el, j)) = rhs(elements(el, j)) + contrib;
  end
end



% build the boundary nodes list:
boundary_nodes = inv_perm([1, n_el_x + 1]);


% apply the dirichlet boundary conditions:
for bcn = 1:size(boundary_nodes, 2)
  mat(boundary_nodes(bcn), boundary_nodes(bcn)) = 1e10;
  rhs(boundary_nodes(bcn)) = 1e10 * 0;
end


u = linsolve(mat, rhs);

u_exact_f = @(x) -(x + l_x/2) .* (x - l_x/2);
u_exact = u_exact_f(nodes);

figure(1); cla; hold on
plot(nodes(inv_perm), u(inv_perm), 'or');
plot(nodes(inv_perm), u_exact(inv_perm), '-g')

