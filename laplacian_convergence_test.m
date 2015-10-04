% test de convergeance du schema:
h = zeros(1, 4);
err = zeros(1, 4);

for i = 1:length(h)
    h(i) = 2^(-2 - i); 
    [nodes, ~, u, u_exact, err(i)] = laplacian(h(i));
end
order = polyfit(log(h), log(err), 1)(1);

printf('The convergence order is: %f\n', order);
