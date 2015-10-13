function [phiphidphi, dphidphi] = test_numerical_integration()

basis    = {@(x, y) (1 - x - y),  @(x, y) x, @(x, y) y,  @(x, y) 27 * (1 - x - y) .* x .* y};
d_basis =  {@(x, y) -1,           @(x, y) 1, @(x, y) 0,  @(x, y) 27 * y .* (1 - 2*x - y),
	    @(x, y) -1,           @(x, y) 0, @(x, y) 1,  @(x, y) 27 * x .* (1 - x - 2*y)};

tol = 1.e-9;

phiphidphi = zeros(4, 4, 2, 4);
dphidphi = zeros(2, 4, 2, 4);

for i = 1:4
  for j = 1:4
    for k = 1:2
      for l = 1:4
	printf('%i %i %i %i\n', i, j, k, l);
	phiphidphi(i, j, k, l) = quadgk(@(x) f_phiphidphi(x, i, j, k, l, basis, d_basis), 0, 1, tol);
      end
    end
  end
end

for i = 1:2
  for j = 1:4
    for k = 1:2
      for l = 1:4
	printf('%i %i %i %i\n', i, j, k, l);
	dphidphi(i, j, k, l) = quadgk(@(x) f_dphidphi(x, i, j, k, l, basis, d_basis), 0, 1, tol);
      end
    end
  end
end



function q = f_phiphidphi(x, i, j, k, l, basis, d_basis)
  tol = 1.e-9;
  q = ones(size(x));
  for p = 1:length(q)
    f = @(y) basis{i}(x(p), y) .* basis{j}(x(p), y) .* d_basis{k, l}(x(p), y);
    q(p) = quadgk(f, 0, 1 - x(p), tol);
  endfor
end


function q = f_dphidphi(x, i, j, k, l, basis, d_basis)
  tol = 1.e-9;
  q = ones(size(x));
  for p = 1:length(q)
    f = @(y) d_basis{i, j}(x(p), y) .* d_basis{k, l}(x(p), y);
    q(p) = quadgk(f, 0, 1 - x(p), tol);
  endfor
end

end
