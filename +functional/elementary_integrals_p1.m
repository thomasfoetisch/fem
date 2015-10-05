function ints = elementary_integrals()

  % elementary integrals for p1 lagrange
  ints = struct();
  ints.one = 1/2;
  ints.phi = [1/6, 1/6, 1/6];
  ints.dphi = [-1/2, 1/2, 0,
	       -1/2, 0, 1/2];
  ints.dphidphi = zeros(2,3,2,3);
  for i = 1:2
    for j = 1:3
	ints.dphidphi(:, :, i, j) = 2 * ints.dphi * ints.dphi(i, j);
    end
  end

  ints.phidphi = zeros(3, 2, 3);
  for i = 1:3
    for j = 1:2
      for k = 1:3
	  ints.phidphi(i, j, k) = ints.phi(i) * 2 * ints.dphi(j, k);
      end
    end
  end
  
  ints.phiphi = [1/12, 1/24, 1/24;
		 1/24, 1/12, 1/24;
		 1/24, 1/24, 1/12];

  ints.phiphidphi = zeros(3, 3, 2, 3);
  for i = 1:2
    for j = 1:3
	ints.phiphidphi(:, :, i, j) = ints.phiphi * 2 * ints.dphi(i, j);
    end
  end

  ints.phiphidphidphi = zeros(3, 3, 2, 3, 2, 3);
  for i = 1:2
    for j = 1:3
	ints.phiphidphidphi(:, :, :, :, i, j) = ints.phiphidphi * 2 * ints.dphi(i, j);
    end
  end
end
