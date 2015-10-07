function ints = elementary_integrals_p1bulle()

  % elementary integrals for p1 lagrange
  ints = struct();
  ints.one = 1/2;
  ints.phi = [1/6, 1/6, 1/6, 0.225];
  ints.dphi = [-1/2, 1/2,   0, 0,
	       -1/2,   0, 1/2, 0];
ints.dphidphi = reshape([ 0.5  ,  0.5  , -0.5  ,  0.   ,  0.   , -0.5  ,  0.   ,  0.   , ...
                           0.5  ,  0.5  , -0.5  ,  0.   ,  0.   , -0.5  ,  0.   ,  0.   , ...
                           -0.5  , -0.5  ,  0.5  ,  0.   ,  0.   ,  0.5  ,  0.   ,  0.   , ...
                           0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   , ...
                           0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   , ...
                           -0.5  , -0.5  ,  0.5  ,  0.   ,  0.   ,  0.5  ,  0.   ,  0.   , ...
                           0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  4.05 ,  2.025, ...
                           0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  2.025,  4.05 ], [2, 4, 2, 4]);


ints.phidphi = reshape([-0.16666667, -0.16666667, -0.16666667, -0.225     , -0.16666667, ...
                          -0.16666667, -0.16666667, -0.225     ,  0.16666667,  0.16666667, ...
                          0.16666667,  0.225     ,  0.        ,  0.        ,  0.        , ...
                          0.        ,  0.        ,  0.        ,  0.        ,  0.        , ...
                          0.16666667,  0.16666667,  0.16666667,  0.225     ,  0.225     , ...
                          -0.225     ,  0.        ,  0.        ,  0.225     ,  0.        , ...
                          -0.225     ,  0.        ], [4, 2, 4]);
  
  
ints.phiphi = reshape([ 0.08333333,  0.04166667,  0.04166667,  0.075     ,  0.04166667, ...
                         0.08333333,  0.04166667,  0.075     ,  0.04166667,  0.04166667, ...
                         0.08333333,  0.075     ,  0.075     ,  0.075     ,  0.075     , ...
                                0.14464286], [4, 4]);