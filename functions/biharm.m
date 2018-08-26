function BH = biharm(l,m,bctype)
  % BIHARM Generate Biharmonic Matrix
  %   BIHARM(l,m,bctype) A function to generate the biharmonic coefficients for
  %   2D given the number of internal grid points given by l and m.
  %   bctype denotes boundary condition
  %
  %         Returns a Sparse matrix BH
  %
  %         bctype  % boundary condition type: 1: simply supported, 2: clamped
  %         l       % number of internal grid points Y axis
  %         m       % number of internal grid points X axis
  %

  l = l + (bctype-1)*2;
  m = m + (bctype-1)*2;
  Iy = speye(l); % identity matrix for y axis
  Ix = speye(m); % identity matrix for x axis

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Building Bi-Harmonic
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  XX = spdiags([ones(m,1),-2*ones(m,1),ones(m,1)],-1:1,Ix);
  XX(1,2) = bctype; XX(m,m-1) = bctype;
  YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,Iy);
  YY(1,2) = bctype; YY(l,l-1) = bctype;
  LA = kron(XX,Iy) + kron(Ix,YY);
  BH = LA*LA;
  % BH for clamped conditions will be larger which will cause unnecessary computaion
  % costs
