function LA = laplace2(l,m,bctype)
  % LAPLACE2 Generate Laplacian Matrix
  %   LAPLACE2(l,m,bctype) A function to generate the biharmonic coefficients for
  %   2D given the number of ALL grid points given by l and m.
  %   bctype denotes boundary condition
  %
  %         Returns a Sparse matrix LA
  %
  %         bctype  % boundary condition type: 1: simply supported, 2: clamped
  %         l       % number of internal grid points Y axis
  %         m       % number of internal grid points X axis
  %

  l = l;
  m = m;
  Iy = speye(l); % identity matrix for y axis
  Ix = speye(m); % identity matrix for x axis
  ss = l*m;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Building Bi-Harmonic
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  XX = spdiags([ones(m,1),-2*ones(m,1),ones(m,1)],-1:1,Ix);
  XX(1,2) = (bctype-1)*2; XX(m,m-1) = (bctype-1)*2;
  YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,Iy);
  YY(1,2) = (bctype-1)*2; YY(l,l-1) = (bctype-1)*2;
  LA = kron(XX,Iy) + kron(Ix,YY);

  LA(:,[1:l+1, ss-l:ss]) = 0; LA([1:l+1, ss-l:ss],:) = 0;
  LA(:,[2*l:l:ss, (2*l)+1:l:ss]) = 0; LA([2*l:l:ss, (2*l)+1:l:ss],:) = 0;
