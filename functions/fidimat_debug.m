% function outM = fidimat(arg1,arg2,arg3,arg4)
%   % FIDIMAT Generate Finite Difference Spatial Sparcity Matrices
%   %   FIDIMAT(l,m,ord,bctype) A function to generate the biharmonic coefficients for
%   %   2D given the number of ALL grid points given by l and m.
%   %   bctype denotes boundary condition
%   %
%   %         Returns a Sparse matrix BH
%   %
%   %         bctype  % boundary condition type: 1: simply supported, 2: clamped
%   %         l       % number of total grid points Y axis
%   %         m       % number of total grid points X axis
%   %         ord     % order of the matrix (string)
%   %
%   %                 % Valid order inputs
%   %                 % 'x-', 'x+', 'xx', 'xxxx', 'laplace', 'biharm'
%
%   % Arguement check
%   if nargin<2
%     error('Not enough input arguements')
%
%   elseif nargin==2
%     if ischar(arg2)
%       ord = arg2;
%       l = arg1;
%       m = 1;
%       bctype = 1;
%     else
%       bctype = arg2;
%       ord = 'xx';
%       l = arg1;
%       m = 1;
%     end
%
%   elseif nargin==3
%
%     if ischar(arg3)
%       ord = arg2;
%       l = arg1;
%       m = arg2;
%       bctype = 1;
%     else
%       bctype = arg3;
%       ord = 'xx';
%       l = arg1;
%       m = arg2;
%     end
%
%   elseif nargin==4
%     l = arg1;
%     m = arg2;
%     ord = arg3;
%     bctype = arg4;
%
%   elseif nargin>1
%     error('Too many input arguements')
%
%   end
clear all
close all

l = 4;
m = 4;
ord = 'xx'
bctype = 1
Iy = speye(l); % identity matrix for y axis
Ix = speye(m); % identity matrix for x axis
ss = l*m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 1D Case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% if only one dimension is stipulated the function will return only
%% a 1D FD matrix. In this case input is by the letter x

%% NOTE May add in that y can be used but, for now this is left
% limited to x to hopefully avoid potential clash or ambiguity.

if m == 1

  switch ord

  case 'x-'
    outM = spdiags([-ones(l,1),ones(l,1)],-1:0,speye(l));

  case 'x+'
    outM = spdiags([ones(l,1),-ones(l,1)],0:1,speye(l));

  case 'x.'
    outM = spdiags([-ones(l,1),zeros(l,1),ones(l,1)],[-1:1],speye(l));

  case 'xx'
    XX = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
    XX(1,2) = (bctype-1)*2; XX(l,l-1) = (bctype-1)*2;
    outM = XX;

  case 'xxxx'
    XX = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
    XX(1,2) = (bctype-1)*2; XX(l,l-1) = (bctype-1)*2;
    outM = XX;
    outM = XX^2;

  otherwise
    error('something went wrong, check your arguements');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% 2D Case
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % if two dimension arguements are given then the matri will out put the
  % corresponding 2D Matrix.
  %
  % y is for columns of a matrix whilst x is for rows.
  % x matrices in this case will contain 'off centre' diagonals corresponding to
  % the number of elements in the row.
  %
  % If that didn't make sense don't worry, its my fault entirely.

  % set the correct points to 0
  % BH(:,[1:l+1, ss-l:ss]) = 0;BH([1:l+1, ss-l:ss],:) = 0;
  % BH(:,[2*l:l:ss, (2*l)+1:l:ss]) = 0; BH([2*l:l:ss, (2*l)+1:l:ss],:) = 0;

else
  switch ord

  case 'y-'
    Y = spdiags([-ones(l,1),ones(l,1)],-1:0,speye(l));
    outM = kron(Ix,Y);

  case 'y+'
    Y = spdiags([ones(l,1),-ones(l,1)],0:1,speye(l));
    outM = kron(Ix,Y);

  case 'y.'
    Y = spdiags([-ones(l,1),zeros(l,1),ones(l,1)],[-1:1],speye(l));
    outM = kron(Ix,Y);

  case 'yy'
    YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
    YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
    YY(1,2) = (bctype-1)*2; YY(l,l-1) = (bctype-1)*2;
    outM = kron(Ix,Y);

  case 'yyyy'
    YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
    YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
    YY(1,2) = (bctype-1)*2; YY(l,l-1) = (bctype-1)*2;
    outM = kron(Ix,YY)^2;

  case 'x-'
    X = spdiags([-ones(m,1),ones(m,1)],-1:0,speye(m));
    outM = kron(X,Iy);

  case 'x+'
    X = spdiags([ones(m,1),-ones(m,1)],0:1,speye(ss));
    outM = kron(X,Iy);

  case 'x.'
    X = spdiags([-ones(m,1),zeros(m,1),ones(m,1)],[-1:1],speye(m));
    outM = kron(X,Iy);

  case 'xx'
    XX = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
    XX = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
    XX(1,2) = (bctype-1)*2; XX(l,l-1) = (bctype-1)*2;
    outM = kron(XX,Iy);

  case 'xxxx'
    XX = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
    XX = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
    XX(1,2) = (bctype-1)*2; XX(l,l-1) = (bctype-1)*2;
    outM = kron(XX,Iy)^2;

  case 'xxyy'
    YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
    YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
    YY(1,2) = (bctype-1)*2; YY(l,l-1) = (bctype-1)*2;

    XX = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
    XX = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,speye(l));
    XX(1,2) = (bctype-1)*2; XX(l,l-1) = (bctype-1)*2;

    outM = kron(Ix,YY)*kron(XX,Iy);

  case 'laplace'

    XX = spdiags([ones(m,1),-2*ones(m,1),ones(m,1)],-1:1,Ix);
    XX(1,2) = (bctype-1)*2; XX(m,m-1) = (bctype-1)*2;
    YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,Iy);
    YY(1,2) = (bctype-1)*2; YY(l,l-1) = (bctype-1)*2;
    LA = kron(XX,Iy) + kron(Ix,YY);
    outM = LA;

  case 'biharm'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Building Bi-Harmonic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    XX = spdiags([ones(m,1),-2*ones(m,1),ones(m,1)],-1:1,Ix);
    XX(1,2) = (bctype-1)*2; XX(m,m-1) = (bctype-1)*2;
    YY = spdiags([ones(l,1),-2*ones(l,1),ones(l,1)],-1:1,Iy);
    YY(1,2) = (bctype-1)*2; YY(l,l-1) = (bctype-1)*2;
    LA = kron(XX,Iy) + kron(Ix,YY);
    BH = LA*LA;

    % set the correct points to 0
    BH(:,[1:l+1, ss-l:ss]) = 0;BH([1:l+1, ss-l:ss],:) = 0;
    BH(:,[2*l:l:ss, (2*l)+1:l:ss]) = 0; BH([2*l:l:ss, (2*l)+1:l:ss],:) = 0;

    outM = BH;

  otherwise
    error('something went wrong, check your arguements');

  end
end
