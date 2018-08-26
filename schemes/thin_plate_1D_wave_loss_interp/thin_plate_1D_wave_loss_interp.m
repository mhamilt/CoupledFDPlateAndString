%%%%% FDTD 2D Wave Model
%%%%% Matthew Hamilton s0674653
%%%%% Description:
%%%%%
%%%%% A Kirchhoff Thin Plate FDTD Model

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Instrument File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plate_and_string_vars % file containing all relevent parameters for this set of schemes
SR = OSR*SR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Error Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Array of variable names and conditions they must meet.

if ~skpErr
  paraLs = {'SR'};
  paraCons = [(SR<1e3),];

  if any(paraCons)

    for paraError = paraLs(paraCons)

      switch char(paraError)

        %%% WARNING ERRORS %%%
      case 'string'
        % case

        %%% FATAL ERRORS %%%

      otherwise
        error([sprintf('The following variables are invalid\n\n')...
        sprintf('''%s'' \n',paraLs{paraCons}) ...
        sprintf('\nCheck their values and try again\n')])

      end

    end

  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Error checking complete')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plate Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Motion Coefficients
pl_D = (pl_E*(pl_H)^3)/(12*(1-(nu^2)));
pl_kappa = sqrt(pl_D / (pl_rho*pl_H) );

%%%%% Scheme Spacing
k = 1/SR;                           % time step
pl_hmin = 2*sqrt(k*pl_kappa);       % NSS (equation 12.5)
pl_N = floor(pl_L/pl_hmin);         % number of segments
pl_h = pl_L/(pl_N);                 % adjusted grid spacing
pl_mu = (pl_kappa * k)/(pl_h^2);    % scheme parameter
Nf = floor(SR*Tf);                  % number of time steps

pl_N = pl_N+1;                      % for includng only internal grid points
ss = pl_N*pl_N;                     % total grid size.

%%% Loss coefficients


%%%%% I/O
Nf = floor(SR*Tf);          % number of time steps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% String Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% derived parameters

st_A = pi*st_r^2;                       % string cross-sectional area
st_I = 0.25*pi*st_r^4;                  % string moment of intertia

st_c = sqrt(st_T/(st_rho*st_A));        % wave speed
st_K = sqrt(st_E*st_I/(st_rho*st_A));   % stiffness constant

%%%%% grid

st_hmin = st_c*k;
st_N = floor(st_L/st_hmin);          % number of segments (N+1 is number of grid points)
st_h = st_L/st_N;                    % adjusted grid spacing
st_lambda = st_c*k/st_h;             % Courant number
st_mu = st_K*k/st_h^2;               % numerical stiffness constant

st_N = st_N+1;                       % change st_N to be number of grid points
%%%%% I/O
li = floor(st_xi*st_N);         % grid index of excitation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % String Loss coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if losstype ==0
  % frequency independant loss
  st_sigma0 = 0;
  st_sigma1 = 0;

  pl_sigma0 = 0;
  pl_sigma1 = 0;
end

if losstype ==1
  % frequency independant loss
  st_sigma0 = 6*log(10)/st_loss(1,2);
  % st_sigma0 = 0
  st_sigma1 = 0;

  pl_sigma0 = 6*log(10)/pl_loss(1,2);
  % pl_sigma0 = 0
  pl_sigma1 = 0;
end

if losstype == 2
  %% this is simply loss from NSS with 'c' removed as tension is not present.
  st_sigma0 = 6*log(10)/st_loss(1,2);


  pl_z1 = 2*pl_kappa*(2*pi*pl_loss(1,1))/(2*pl_kappa.^2); %from 1D
  pl_z2 = 2*pl_kappa*(2*pi*pl_loss(2,1))/(2*pl_kappa.^2);

  pl_sigma0 = 6*log(10)*(-pl_z2/pl_loss(1,2) + pl_z1/pl_loss(2,2))./(pl_z1-pl_z2);
  pl_sigma1 = 6*log(10)*(1/pl_loss(1,2) - 1/pl_loss(2,2))./(pl_z1-pl_z2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % Read In/Out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lo = rp*pl_N;
lo = floor(sub2ind([pl_N pl_N],lo(1), lo(2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IV Read In/Out Error Check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check that xo and xi are greater than h distance away from boundaries
if any([(st_xi*st_L)<st_h (st_L-(st_xi*st_L))<st_h])
  warning([sprintf('xo is too close to string boundary\n'),...
  sprintf('changing xo to equal h\n')]);
  st_xi = abs(round(st_xi)-st_h);
  if li<=2
    li = 3;
  elseif li>=st_N
    li = st_N-1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Force Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st_f = sparse(Nf,1);                        % input force signal
durint = floor(st_dur*SR);                 % duration of force signal, in samples
exc_st_int = (floor(st_exc_st*SR))+1;      % start time index for excitation
durf = exc_st_int:exc_st_int+durint-1;     % sample values of force
st_d0 = (k^2/(st_h*st_rho*st_A));          % force foefficient

%%% reevaluate relevant points in force vector
st_f(durf) = st_famp*0.5*(1-cos((2/itype)*pi.*(durf/durint)))*st_d0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plate Coefficient Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LA = laplace2(pl_N,pl_N,bctype); % Laplacian Matrix
BH = biharm2(pl_N,pl_N,bctype); % biharmonic matrix
pl_mI = fidimat(pl_N,pl_N,'I');

pl_mA = (1/(1+k*pl_sigma0))*speye(ss);
pl_mB = (-(pl_mu^2)*BH + 2*pl_mI + (2*k*pl_sigma1/pl_h^2)*LA) * pl_mA;
pl_mC = (-(1-pl_sigma0*k)*pl_mI - (2*k*pl_sigma1/pl_h^2)*LA) * pl_mA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% String Coefficients Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% scheme coefficients
Dn = fidimat(st_N,'x+',bctype);
Dnn = fidimat(st_N,'xx');
st_mI = speye(st_N);

st_mA = (1/(1+k*st_sigma0))*speye(st_N);

%%% String sparsity matrix, uncoupled.
% from definition of coupling point
% a forward difference instead of 0 or a continuation of Dnn
% w(1) = lambda^2*Dn(1,:)*w1 + 2*w1(1) - w2(1) - tau_s*cpl_v
Dnn(1,:) = Dn(1,:);
st_mB = (2*speye(st_N)+(st_lambda^2)*Dnn) * st_mA ;

st_mC = -(1-st_sigma0*k)*speye(st_N) * st_mA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Couple Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Total Number of Grid points
Nt = ss+ st_N;

%%% floored coupling points
pl_w0 = sub2ind([pl_N, pl_N],floor(pl_ctr(1,1)*pl_N),floor(pl_ctr(1,2)*pl_N));
pl_wL = sub2ind([pl_N, pl_N],floor(pl_ctr(2,1)*pl_N),floor(pl_ctr(2,2)*pl_N));


%%% Interpolation off set
pl_ax = (pl_ctr(1,1)*pl_N) - floor(pl_ctr(1,1)*pl_N);
pl_ay = (pl_ctr(1,1)*pl_N) - floor(pl_ctr(1,2)*pl_N);


%%% Mass Ratio Coefficients
% st_Mr = (pl_rho*pl_H*pl_h^2)*(1+k*pl_sigma0)/((pl_rho*pl_H*(pl_h^2)*(1+k*pl_sigma0))+(st_rho*st_A*st_h*(1+k*st_sigma0)))
% pl_Mr = (st_rho*st_A*st_h)*(1+k*st_sigma0)/((pl_rho*pl_H*(pl_h^2)*(1+k*pl_sigma0))+(st_rho*st_A*st_h*(1+k*st_sigma0)))

pl_Mr = 1/((pl_rho*pl_H*pl_h^2)*(1+k*pl_sigma0));
st_Mr = 1/((st_rho*st_A*st_h)*(1+k*st_sigma0));


%%% Spreading vector
if interpJ == 0
  J = sparse(zeros(Nt,1));
  % J([pl_w0 ss+1]) = [pl_Mr -st_Mr];
  pl_J = sparse(ss,1);pl_J(pl_w0) = 1;
  st_J = sparse(st_N,1); st_J(1) = 1;
  J([pl_w0 ss+1]) = [pl_Mr -st_Mr];
end

if interpJ == 1
  %%% Linear interpolation indeces
  lii = [pl_w0, pl_w0+1, pl_w0+pl_N,pl_w0+pl_N+1];
  %%% Linear interpolation coeffcients
  lic = [(1-pl_ax)*(1-pl_ay),(1-pl_ax)*pl_ay,pl_ax*(1-pl_ay),pl_ax*pl_ay];
  J = sparse(Nt,1);
  J([lii, ss+1]) = [pl_Mr*lic, -st_Mr]

  pl_J = sparse(ss,1); pl_J(lii) = lic;
  st_J = sparse(st_N,1); st_J(1) = 1;

end

Mr = inv( 1/(st_rho*st_A*st_h*(1+k*st_sigma0)) +...
      (1/(pl_rho*pl_H*(pl_h^2)*(1+k*pl_sigma0)))*(pl_J'*pl_J));

% Mr = (pl_rho*pl_H*pl_h^2)*(st_rho*st_A*st_h)*(1+k*st_sigma0)*(1+k*pl_sigma0)/(()+(st_rho*st_A*st_h*(1+k*st_sigma0)))

%%% Coupling Vector
% F = [-pl_mB(pl_w0,:),st_mB(1,:)];        % coupling vector for B matrix
% F = [J*B];        % coupling vector for B matrix
F = Mr*[-pl_J'*pl_mB,st_J'*st_mB];        % coupling vector for B matrix
% F1 = Mr*[(2*pl_mI(pl_w0,:) + (2*pl_sigma1*k/pl_h^2)*LA(pl_w0,:))/(1+k*pl_sigma0), -2*st_mI(1,:)/(1+k*st_sigma0)];   % coupling vector for C matrix
% F1 = Mr*[(2*pl_J'*pl_mI + (2*pl_sigma1*k/pl_h^2)*pl_J'*LA)/(1+k*pl_sigma0), -2*st_J'*st_mI/(1+k*st_sigma0)];   % coupling vector for C matrix
F1 = Mr*[pl_J'*(2*pl_mI + (2*pl_sigma1*k/pl_h^2)*LA)/(1+k*pl_sigma0),...
        -2*st_J'*st_mI/(1+k*st_sigma0)];   % coupling vector for C matrix

%%% Main Sparisty Matrices
B = blkdiag(pl_mB,st_mB) + J*F;
C = blkdiag(pl_mC,st_mC) + J*F1;


%%%% Add coupling to Matrices
% B = B + J*F;
% B(pl_w0,:) = B(pl_w0,:) + pl_Mr*F;
% B(ss+1,: ) = B(ss+1,: ) - st_Mr*F;

% OR
%%%careful of signs, exploding systems tend to be due to backward time step
% C(pl_w0,:) = C(pl_w0,:) + pl_Mr*F1;
% C(ss+1,: ) = C(ss+1,: ) - st_Mr*F1;
% C = C + J*F1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialise I/O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% initialise output
u = zeros(ss,1);
w = zeros(st_N,1);

%%% Joined vectors
uw = [u;w];
uw1 = uw;
uw2 = uw;

%%%% input
fvect = [u;w];                     % force vector

%%%% output
y = zeros(Nf,1);


%%%% Energy %%%NOTE%%% Following NSS, still not sure about Potential energy
if energyAn
  % energy coefficients = [Kinetic, Potential];
  pl_coE = [.5*(pl_h^2)/k^2 , .5*(pl_kappa^2)/(pl_h^2)];
  st_coE = [.5*(st_h)/k^2 , .5*(st_c^2)/st_h];

  pEnergy = zeros(Nf,1); % total
  pKE = pEnergy;          % kineteic
  pPE = pEnergy;          % potential

  sEnergy = zeros(Nf,1); % total
  sKE = sEnergy;          % kineteic
  sPE = sEnergy;          % potential

  fEnergy = sEnergy;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('All variables initialised')

%%%%% Main loop
tic
if run
  for n = 1:Nf

    %%% update input forces
    fvect(ss + li) = st_f(n);

    %%% main operation
    uw = B*uw1 + C*uw2 + fvect;


    % read output
    if (outtype==1)
      y(n,:) = uw(lo);

    elseif (outtype==2)
      y(n,:) = (SR*(uw(lo)-uw1(lo)));

    end

    if plot_on
      if ~mod(n,FrameDrop)

        figure(1)
        %%% plate
        subplot(2,1,1)
        mesh(reshape(uw(1:ss) ,pl_N,pl_N));
        zlim([-pTr pTr])
        % shading interp
        % view(2)
        % axis square

        %%% string
        subplot(2,1,2)
        plot([1:st_N]'*st_h, uw(ss+1:end), 'k');
        axis([0 st_L -sTr sTr])

        % surf(sZ+3,sY+3,sX+uw(ss+1:end))
        % axis([0 st_L -sTr sTr -sTr sTr])
        % axis square
        set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen

        % F(n) = getframe(gca);

        drawnow
      end
    end

    if energyAn % shift energy state
      u = uw(1:ss); u1 = uw1(1:ss);
      w = uw(ss+1:end); w1 = uw1(ss+1:end);

      pKE(n) = pl_coE(1)*((u-u1)'*(u-u1)); pPE(n) = pl_coE(2)*((LA*u)' * (LA*u1));
      pEnergy(n) = pl_coE*[sum(((u-u1).^2)); (LA*u)' * (LA*u1)];

      sKE(n) = st_coE(1)*((w-w1)'*(w-w1)); sPE(n) = st_coE(2)*((Dn*w)' * (Dn*w1));
      sEnergy(n) = st_coE*[sum(((w-w1).^2)); (Dn*w)' * (Dn*w1)];
      fEnergy(n) = (k*st_c^2/(2*st_h^2))*(w(2)-w(1)) + SR*(w(1)-w1(1)) + (pl_kappa^2*k*.5)*pl_mB(pl_w0,:)*u - SR*(u(pl_w0) - u1(pl_w0));
    end

    %%% shift state
    uw2 = uw1; uw1 = uw;

  end
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis Module for Thin Plate and string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plate_string_analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if benchtest
  profile viewer
  profile off
end

% EOF
