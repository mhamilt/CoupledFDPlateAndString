%%%%% FDTD 2D Wave Model
%%%%% Matthew Hamilton s0674653
%%%%% Description:
%%%%%
%%%%% A Coupled Kirchhoff Thin Plate and Stiff String FDTD Model

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Instrument File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plate_and_string_vars % file containing all parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plate Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Motion Coefficients
pl_D = (pl_E*(pl_H)^3)/(12*(1-(nu^2)));
pl_kappa = sqrt(pl_D / (pl_rho*pl_H) );

%%%%% Scheme Spacing
k = 1/SR;                           % time step
pl_hmin = 2*sqrt(k*pl_kappa);       % NSS (equation 12.5)
pl_Nx = floor((pl_Lx)/pl_hmin);          % number of segments
pl_Ny = floor((pl_Ly)/pl_hmin);          % number of segments
pl_h = sqrt(pl_Lx*pl_Ly)/sqrt(pl_Nx*pl_Ny);    % adjusted grid spacing
pl_mu = (pl_kappa * k)/(pl_h^2);    % scheme parameter
Nf = floor(SR*Tf);                  % number of time steps
pl_Nx = pl_Nx +1;
pl_Ny = pl_Ny +1;
ss = pl_Nx*pl_Ny;                     % total grid size.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% String Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% derived parameters
st_A = pi*st_r^2;                       % string cross-sectional area
st_I = 0.25*pi*st_r^4;                  % string moment of intertia
st_c = sqrt(st_T/(st_rho*st_A));        % wave speed
st_K = sqrt(st_E*st_I/(st_rho*st_A));   % stiffness constant

%%%%% grid
st_hmin = sqrt(0.5* (st_c.^2*k^2+sqrt(st_c.^4*k^4+16*st_K.^2.*k.^2)) );
st_N = floor(st_L/st_hmin);          % number of segments (N+1 is number of grid points)
st_h = st_L/st_N;                    % adjusted grid spacing
st_lambda = st_c*k/st_h;             % Courant number
st_mu = st_K*k/st_h^2;               % numerical stiffness constant
st_N = st_N+1;                       % change st_N to be number of grid points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % String Loss coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if losstype ==0  %lossless
  st_sigma0 = 0;  st_sigma1 = 0;
  pl_sigma0 = 0;  pl_sigma1 = 0;
end

if losstype ==1
  % frequency independant loss
  st_sigma0 = 6*log(10)/st_loss(1,2);
  st_sigma1 = 0;

  pl_sigma0 = 6*log(10)/pl_loss(1,2);
  pl_sigma1 = 0;
end

if losstype == 2
  %%% String
  st_z1 = (-st_c.^2 + sqrt(st_c.^4 + 4*st_K.^2.*(2*pi*st_loss(1,1))^2))./(2*st_K.^2);
  st_z2 = (-st_c.^2 + sqrt(st_c.^4 + 4*st_K.^2.*(2*pi*st_loss(2,1))^2))./(2*st_K.^2);
  st_sigma0 = 6*log(10)*(-st_z2/st_loss(1,2) + st_z1/st_loss(2,2))./(st_z1-st_z2);
  st_sigma1 = 6*log(10)*(1/st_loss(1,2) - 1/st_loss(2,2))./(st_z1-st_z2);

  %%% Plate
  pl_z1 = 2*pl_kappa*(2*pi*pl_loss(1,1))/(2*pl_kappa.^2); %from 1D
  pl_z2 = 2*pl_kappa*(2*pi*pl_loss(2,1))/(2*pl_kappa.^2);
  pl_sigma0 = 6*log(10)*(-pl_z2/pl_loss(1,2) + pl_z1/pl_loss(2,2))./(pl_z1-pl_z2);
  pl_sigma1 = 6*log(10)*(1/pl_loss(1,2) - 1/pl_loss(2,2))./(pl_z1-pl_z2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % Read In/Out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lo = rp.*[pl_Nx pl_Ny];
lo = [floor(sub2ind([pl_Nx pl_Ny],lo(1), lo(3))), floor(sub2ind([pl_Nx pl_Ny],lo(2), lo(4)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IV Read In/Out Error Check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Force Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st_f = zeros(Nf,1);                        % input force signal
durint = floor(st_dur*SR);                 % duration of force signal, in samples
exc_st_int = (floor(st_exc_st*SR))+1;      % start time index for excitation
durf = exc_st_int:exc_st_int+durint-1;     % sample values of force
st_d0 = (k^2/(st_h*st_rho*st_A));          % force foefficient

%%% reevaluate relevant points in force vector
st_f(durf) = st_famp*0.5*(1-cos((2/itype)*pi.*(durf/durint)))*st_d0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plate Coefficient Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LA = laplace2(pl_Ny,pl_Nx,bctype); % Laplacian Matrix
BH = biharm2(pl_Ny,pl_Nx,bctype); % biharmonic matrix
pl_mI = fidimat(pl_Ny,pl_Nx,'I');

pl_mA = (1/(1+k*pl_sigma0))*speye(ss);
pl_mB = (-(pl_mu^2)*BH + 2*pl_mI + (2*k*pl_sigma1/pl_h^2)*LA) * pl_mA;
pl_mC = (-(1-pl_sigma0*k)*pl_mI - (2*k*pl_sigma1/pl_h^2)*LA) * pl_mA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% String Coefficients Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% scheme coefficients
Dnp = fidimat(st_N,'x+');
Dnm = fidimat(st_N,'x-');
Dnn = fidimat(st_N,'xx');
Dnnnn = fidimat(st_N,'xxxx',st_bctype);
Dnnnn(2,1) = -2; Dnnnn(st_N-1,st_N) = -2;
st_mI = speye(st_N);

%%% String sparsity matrix, uncoupled.
st_mA = (1/(1+k*st_sigma0));
st_mB = ((st_lambda^2)*Dnn - (st_mu^2)*Dnnnn + (2*st_mI) + ((2*st_sigma1*k/st_h^2)*Dnn)) * st_mA ;
st_mC = -((1-st_sigma0*k)*st_mI + ((2*st_sigma1*k/st_h^2)*Dnn) ) * st_mA;

%%%%% Alter matrices for coupling string boundary
st_mB([1 end],:) = 0;st_mC([1 end],:) = 0;
st_mB([1 st_N],:) = ((st_lambda^2)*[Dnp(1,:);Dnm(st_N,:)] - (st_mu^2)*Dnn([2 st_N-1],:) + (2*st_mI([1 st_N],:))) * st_mA;
st_mC([1 st_N],:) = -((1-st_sigma0*k)*st_mI([1 st_N],:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Couple Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Total Number of Grid points
Nt = ss+ st_N;

%%% floored coupling points
pl_w0 = sub2ind([pl_Nx, pl_Ny],floor(pl_ctr(1,1)*pl_Nx),floor(pl_ctr(1,2)*pl_Ny));
pl_wL = sub2ind([pl_Nx, pl_Ny],floor(pl_ctr(2,1)*pl_Nx),floor(pl_ctr(2,2)*pl_Ny));

%%% Mass Ratio Coefficients
pl_Mr = 1/((pl_rho*pl_H*pl_h^2)*(1+k*pl_sigma0));
st_Mr = 1/((st_rho*st_A*st_h)*(1+k*st_sigma0));

%%% Spreading vector
if interpJ == 0
  %%% Spreading operators
  J = sparse(zeros(Nt,2));
  J([pl_w0 pl_wL+Nt ss+1 2*Nt]) = [pl_Mr -pl_Mr -st_Mr st_Mr];
  pl_J = sparse(ss,2);pl_J([pl_w0 pl_wL+Nt]) = 1;
  st_J = sparse(st_N,2); st_J([1 2*st_N]) = 1;
end

if interpJ == 1
  pl_w0 = [pl_w0; pl_wL];

  %%% Interpolation off set
  pl_ax = (pl_ctr(:,1)*pl_Nx) - floor(pl_ctr(:,1)*pl_Nx);
  pl_ay = (pl_ctr(:,2)*pl_Ny) - floor(pl_ctr(:,2)*pl_Ny);

  %%% Linear interpolation indeces
  lii = [pl_w0, pl_w0+1, pl_w0+pl_Ny,pl_w0+pl_Ny+1];
  %%% Linear interpolation coeffcients
  lic = [(1-pl_ax).*(1-pl_ay),(1-pl_ax).*pl_ay,pl_ax.*(1-pl_ay),pl_ax.*pl_ay];

  %%% Spreading operators
  J = sparse(Nt,2); J([lii+[0;1]*Nt]) = [pl_Mr*lic(1,:);-pl_Mr*lic(2,:)];
  J(st_w0+ss+[0;1]*Nt) = [-st_Mr;st_Mr];
  pl_J = sparse(ss,2); pl_J(lii+[0;1]*ss) = lic;
  st_J = sparse(st_N,2); st_J([1 2*st_N]) = 1;
end

%%% Mass ratio
Mr = inv( 1/(st_rho*st_A*st_h*(1+k*st_sigma0)) +...
      (1/(pl_rho*pl_H*(pl_h^2)*(1+k*pl_sigma0)))*(pl_J'*pl_J));

%%% Coupling Vector
F = Mr*[-pl_J'*pl_mB,st_J'*st_mB];        % coupling vector for B matrix
F1 = Mr*[pl_J'*(2*pl_mI + (2*pl_sigma1*k/pl_h^2)*LA)/(1+k*pl_sigma0),...
        -2*st_J'*st_mI/(1+k*st_sigma0)];   % coupling vector for C matrix

%%% Main Sparisty Matrices
B = blkdiag(pl_mB,st_mB) + J*F;
C = blkdiag(pl_mC,st_mC) + J*F1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialise I/O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% initialise output
u = zeros(ss,1);
w = zeros(st_N,1);

%%% Joined vectors
uw = [u;w];
uw1 = uw;
uw2 = uw;

%%%% input
fvect = [u;w];                % force vector
li = floor(st_xi*st_N) + ss;  % grid index of excitation

%%%% output
y = zeros(Nf,2);


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('All variables initialised')

%%%%% Main loop
tic
if run
  %%%%% Load Bar
  loadBar = waitbar(0, 'Processing...');
  loadDrop = floor(Nf*.01);
  oNf = 1/Nf;

  for n = 1:Nf

    %%% update input forces
    fvect(li) = st_f(n);

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
        surf(reshape(uw(1:ss) ,pl_Ny,pl_Nx));
        axis([0 max([pl_Nx pl_Ny]) 0 max([pl_Nx pl_Ny]) -pTr pTr]);

        %%% string
        subplot(2,1,2)
        plot(uw(ss+1:end), 'k');
        ylim([-sTr sTr])

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

    if ~mod(n,loadDrop)
      waitbar(n*oNf, loadBar);
    end

  end

end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plate_string_analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EOF
