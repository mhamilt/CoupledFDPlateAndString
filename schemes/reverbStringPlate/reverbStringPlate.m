%%%%% FDTD 2D Wave Model
%%%%% Matthew Hamilton s0674653
%%%%% Description:
%%%%%
%%%%% A Couple String and Plate Reverb Effect

clear all
close all hidden

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Instrument File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reverbParams % file containing all relevent parameters for this set of schemes

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
pl_Nx = floor((pl_Lx)/pl_hmin);          % number of segments
pl_Ny = floor((pl_Ly)/pl_hmin);          % number of segments
pl_h = sqrt(pl_Lx*pl_Ly)/sqrt(pl_Nx*pl_Ny);    % adjusted grid spacing
pl_mu = (pl_kappa * k)/(pl_h^2);    % scheme parameter
Nf = floor(SR*Tf);                  % number of time steps
pl_Nx = pl_Nx +1;
pl_Ny = pl_Ny +1;
ss = pl_Nx*pl_Ny;                     % total grid size.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% String Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% derived parameters

st_A = pi*st_r.^2;                       % string cross-sectional area
st_I = 0.25*pi*st_r.^4;                  % string moment of intertia

st_c = sqrt(st_T./(st_rho.*st_A));        % wave speed
st_K = sqrt(st_E.*st_I./(st_rho*st_A));   % stiffness constant

%%%%% grid

st_hmin = sqrt(0.5* (st_c.^2*k^2+sqrt(st_c.^4*k^4+16*st_K.^2.*k.^2)) );
st_N = floor(st_L./st_hmin);          % number of segments (N+1 is number of grid points)
st_h = st_L./st_N;                    % adjusted grid spacing
st_lambda = st_c*k./st_h;             % Courant number
st_mu = st_K*k./st_h.^2;               % numerical stiffness constant
st_N = st_N+1;                       % change st_N to be number of grid points

%%%%% I/O
li = floor(st_xi*st_N);         % grid index of excitation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Error Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Array of variable names and conditions they must meet.

if ~skpErr
  paraLs = {'SR' 'strNum'};
  paraCons = [(SR<1e3), strNum~=length(unique(floor(pl_ctr*pl_Nx)))];

  if any(paraCons)

    for paraError = paraLs(paraCons)

      switch char(paraError)

        %%% WARNING ERRORS %%%
        case 'strNum'
        error([sprintf('two string share the same grid point, make pl_N bigger')])

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
%%%% % String Loss coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if losstype ==0
  st_sigma0 = zeros(1,strNum);
  st_sigma1 = zeros(1,strNum);

  pl_sigma0 = 0;
  pl_sigma1 = 0;
end

if losstype ==1
  % frequency independant loss
  st_sigma0 = (6*log(10)/st_loss(1,2))*ones(1,strNum);
  st_sigma1 = zeros(1,strNum);

  pl_sigma0 = 6*log(10)/pl_loss(1,2);
  pl_sigma1 = 0;
end

if losstype == 2

  %%% String Loss
  st_z1 = (-st_c.^2 + sqrt(st_c.^4 + 4*st_K.^2.*(2*pi*st_loss(1,1)).^2))./(2*st_K.^2);
  st_z2 = (-st_c.^2 + sqrt(st_c.^4 + 4*st_K.^2.*(2*pi*st_loss(2,1)).^2))./(2*st_K.^2);

  st_sigma0 = 6*log(10)*(-st_z2./st_loss(1,2) + st_z1./st_loss(2,2))./(st_z1-st_z2);
  st_sigma1 = 6*log(10)*(1./st_loss(1,2) - 1./st_loss(2,2))./(st_z1-st_z2);

  %%% Plate Loss
  pl_z1 = 2*pl_kappa*(2*pi*pl_loss(1,1))/(2*pl_kappa.^2); %from 1D
  pl_z2 = 2*pl_kappa*(2*pi*pl_loss(2,1))/(2*pl_kappa.^2);

  pl_sigma0 = 6*log(10)*(-pl_z2/pl_loss(1,2) + pl_z1/pl_loss(2,2))./(pl_z1-pl_z2);
  pl_sigma1 = 6*log(10)*(1/pl_loss(1,2) - 1/pl_loss(2,2))./(pl_z1-pl_z2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % Read In/Out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st_lo = floor(st_xo*st_N);         % grid index of output


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

if any([(st_xo*st_L)<st_h (st_L-(st_xo*st_L))<st_h])
  warning([sprintf('xo is too close to string boundary\n'),...
  sprintf('changing xo to equal h\n')]);
  st_xo = abs(round(st_xo)-st_h);
  if any(st_lo<=2)
    st_lo(st_lo<=2) = 3;
  end
  if any(st_lo>=st_N)
    st_lo(st_lo>=st_N) = st_N(st_lo>=st_N)-1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Force Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pl_d0 = (k^2)/(pl_rho*pl_H*(pl_h^2))*(1/(1+k*pl_sigma0)); % input force coefficient
force = [force;zeros(zpad*SR,1)]*pl_d0;       % input oscillator
pl_fi = floor(((f_ctr(1)*pl_Nx)-1)*pl_Nx + f_ctr(2)*pl_Ny);    % linear index for force

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plate Coefficient Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LA = laplace2(pl_Ny,pl_Nx,bctype); % Laplacian Matrix
BH = biharm2(pl_Ny,pl_Nx,bctype); % biharmonic matrix
pl_mI = fidimat(pl_Ny,pl_Nx,'I');

pl_mA = (1/(1+k*pl_sigma0))*speye(ss);
pl_mB = (-(pl_mu^2)*BH + 2*pl_mI + (2*k*pl_sigma1/pl_h^2)*LA) * pl_mA;
pl_mC = (-(1-pl_sigma0*k)*pl_mI - (2*k*pl_sigma1/pl_h^2)*LA) * pl_mA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% String Coefficients Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

countN = zeros(1,strNum+1);

for s = 1:strNum
  %%%%% scheme coefficients
  Dn = fidimat(st_N(s),'x+',bctype);
  Dnn = fidimat(st_N(s),'xx',1);
  Dnnnn = fidimat(st_N(s),'xxxx',st_bctype);
  Dnnnn(2,1) = -2;
  st_mI = speye(st_N(s));

  %%% String sparsity matrix, uncoupled.
  st_mA = (1/(1+k*st_sigma0(s)));
  st_mB{s} = ( (st_lambda(s)^2)*Dnn - (st_mu(s)^2)*Dnnnn + (2*st_mI) +...
  ((2*st_sigma1(s)*k/st_h(s)^2)*Dnn)) * st_mA ;
  st_mC{s} = -((1-st_sigma0(s)*k)*st_mI + (2*st_sigma1(s)*k/st_h(s)^2)*Dnn) * st_mA;

  %%%%% Alter matrices for coupling string boundary
  st_mB{s}([1 end],:) = 0;st_mC{s}([1 end],:) = 0;
  st_mB{s}(:,end) = 0; st_mC{s}(:,end) = 0;
  st_mB{s}(1,:) = ((st_lambda(s)^2)*Dn(1,:) - (st_mu(s)^2)*Dnn(2,:) + (2*st_mI(1,:))) * st_mA;
  st_mC{s}(1,:) = -((1-st_sigma0(s)*k)*st_mI(1,:));

  countN(s+1) = countN(s)+st_N(s); % increment count of grid size

end

st_mB = blkdiag(st_mB{:});
st_mC = blkdiag(st_mC{:});
st_mI = speye(sum(st_N));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Couple Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% coupling points
pl_w0 = sub2ind([pl_Nx, pl_Ny],floor(pl_ctr*pl_Nx),floor(.25*pl_Ny)*ones(1,strNum));
% pl_wL = sub2ind([pl_N, pl_N],floor(pl_ctr(2,1)*pl_N),floor(pl_ctr(2,2)*pl_N));

%%% NOTE Very long winded an in the process of being tidied up

%%% Mass Ratio Coefficients
st_Mr = (pl_rho*pl_H*pl_h^2)*(1+k*pl_sigma0)./((pl_rho*pl_H*(pl_h^2)*(1+k.*pl_sigma0))+(st_rho.*st_A.*st_h.*(1+k.*st_sigma0)));
pl_Mr = (st_rho*st_A.*st_h).*(1+k*st_sigma0)./((pl_rho*pl_H*(pl_h^2)*(1+k.*pl_sigma0))+(st_rho.*st_A.*st_h.*(1+k.*st_sigma0)));

%%% Spreading vector
J = sparse(zeros(ss,1)); J(pl_w0) = 1;

% Main Sparisty Matrices
B = blkdiag(pl_mB,st_mB);
C = blkdiag(pl_mC,st_mC);

%%% Coupling Vector
for s = 1:strNum

  %%% Coupling Vector
  cpl_v = [-pl_mB(pl_w0(s),:),st_mB(1+countN(s),:)];        % coupling vector for B matrix
  cpl_vc = [(2*pl_mI(pl_w0(s),:) + (2*pl_sigma1*k/pl_h^2)*LA(pl_w0(s),:))/(1+k*pl_sigma0), -2*st_mI(1+countN(s),:)/(1+k*st_sigma0(s))];   % coupling vector for C matrix

  % %%%% Add coupling to Matrices
  B(pl_w0(s),:) = B(pl_w0(s),:) + pl_Mr(s)*cpl_v;
  B(ss+1+countN(s),:) = B(ss+1+countN(s),: ) - st_Mr(s)*cpl_v;
  %
  % %%%careful of signs, exploding systems tend to be due to backward time step
  C(pl_w0(s),:) = C(pl_w0(s),:) + pl_Mr(s)*cpl_vc;
  C(ss+1+countN(s),: ) = C(ss+1+countN(s),: ) - st_Mr(s)*cpl_vc;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialise I/O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% initialise output
u = zeros(ss,1);
w = zeros(sum(st_N),1);

%%% Joined vectors
uw = [u;w];
uw1 = uw;
uw2 = uw;

%%%% input
fvect = [u;w];                     % force vector
li = li + countN(1:end-1);


%%%% output
y = zeros(Nf,2);
st_lo = ss+st_lo+countN(1:end-1)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    fvect(li) = force(n);

    %%% main operation
    uw = B*uw1 + C*uw2 + fvect;


    % read output
    if (outtype==1)
      y(n,:) = panLR*uw(st_lo);

    elseif (outtype==2)
      y(n,:) = panLR*(SR*(uw(st_lo)-uw1(st_lo)));

    end


    if plot_on
      if ~mod(n,FrameDrop)

        figure(1)
        %%% plate
        subplot(2,1,1)
        surf(reshape(uw(1:ss) ,pl_Ny,pl_Nx));
        axis([0 max([pl_Nx pl_Ny]) 0 max([pl_Nx pl_Ny]) -pTr pTr]);
        % shading interp
        % view(2)
        axis square

        %%% string
        subplot(2,1,2)
        % plot([1:sum(st_N)]'*st_h, uw(ss+1:end), 'k');
        % axis([0 sum(st_L) -sTr sTr])
        plot(uw(ss+1:end), 'k');
        ylim([-sTr sTr])
        line([countN(:), countN(:)], .25*[-sTr sTr],'Color',[0 0 0],'LineWidth',2)
        % surf(sZ+3,sY+3,sX+uw(ss+1:end))
        % axis([0 st_L -sTr sTr -sTr sTr])
        % axis square
        set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen

        % F(n) = getframe(gca);

        drawnow
      end
    end


    %%% shift state
    uw2 = uw1; uw1 = uw;

    %%% upadate progress bar
    if ~mod(n,loadDrop)
      waitbar(n*oNf, loadBar);
    end

    %%%% Kill the script if it blows up
    if any(isnan(uw))
      error('UNSTABLE script terminated')
    end

  end

  close(loadBar)
end

toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis Module for Thin Plate and string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = decimate(y,OSR);
save_file = sprintf('%s_stringVerb.wav',sfile);
audiowrite(save_file, (y/abs(max(y(:)))), SR/OSR);

if benchtest
  profile viewer
  profile off
end

% EOF
