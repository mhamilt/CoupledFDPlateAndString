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

st_f = zeros(Nf,1);                        % input force signal
durint = floor(st_dur*SR);                 % duration of force signal, in samples
exc_st_int = (floor(st_exc_st*SR))+1;          % start time index for excitation
for n = 1:length(st_exc_st)
  hit = exc_st_int(n):exc_st_int(n)+durint-1;
  st_f(hit) = st_famp*0.5*(1-cos((2/itype)*pi.*...
  (hit-exc_st_int(n))/durint))*(k^2/(st_h*st_rho*st_A));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plate Coefficient Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LA = laplace2(pl_N,pl_N,bctype); % Laplacian Matrix
BH = biharm2(pl_N,pl_N,bctype); % biharmonic matrix
pl_I = speye(ss);

pl_mB = -(pl_mu^2)*BH + 2*fidimat(pl_N,pl_N,'I');
pl_mC = -pl_I;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% String Coefficients Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% scheme coefficients
Dn = fidimat(st_N,'x+',bctype);
Dnn = fidimat(st_N,'xx');
st_I = speye(st_N);



%%% String sparsity matrix, uncoupled.
st_mB = 2*speye(st_N)+(st_lambda^2)*Dnn;

% from definition of coupling point
% a forward difference instead of 0 or a continuation of Dnn
% w(1) = lambda^2*Dn(1,:)*w1 - tau_s*cpl_v
% st_mB(1,[1 2]) = st_lambda^2*[-1 1] + [2 0];
st_mB(1,1) = 1;
st_mC = -st_I;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Couple Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% coupling points
pl_w0 = sub2ind([pl_N, pl_N],floor(pl_ctr(1,1)*pl_N),floor(pl_ctr(1,2)*pl_N));
pl_wL = sub2ind([pl_N, pl_N],floor(pl_ctr(2,1)*pl_N),floor(pl_ctr(2,2)*pl_N));

%%% NOTE Very long winded an in the process of being tidied up


%%% Toque coefficients
tau_p = st_T*st_h/((pl_rho*pl_H*pl_h^2*st_c^2)+(st_T*st_h));
tau_s = st_c^2*pl_rho*pl_H*pl_h^2/((pl_rho*pl_H*pl_h^2*st_c^2)+(st_T*st_h));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Spreading vector
J = sparse(zeros(ss,1)); J(pl_w0) = 1;

% Main Sparisty Matrices
B = blkdiag(pl_mB,st_mB);
C = blkdiag(pl_mC,st_mC);

%%% Coupling Vector
cpl_v = [-pl_mB(pl_w0,:), st_mB(1,:)]; % coupling vector

%%%% Add coupling to Matrices
B(pl_w0,:) = B(pl_w0,:) + tau_p*cpl_v;
B(ss+1,:) = B(ss+1,:) - tau_s*cpl_v;

C(pl_w0,:) = C(pl_w0,:) + tau_p*[2*J', -2*st_I(1,:)];
C(ss+1,:) = C(ss+1,:) - tau_s*[2*J', -2*st_I(1,:)];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialise I/O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% initialise output
u = zeros(ss,1);
u2 = u;
u1 = u;

w  = zeros(st_N,1);
w1 = w;
w2 = w;

y = zeros(Nf,1);

%%%% input
fvect = [w];                     % force vector
% li = ss + li;

%%% Joined vectors
% wu = fvect;
% wu1 = fvect;
% wu2 = fvect;




%%%% Energy %%%NOTE%%% Following NSS, still not sure about Potential energy
if energyAn

  % energy coefficients = [Kinetic, Potential];
  coE = [.5*(h^2)/k^2 , .5*(pl_kappa^2)*(1/h^2)];

  Energy = zeros(Nf,1); % total
  KE = Energy;          % kineteic
  PE = Energy;          % potential

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
    fvect(li) = st_f(n);

    %%% main operation
    % wu = B*wu1 - C*wu2 + fvect;
    w(1) = st_lambda^2*Dn(1,:)*w1 + 2*w1(1) - w2(1) - ...
    tau_s*(st_lambda^2*Dn(1,:)*w1 + 2*w1(1) - 2*w2(1) + J'*(pl_mu^2*BH*u1 - 2*u1 + 2*u2 ));
    w = st_mB*w1 + st_mC*w2 + fvect;
    u = pl_mB*u1 + pl_mC*u2 + J*(st_lambda^2*Dn(1,:)*w1 + 2*w1(1) - 2*w2(1) + J'*(pl_mu^2*BH*u1 - 2*u1 + u2 ));


    % read output
    if (outtype==1)
      y(n,:) = u(lo);

    elseif (outtype==2)
      y(n,:) = (SR*(u(lo)-u1(lo)));

    end

    if plot_on
      if ~mod(n,FrameDrop)
        figure(1)
        subplot(2,1,1)
        % mesh(reshape(wu(1:ss) ,pl_N,pl_N));
        mesh(reshape(u ,pl_N,pl_N));
        zlim([-1e-4 1e-4])
        % shading interp
        % view(2)

        subplot(2,1,2)
        % plot([1:st_N]'*st_h, wu(ss+1:end), 'k');
        plot([1:st_N]'*st_h, w, 'k');
        axis([0 st_L -1e-3 1e-3])
        set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen
        % axis(view);
        % F(n) = getframe(gca);

        drawnow
      end
    end

    %%% shift state
    % wu2 = wu1; wu1 = wu;
    w2 = w1; w1 = w;
    u2 = u1; u1 = u;


    if energyAn % shift energy state
      KE(n) = coE(1)*sum((u1-u2).^2); PE(n) = coE(2)*((LA*u1)' * (LA*u2));
      Energy(n) = coE*[sum((u1-u2).^2); (LA*u1)' * (LA*u2)];
    end

  end
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (Housekeeping)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Modal Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculating plate modes

if modeAn
  p = [1:5];
  q = [1:5]';

  [P,Q] = meshgrid(p,q);

  mfreqs = (pi*pl_kappa/(2*L^2))*((P.^2) + (Q.^2));     % mode frequencies
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
font_s = 14; % font point size
set(0,'DefaultAxesFontSize',font_s);

% Output FFT

if modeAn
  YF = abs(fft(y));
  fy = linspace(0, SR, length(y)); % frequency axis
  figure(1);
  plot(fy, 20*log10(YF));
  ax1 = gca;
  line([mfreqs(:) mfreqs(:)],ax1.YLim,'Color',[1 0 0])
  ax1.XLim = [0, 1500];
  title('Scheme Output FFT', 'FontSize', font_s)
  xlabel('Freq._{Hz}');
  ylabel('magnitude');
  legend('FFT Output_{dB}','Mode Frequencies_{Hz}')
  set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen
end

if energyAn
  %change in energy
  dEdt = (Energy(2:end)-Energy(1:end-1));
  figure(2)
  subplot(2,1,1)
  plot(Energy);
  hold on
  plot(dEdt)
  plot(KE)
  plot(PE)
  ax2 = gca;
  title('Change in Energy at Time Step k', 'FontSize', font_s)
  xlabel('Sample_{k}');
  ylabel('Power_{Watts}');
  set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen
  legend('Total Energy','Power_{dE/dt}','Kinetic','Potential')


  subplot(2,1,2)
  plot(dEdt,'.')
  title('dE/dt rounding error', 'FontSize', font_s)
  xlabel('Sample_{k}');
end

% if plot_on
%   v = VideoWriter('2DWave','Motion JPEG AVI');
%   % v.CompressionRatio = 3;
%   v.FrameRate = 24;
%   open(v)
%   writeVideo(v,F(1:4000))
%   close(v)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if benchtest
  profile viewer
  profile off
end
