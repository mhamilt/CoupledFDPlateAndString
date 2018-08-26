%%%%% FDTD 2D Wave Model
%%%%% Matthew Hamilton s0674653
%%%%% Description:
%%%%%
%%%%% A Kirchhoff Thin Plate FDTD Model

clear all
close all hidden

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Flags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% EDIT THESE %%%%%%%

% Conditions
bctype = 1;          % boundary condition type: 1: simply supported, 2: clamped
outtype = 2;         % output type: 1: displacement, 2: velocity
losstype = 2;        % loss type: 1: independant, 2: dependant

% Scheme Run Options
plot_on = true;     % plot output
play = false;        % play out audio at the end
run = true;          % run scheme
skpErr = true;       % skip error checking for debugging

% Analysis
modeAn = true;        % run modal anlysis
energyAn = false;     % run energy analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% EDIT THESE %%%%%%%

% simulation
Tf = .125;                   % duration
nu = .4;                     % Poisson Ratios (< .5)
ctr = [.5, .5];            % centre point of excitation as percentage
wid = .2;                    % width (m or %)?
u0 = 0; v0 = 10;              % excitation displacement and velocity
rp = [.125, .825; .825, .125];   % readout position as percentage on grid.

% physical parameters
E = 2e11;                    % Young's modulus
rho = 7850;                  % density (kg/m^3)

H = .005;                     % thickness (m)
L = 1;                       % plate length (m)
Lx = 2;                       % plate length X axis (m)
Ly = 2;                       % plate length Y axis (m)
loss = [100, 1; 1000, .5];    % loss [freq.(Hz), T60;...]

% I/O
OSR = 8;                     % Oversampling ratio
SR = 44.1e3;                 % sample rate (Hz)

% plotting
pTr = 5e-3;                   % max trans limit of of plot
FrameDrop = 100;
spin = 0;         % control spin of plotting
%%%%%%%%% EDIT THESE %%%%%%%

%%% oversample
SR = SR*OSR;                        % redefine SR by OSR



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Error Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
skpErr = true; % skip error checking for debugging

% Array of variable names and conditions they must meet.
paraLs = {'SR', 'bctype'};
paraCons = [(SR<1e3),(bctype == 2 && energyAn)];

if ~skpErr

  if any(paraCons)

    for paraError = paraLs(paraCons)

      switch char(paraError)

        %%% WARNING ERRORS %%%

      case 'bctype'
      warning([sprintf('Energy Analysis is only valid\n'),...
      sprintf('for simply supported conditions\n'),...
      sprintf('Energy Analysis Off\n')])
      energyAn = false;

        %%% FATAL ERRORS %%%

        case 'SR'
        error([sprintf('SR is too low\n'),...
        sprintf('Make sure it is > 1000 Hz\n')])
      otherwise
        error([sprintf('The following variables are invalid\n\n')...
        sprintf('''%s'' \n',paraLs{paraCons}) ...
        sprintf('\nCheck their values and try again\n')])

      end

    end

  end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Motion Coefficients
D = (E*(H)^3)/(12*(1-(nu^2)));
kappa = sqrt(D / (rho*  H) );

%%%%% Scheme Spacing
k = 1/SR;                       % time step
hmin = 2*sqrt(k*kappa);         % NSS (equation 12.5)
Nx = floor((Lx)/hmin);          % number of segments
Ny = floor((Ly)/hmin);          % number of segments
h = sqrt(Lx*Ly/(Nx*Ny));        % adjusted grid spacing
Nf = floor(SR*Tf);              % number of time steps
Nx = Nx+1;                      % grid points x
Ny = Ny+1;                      % grid points y
ss = Nx*Ny;                     % total grid size.

mu = (kappa * k)/(h^2);         % scheme parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % Loss coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if losstype ==0 % lossless

  sigma0 = 0;
  sigma1 = 0;
end

if losstype ==1 % frequency independant loss

  sigma0 = 6*log(10)/loss(1,2);
  sigma1 = 0;
end

if losstype == 2 % frequency dependant loss

  z1 = 2*kappa*(2*pi*loss(1,1))/(2*kappa.^2);
  z2 = 2*kappa*(2*pi*loss(2,1))/(2*kappa.^2);

  sigma0 = 6*log(10)*(-z2/loss(1,2) + z1/loss(2,2))./(z1-z2);
  sigma1 = 6*log(10)*(1/loss(1,2) - 1/loss(2,2))./(z1-z2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % Read In/Out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lo = rp.*[Nx Ny;Nx Ny];
lo = floor(sub2ind([Nx Ny],lo(1), lo(3)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Force Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y] = meshgrid([1:Nx]*h, [1:Ny]*h);         % Grid of point in value of meters
dist = sqrt((X-(ctr(1)*Lx)).^2 + (Y-(ctr(2)*Ly)).^2); % distance of points from excitation
ind = sign(max(-dist+(wid*0.5),0));         % displaced grid points (logical)
rc = .5*ind.*(1+cos(2*pi*dist/wid));        % displacement
rca = rc;
rc = rc(:);                                 % 2D plane as vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficient Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BH = biharm2(Ny,Nx,bctype); % biharmonic matrix
LA = laplace2(Ny,Nx,bctype); % Laplacian matrix

A = (1/(1+k*sigma0))*speye(ss);   %NOTE% Currently inverted
B = (-(mu^2)*BH + (2*sigma1*k/(h^2))*LA + 2*speye(ss)) * A;
C = (-(2*sigma1*k/(h^2))*LA - (1-sigma0*k)*speye(ss))  * A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialise I/O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% initialise vectors
u2 = u0*rc;
u1 = (u0+(k*v0))*rc;
u  = u2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_on
  fig = figure(1);
  fig.Color = [0 0 0];
  ax = axes('XLim',[0 Nx],'YLim',[0 Ny],'ZLim',[-pTr pTr]);
  colormap copper
  shading interp
  material SHINY
%   set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Energy Coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if energyAn
  % energy coefficients = [Kinetic, Potential];
  coE = [.5*(h^2)/k^2 , .5*(kappa^2)*(1/h^2)];
  Energy = zeros(Nf,1); % total
  KE = Energy;          % kineteic
  PE = Energy;          % potential
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Main loop
if run
  for n = 1:Nf

    % Main Operation
    u = B*u1 + C*u2;

    % read output
    if (outtype==1)
      y(n,:) = u(lo);

    elseif (outtype==2)
      y(n,:) = (SR*(u(lo)-u1(lo)));

    end

    % Shift State
    u2 = u1; u1 = u;

    if plot_on
      if ~mod(n,FrameDrop)

        %%% plate
        surf(ax,reshape(u,Ny,Nx),'FaceLighting','phong','FaceColor','interp',...
        'AmbientStrength',0.5);
        light('Position',[1 0 1],'Style','local')
        shading interp
        view(45+spin*n/(Nf*FrameDrop),30)
        caxis([-1e-5 3.5e-5])
        ax.XLim =[0 max([Nx Ny])]; ax.YLim =[0 max([Nx Ny])]; ax.ZLim =[-pTr pTr];
        axis off

        drawnow

      end
    end

    %%% Energy Calculation
    if energyAn % shift energy state
      KE(n) = coE(1)*sum((u1-u2).^2); PE(n) = coE(2)*((LA*u1)' * (LA*u2));
      Energy(n) = coE*[sum((u1-u2).^2); (LA*u1)' * (LA*u2)];
    end

  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plate_analysis

% EOF
