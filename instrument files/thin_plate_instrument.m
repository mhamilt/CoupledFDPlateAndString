%%%%% FDTD 2D Wave Model
%%%%% Matthew Hamilton s0674653
%%%%% Description:
%%%%%
%%%%% Instrument file for 2D Plate
%%%%% Contains all the parameters to run the FDTD thin plate scheme


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Flags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% EDIT THESE %%%%%%%

% Conditions
bctype = 1;          % boundary condition type: 1: simply supported, 2: clamped
outtype = 1;         % output type: 1: displacement, 2: velocity
losstype = 0;        % loss type: 1: independant, 2: dependant
scale    = 0;        %  scaling: 0: Scaled  1: unscaled

% Scheme Run Options
plot_on = true;     % plot output
play = false;        % play out audio at the end
run = true;          % run scheme
skpErr = true;       % skip error checking for debugging

% Analysis
modeAn = false;        % run modal anlysis
energyAn = true;     % run energy analysis
benchtest = false;   % turn benchtesting on to view MATLAB profiler

%%%%%%%%% EDIT THESE %%%%%%%
%%START BENCH TEST%%
if benchtest
  profile on
end
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% EDIT THESE %%%%%%%

% simulation
Tf = 10;                      % duration
nu = .5;                     % Poisson Ratios (< .5)
nu_xy = .5;
nu_yx = .5;
G_xy = 13;                  % shear modulus
ctr = [.45, .45];              % centre point of excitation as percentage
wid = .2;                    % width (m or %)?
u0 = 0; v0 = 1;              % excitation displacement and velocity
% rp = [.35, .65; .85, .15];              % readout position as percentage on grid.
rp = [.45, .65];              % readout position as percentage on grid.
% physical parameters
% E = 2e11;                    % Young's modulus
% rho = 7850;                  % density (kg/m^3)

% % wood
E = 11e9;
rho = 480;
% aniso = 2;                 % anisotropy ratio
% Ey = Ex/aniso;

H = .005;                     % thickness (m)
L = .5;                       % plate length (m)
Lx = .7;                       % plate length X axis (m)
Ly = .7;                       % plate length Y axis (m)

loss = [100, 1; 1000, .5];    % loss [freq.(Hz), T60;...]

% I/O
OSR = 1;                     % Oversampling ratio
SR = 44.1e3;                 % sample rate (Hz)

SR = OSR*SR;
% Force
fHz = 258.95*2;                   % frequency of force input


% plotting
pTr = 5e-4;                   % max trans limit of of plot
FrameDrop = 1;
spin = 2e5;         % control spin of plotting


fmode = 396.43;
zpad = loss(2,2)/2;
