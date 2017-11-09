%%%%% FDTD 2D Wave Model
%%%%% Matthew Hamilton s0674653
%%%%% Description:
%%%%%
%%%%% Instrument file for 2D Plate and coupled strings
%%%%% Contains all the parameters to run the FDTD thin plate scheme


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Flags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Conditions
bctype = 2;          % boundary condition type: 1: simply supported, 2: clamped
st_bctype = 2;
outtype = 2;         % output type: 1: displacement, 2: velocity
losstype = 2;        % loss type: 1: independant, 2: dependant
itype = 2;           % type of input: 1: struck, 2: plucked
interpJ = 0;         % J Interpolation type 0: floor, 1: linear

% Scheme Run Options
plot_on = false;      % plot output
play_on = true;        % play out audio at the end
run = true;          % run scheme
skpErr = false;       % skip error checking for debugging

% Analysis
modeAn = false;        % run modal anlysis
energyAn = false;     % run energy analysis
benchtest = false;   % turn benchtesting on to view MATLAB profiler

%%%%%%%%% EDIT THESE %%%%%%%

%%START BENCH TEST%%
if benchtest
  profile on
end
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% EDIT THESE %%%%%%%

% simulation
Tf = 3;                      % duration
nu = .5;                     % Poisson Ratio (< .5)
rp = [.15 .15;.15 .85];      % readout position as percentage on grid.

% I/O
OSR = 1;                     % Oversampling ratio
SR = 44.1e3;                 % sample rate (Hz)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Physical Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Plate

% pl_E  = 2e11;                   % Young's modulus
% pl_rho = 7850;                  % density (kg/m^3)
pl_E  = 11e9;                     % Young's modulus
pl_rho = 480;                     % density (kg/m^3)

pl_H = .005;                    % thickness (m)
pl_L = 1;
pl_Lx = .5;                       % plate length (m)
pl_Ly = .5;                       % plate length (m)
pl_loss = [100, 1; 1000, .75];    % loss [freq.(Hz), T60;...]
pl_u0 = 0; pl_v0 = 1;           % excitation displacement and velocity
pl_wid = .2;                    % width (m or %)?
pl_ctr = [.35 .35; .71 .78];    % coordinates of string coupling [Xs1,Ys1; Xs2,Ys2];

%%%% String

st_gauge = 25;                  % string gauge
st_note = {'A4'};               % string note Pitch Notation 'N#' (N is the note # is the octave, must be a string in a cell array)
st_f0 = note2hz(st_note);       % frequency of note in Hz (see function at EOF)
st_r = (st_gauge * 2.54e-5)*.5; % string radius (m)
st_L = .6;                       % length (m)
st_E = 2e11;                    % Young's modulus (Pa) (GPa = 1e9 Pa) of steel
st_rho = 7850;                  % density (kg/m^3) of steel
st_loss = [100, 8; 1000 7];     % loss [freq.(Hz), T60;...]


%%%% % I/O string
st_xi = 0.75;                   % coordinate of excitation (normalised, 0-1)
st_famp = 10;                  % peak amplitude of excitation (N)
st_dur = 0.001;                  % duration of excitation (s)
st_exc_st = 0;                  % start time of excitation (s)
st_w0 = 1; st_v0 = 0;
st_ctr = .7; st_wid = .1;

st_T = (((2*st_f0.*st_r).*st_L).^2)*pi*st_rho; % Tension in Newtons

%%%% Extra Ploting Variables
pTr = 5e-6;
sTr = 5e-6;
FrameDrop = 10;
spin = 1e3;         % control spin of plotting
