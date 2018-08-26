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
losstype = 1;        % loss type: 1: independant, 2: dependant

% Scheme Run Options
plot_on = false;     % plot output
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

%% input file

% [sfile, spath] =  uigetfile({'*.wav';'*.flac';'*.ogg';...
%                           '*.mp3'; '*.m4a'; '*.mp4'},... % unsupported in OS X
%                           'Choose a Sound File'); %% No Message in OS X
% if ~sfile; error('Load file cancelled'); end
%
% [force, SR] = audioread([spath sfile]);


% simulation
Tf = 1;                      % duration
nu = .5;                     % Poisson Ratio (< .5)
ctr = [.5, .5];              % centre point of excitation as percentage
wid = .1;                    % width (m or %)?
u0 = 0; v0 = 1;              % excitation displacement and velocity
% famp = 1;                  % excitation force amplitude
rp = [.45 .35];              % readout position as percentage on grid.
% mu = .25                   % free parameter


% physical parameters
E = 2e11;                    % Young's modulus
rho = 7850;                  % density (kg/m^3)

% wood
% E = 11e9;
% rho = 4.5;

H = .01;                     % thickness (m)
L = 1;                       % plate length (m)


% T60 = 3;
loss = [100, 8; 1000, 6];    % loss [freq.(Hz), T60;...]

% I/O
OSR = 1;                     % Oversampling ratio
SR = 44.1e3;                 % sample rate (Hz)

% Force
fHz = 258.95*2;                   % frequency of force input


% plotting
pTr = 1e-6;                   % max trans limit of of plot
