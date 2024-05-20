%%%%% FDTD 2D Wave Model MIDI input
%%%%% Matthew Hamilton s0674653
%%%%% Description:
%%%%%
%%%%% Instrument file for 2D Plate
%%%%% Contains all the parameters to run the FDTD thin plate scheme

clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Flags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../functions/');
addpath('../../midi/');
addpath('../../calibration/');
addpath('../../schemes/');
% Conditions
bctype = 2;          % boundary condition type: 1: simply supported, 2: clamped
st_bctype = 1;
outtype = 1;         % output type: 1: displacement, 2: velocity
losstype = 2;        % loss type: 1: independant, 2: dependant
itype = 2;           % type of input: 1: struck, 2: plucked
stiff = true;

% Scheme Run Options
plot_on = false;      % plot output
play = false;        % play out audio at the end
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

% -----------------------READ MIDI DATA----------------------------- %
path = '../../midi/';
file = 'mood.mid';
filename = [path file];
BPM = 156;
MIDI = read_midi_file(filename,BPM);
jiggle = 3000; % make plucks out of time, in samples
% ---------------------------------------------------- %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% EDIT THESE %%%%%%%

% simulation
nu = .5;                     % Poisson Ratio (< .5)
rp = [.15, .85; .85, .15];              % readout position as percentage on grid.

% I/O
OSR = 1;                     % Oversampling ratio
SR = 44.1e3;                 % sample rate (Hz)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Physical Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plate

pl_E  = 2e11;                   % Young's modulus
pl_rho = 7850;                  % density (kg/m^3)
% pl_E  = 11e9;                   % Young's modulus
% pl_rho = 480;                  % density (kg/m^3)

pl_H = .005;                    % thickness (m)
pl_L = 1.4;                       % plate length (m)
pl_Lx = 1;                       % plate length X axis (m)
pl_Ly = .5;                       % plate length Y axis (m)

pl_loss = [100, 3; 1000, 2.5];    % loss [freq.(Hz), T60;...]
% pl_ctr = [.25 .25; .75 .75];    % coordinates of string coupling [Xs1,Ys1; Xs2,Ys2];

%%%% String

tuning = unique(MIDI(:,1)); % number and value of unique notes in MIDI sequence
st_f0 = 440*2.^((tuning - 69)/12)'; % MIDI note number to Hz
strNum = length(st_f0); % number of strings

st_gauge = 20;                  % string gauge
st_T = 850; % Tension in Newtons
st_r = (st_gauge * 2.54e-5)*.5; % string radius (m)
st_E = 2e11;                    % Young's modulus (Pa) (GPa = 1e9 Pa) of steel
st_rho = 7850;                  % density (kg/m^3) of steel

%% Nylon
st_E = 2e11;
st_rho = 7850;

st_loss = [100, 6; 1000 5];     % loss [freq.(Hz), T60;...]
st_L = sqrt(st_T./(pi*st_rho))./(2*st_f0.*st_r); % string lengths

%%%% % I/O string
st_xi = 0.75;                   % coordinate of excitation (normalised, 0-1)
st_famp = 1;                  % peak amplitude of excitation (N)
st_dur = 0.001;                  % duration of excitation (s)
st_exc_st = 0;                  % start time of excitation (s)
st_w0 = 1; st_v0 = 0;
st_ctr = .7; st_wid = .1;

exc_st = [];
for n = 1:length(tuning)
  times = MIDI(find(MIDI(:,1) == tuning(n)),3);        % start times of excitation (s)
  times(:,2) = n;
  exc_st = [exc_st; times];
end

Tf = max(st_loss(:,2)) + max(exc_st(:,1));             % duration of simulation (s)


pl_ctr = linspace(.1,.9,strNum);
pTr = 1e-5;
sTr = 2e-3;
FrameDrop = 10;
