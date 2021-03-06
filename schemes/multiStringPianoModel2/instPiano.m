%%%%% FDTD 2D Wave Model MIDI input
%%%%% Matthew Hamilton s0674653
%%%%% Description:
%%%%%
%%%%% Instrument file for 2D Plate
%%%%% Contains all the parameters to run the FDTD thin plate scheme


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Flags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Conditions
bctype = 2;          % boundary condition type: 1: simply supported, 2: clamped
st_bctype = 1;
outtype = 2;         % output type: 1: displacement, 2: velocity
losstype = 2;        % loss type: 1: independant, 2: dependant
itype = 2;           % type of input: 1: struck, 2: plucked

% Scheme Run Options
plot_on = true;      % plot output
play = false;        % play out audio at the end
run = true;          % run scheme
skpErr = false;       % skip error checking for debugging
interpJ = true;       %% J Interpolation type 0: floor, 1: linear

% Analysis
modeAn = false;        % run modal anlysis
energyAn = false;     % run energy analysis
benchtest = false;   % turn benchtesting on to view MATLAB profiler

%%%%%%%%% EDIT THESE %%%%%%%

% -----------------------READ MIDI DATA----------------------------- %
path = 'midi/';
file = 'chord.mid';
filename = [path file];
BPM = 156;
MIDI = read_midi_file(filename,BPM);
% MIDI extends the piano range, -20 aligns it with the
% values on a piano.
MIDI(:,1)= MIDI(:,1)-20;
% make plucks out of time, in samples
jiggle = 1000;
% ---------------------------------------------------- %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% EDIT THESE %%%%%%%

% simulation
nu = .4;                     % Poisson Ratio (< .5)
rp = [.85, .15; .85, .75];   % readout position as percentage on grid.

% I/O
OSR = 1;                     % Oversampling ratio
SR = 44.1e3;                 % sample rate (Hz)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Physical Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plate

% pl_E  = 2e11;                   % Young's modulus
% pl_rho = 7850;                  % density (kg/m^3)
pl_E  = 11e9;                   % Young's modulus
pl_rho = 480;                  % density (kg/m^3)

pl_H = .005;                    % thickness (m)
pl_Ly = 1.8;                       % plate length Y axis (m)
pl_Lx = 1.2;                       % plate length X axis (m)

pl_loss = [100 .4; 1000 .2];
% pl_ctr = [.25 .25; .75 .75];    % coordinates of string coupling [Xs1,Ys1; Xs2,Ys2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% String
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Parameters derived from calibrating strings.
GrandPiano

st_f0 = note2hz(tuning);                            % Hz value of note based on function
st_r = (st_gauge * 2.54e-5)*.5;                  % string core radius (m)
st_T = ((2*st_f0.*st_r.*st_L).^2)*pi.*st_rho; % Tension in Newtons
% st_loss = [100, 8, 1000 7]';                      % loss [freq.(Hz), T60;...]
strNum = 88;                                     % number of strings

%%%% % I/O string
st_xi = 0.85;                   % coordinate of excitation (normalised, 0-1)
st_famp = .1;                  % peak amplitude of excitation (N)
st_dur = 0.0018;                  % duration of excitation (s)
st_dur = logspace(log10(1.5e-3),log10(.5e-3),88);
% st_dur = linspace(2.5e-3,1e-3,88);
st_exc_st = 0;                  % start time of excitation (s)
% st_w0 = 1; st_v0 = 0;
st_ctr = .7; st_wid = .1;

%%%% Parse MIDI notes
exc_st = [];
notes = unique(MIDI(:,1));
note_i = 1:length(notes);         % note index for force vector vector.
for n = 1:length(notes)
  times = MIDI(find(MIDI(:,1) == notes(n)),3);        % start times of excitation (s)
  times(:,2) = n;
  times(:,3) = MIDI(find(MIDI(:,1) == notes(n)),2);    % note velocities
  exc_st = [exc_st; times];
end

%%%% Extra Details
Tf = 8 + max(exc_st(:,1));            % duration of simulation in seconds
pl_ctr = linspace(.1,.9,strNum);      % connection point of strings
pl_rx = .2;                           % depth of connection into x-axis

%%%% Visualising
pTr = 1e-4;           %%% Z limit for plate
sTr = 6e-5;           %%% y limit for strings.
FrameDrop = 2000;     %%% number of frames to skip
spin = 3e5;         % control spin of plotting
