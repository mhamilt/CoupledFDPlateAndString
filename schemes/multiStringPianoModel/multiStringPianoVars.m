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
itype = 1;           % type of input: 1: struck, 2: plucked


% Scheme Run Options
plot_on = false;      % plot output
play = false;        % play out audio at the end
run = true;          % run scheme
skpErr = false;       % skip error checking for debugging
interpJ = true;      %

% Analysis
modeAn = false;        % run modal anlysis
energyAn = false;     % run energy analysis
benchtest = false;   % turn benchtesting on to view MATLAB profiler

%%%%%%%%% EDIT THESE %%%%%%%

% -----------------------READ MIDI DATA----------------------------- %
path = 'midi/';
file = 'thirdman.mid';
filename = [path file];
BPM = 156;
MIDI = read_midi_file(filename,BPM);
% MIDI(:,1)= MIDI(:,1)-20;    %%%NOTE THIS NEEDS REVERSED CORRECTLY%%%
MIDI(:,1)= 89-(MIDI(:,1)-20);

jiggle = 20; % make plucks out of time, in samples
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

%%%% Steel
% pl_E  = 2e11;                   % Young's modulus
% pl_rho = 7850;                  % density (kg/m^3)

%%%% Generic wood
pl_E  = 11e9;                   % Young's modulus
pl_rho = 480;                    % density (kg/m^3)

pl_H = .005;                    % thickness (m)
pl_Ly = 2.88;                       % plate length Y axis (m)
pl_Lx = 1;                       % plate length X axis (m)

pl_loss = [100, 1.3; 1000, 1];    % loss [freq.(Hz), T60;...]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% String
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Parameters derived from calibrating strings.
multiStringGrandPiano

st_f0{1} = note2hz(tuning);                            % Hz value of note based on function
st_r{1} = (st_gauge * 2.54e-5)*.5;                  % string core radius (m)
st_T{1} = ((2*st_f0{1}.*st_r{1}.*st_L).^2)*pi.*st_rho; % Tension in Newtons
% st_loss = [100, 8, 1000 7]';                      % loss [freq.(Hz), T60;...]
strNum = 88;                                     % number of strings
Ts = {88,75,62};                % number of string on each tier
fStrNum = 88 +75 +62;           % full number of strings
To = {0,88,163};                % offset for each tier

%%%% % I/O string
st_xi = 0.8;                   % coordinate of excitation (normalised, 0-1)
st_famp = 1;                  % peak amplitude of excitation (N)
st_dur = 0.0018;                  % duration of excitation (s)
st_dur = logspace(log10(.5e-3),log10(2.25e-3),88);
% st_dur = linspace(1e-3,2.5e-3,88);
st_exc_st = 0;                  % start time of excitation (s)
st_w0 = 1; st_v0 = 0;
st_ctr = .7; st_wid = .1;


TT = [1,.98,1.01];
%%% tiered strings
for t = 2:3

  %%% range of string tier
  Tr = 1:Ts{t};
  st_f0{t} = st_f0{1}(Tr);
  st_T{t} = st_T{1}(Tr)*TT(t);
  st_r{t} = st_r{1}(Tr);

end



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
% Tf = max(st_loss(:,2)) + max(exc_st(:,1));  % duration of simulation (s)
Tf = 2 + max(exc_st(:,1));  % duration of simulation (s)
pl_ctr = linspace(.1,.9,strNum);            % connection point of strings
pl_rx = {.2 .21 .22};                                 % depth of connection into x-axis

%%%% Visualising
pTr = 5e-4;           %%% Z limit for plate
sTr = 2e-3;           %%% y limit for strings.
FrameDrop = 100;    %%% number of frames to skip
