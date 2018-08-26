%%%%% FDTD 2D Wave Model MIDI input
%%%%% Matthew Hamilton s0674653
%%%%% Description:
%%%%%
%%%%% Instrument file for Coupled String and Plate Reverb
%%%%% NOTE: Since the string is stiff, Length of all strings
%%%%% is the same to emphasise the stiffness


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Flags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Conditions
bctype = 2;          % boundary condition type: 1: simply supported, 2: clamped
st_bctype = 1;
outtype = 2;         % output type: 1: displacement, 2: velocity
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

filetypes = {'*.wav';'*.flac';'*.ogg';'*.mp3'; '*.m4a'; '*.mp4'};

[sfile, spath] =  uigetfile(filetypes,'Choose a Sound File'); %% No Message in OS X
if ~sfile; error('Load file cancelled'); end

[force, SR] = audioread([spath sfile]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% EDIT THESE %%%%%%%

% simulation
nu = .4;                     % Poisson Ratio (< .5)
rp = [.45, .85; .45, .15];              % readout position as percentage on grid.

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
pl_L = 1.4;                       % plate length (m)
pl_Lx = 1;                       % plate length X axis (m)
pl_Ly = .5;                       % plate length Y axis (m)

pl_loss = [100, 1; 1000, .5];    % loss [freq.(Hz), T60;...]
% pl_ctr = [.25 .25; .75 .75];    % coordinates of string coupling [Xs1,Ys1; Xs2,Ys2];

%%%% String

% MIDI note number to Hz
% tuning = [25 27 28 30 32 33 35 37];
% st_f0 = 440*2.^((tuning - 69)/12)';
tuning = {'A4' 'A4' 'A3'};               % string note Pitch Notation 'N#' (N is the note # is the octave, must be a string in a cell array)
st_f0 = note2hz(tuning)';       % frequency of note in Hz (see function at EOF)

strNum = length(st_f0); % number of strings
st_L = 1;           % string lengths
st_E = 2e11;        % Young's modulus (Pa) (GPa = 1e9 Pa) of steel
st_rho = 7850;      % density (kg/m^3) of steel
st_T = 400;         % Tension in Newtons
st_r = sqrt(st_T/(pi*st_rho))./(2*st_f0*st_L); % string radius (m)

st_loss = [100, 6; 1000 5];     % loss [freq.(Hz), T60;...]
%%%% % I/O string
st_xi = 0.75;                 % coordinate of excitation (normalised, 0-1)
st_xo = .2                    % coordinate of output (normalised, 0-1)
st_famp = 1;                  % peak amplitude of excitation (N)
st_dur = 0.001;               % duration of excitation (s)
st_exc_st = 0;                % start time of excitation (s)
st_w0 = 1; st_v0 = 0;
st_ctr = .7; st_wid = .1;
pl_ctr = linspace(.1,.9,strNum);  %% coupling points on string

%% panning
pan = linspace(-1, 1, length(st_f0));
panLR = [cos((pi*(pan+1)/4)); sin((pi*(pan+1)/4))];

Tf = length(force)/SR;

%%%% Plotting
pTr = 1e-5;
sTr = 2e-3;
FrameDrop = 10;

%%%% Reverb extras
% extra variables relating to the reverb nature of the script
zpad = pl_loss(2,2);               % pad the end of the input with some zeroes
Nf = length(force)+ zpad*SR;      % number of time steps
f_ctr = [.43 .43];                % cordinates of input force

%%%%% Oversampling
SR = OSR*SR;
f_chan = 1;
force = interp(force(:,f_chan),OSR);

%%%% Readout From string
%panning
pan = linspace(-1, 1, length(tuning));
panLR = [cos((pi*(pan+1)/4)); sin((pi*(pan+1)/4))];
