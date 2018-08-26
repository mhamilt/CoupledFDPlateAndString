%%%%% ASSIGNMENT 6: PMMI
%%%%% Matthew Hamilton s0674653
%%%%% Description:
%%%%%
%%%%% Finite difference string including panning, frequency dpendant loss,
%%%%% implicit computation, string overwinding simulation

clear all
% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Flags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% EDIT THESE %%%%%%%
itype = 1;                  % type of input: 1: struck, 2: plucked
bctype = 1;                  % boundary condition type: 1: simply supported, 2: clamped
outtype = 2;                % output type: 1: string displacement, 2: string velocity
losstype = 2;
plot_string = false;
plot_waves = true;          % plot output waveform and Force waveform
play = true;                 % play out audio at the end
benchtest = false;           % turn benchtesting on to view MATLAB profiler
%%%%%%%%% EDIT THESE %%%%%%%

%%START BENCH TEST%%
if benchtest
  profile on
end
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% EDIT THESE %%%%%%%
% guitar parameters
gauge = [20]; % list of string gauges
tuning = {'F3'}; % list of notes in scientific notation

%%% comparison
path = 'PianoSamples/';
file = '045.wav';
[y1, SR1] = audioread([path file]);

% physical string parameters
L = .44;                     % length (m)
E = 2e11;                  % Young's modulus (Pa) (GPa = 1e9 Pa) of steel
rhoP = 7850;                 % density of plain string Material
% TT = .895                     %%TensionTweak
loss = [100, 8; 1000 7];   % loss [freq.(Hz), T60;...]


% I/O
OSR = 1;                            % Oversampling ratio
SR = 44.1e3;                        % sample rate (Hz)
SR = SR*OSR;                        % redefine SR by OSR
xi = 0.8;                           % coordinate of excitation (normalised, 0-1)
famp = 1;                           % peak amplitude of excitation (N)
dur = 0.001;                       % duration of excitation (s)
xo = 0.45;                    % coordinate of output (normalised, 0-1)
window_dur = 0.02;          % duration of fade-out window (s). See Q6.
%%%%%%%%% EDIT THESE %%%%%%%


%%%%%%%%% DONT EDIT THESE %%%%%%%
r = (gauge * 2.54e-5)*.5;   % string core radius (m)
f0 = note2hz(tuning);     % Hz value of note based on function at EOF
strNum = length(f0);            % number of strings
rho = repmat(rhoP,1,strNum);       % density for each string: Deafult is rhoP
exc_st = ([1:strNum]*.5)-.5;        % start times of excitation (s)
Tf = max(loss(:,2)) + max(exc_st);    % duration of simulation (s)
T = pi*rho.*((f0*2*L).^2).*(r.^2);    % Tension Newtons
T = ((2*f0.*r.*L).^2)*pi*rho; % Tension in Newtons
%%%%%%%%% DONT EDIT THESE %%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Overwound Strings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Error Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
skpErr = false; % skip error checking for debugging

% Array of variable names and conditions they must meet.
paraLs = {'SR' 'r' 'xi' 'dur' 'exc_st' 'xo' 'Tf' 'loss' ,'gauge/tuning'};
paraCons = [(SR<1e3), any((r<=0)), (xi<0 || xi>1), (dur<(1/SR)),...
(any(sum([(exc_st>Tf) (exc_st<0)]))),(xo<0 || xo>1),...
(Tf < max(loss(:,2))+exc_st), any(loss(:,2)<=0), (L<=0), (length(gauge)~=length(tuning))];

if ~skpErr
  if any(paraCons)
    for paraError = paraLs(paraCons)
      switch char(paraError)

        %%% WARNING ERRORS %%%
      case 'Tf'
        warning([sprintf('Tf of %.2f secs is too short for a T60 of %.2f secs\n',...
        Tf, max(loss(:,2)))...
        sprintf('Don''t forget the excitation time!\nChanging Tf to equal T60 + exc_st\n')])
        Tf = max(loss(:,2))+exc_st;

      case 'r'
        if r<0
          warning([sprintf('''r'' is negative for some reason\n'),...
          sprintf('Changing r to absolute value\n')])
        else
          error([sprintf('r is equal to zero\n')])
        end

      case 'dur'
        warning([sprintf('Excitation duration of %.2f ms is too short\n',...
        dur*10^3)...
        sprintf('Changing dur to equal time step 1/SR\n')])
        dur = 1/SR;

        %%% FATAL ERRORS %%%
        % basically case {'SR' 'xi' 'exc_st' 'xo' 'rw' 'T60'}
        % if there isn't a custom error for a parameter
        % it will default to a fatal error
      otherwise
        error([sprintf('The following variables are invalid\n\n')...
        sprintf('''%s'' \n',paraLs{paraCons}) ...
        sprintf('\nCheck their values and try again\n')])
      end

    end
    %% Check N will be greater than 10000
    k = 1/SR;                   % time step
    Nlim = L./sqrt( (.5*(k^2)./(pi*rho.*(r.^2))).*(T+sqrt( ( ( T.^2 + ((4*E.*rho.*(r.^3*pi*SR).^2 )))))));
    if Nlim > 10e3
      error('N too high, change the string length or make T or E bigger')
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%panning
pan = linspace(-1, 1, length(f0));
panLR = [cos((pi*(pan+1)/4)); sin((pi*(pan+1)/4))];


A = pi*r.^2;                 % string cross-sectional area of whole string
I = 0.25*pi*r.^4;            % string moment of intertia
% T = T * (1 - E*I/T)*TT;            %% compensate tuning


c = sqrt(T./(rho.*A));        % wave speed based on internal core
K = sqrt(E*I./(rho.*A));      % stiffness constant based on whole string


%%% loss coefficients


if losstype ==0
  % frequency independant loss
  sig0 = 0;
  sig1 = 0;
end

if losstype ==1
  % frequency independant loss
  sig0 = (6*log(10)/loss(1,2));
  sig1 = 0;

end

if losstype == 2
  z1 = (-c.^2 + sqrt(c.^4 + 4*K.^2.*(2*pi*loss(1,1))^2))./(2*K.^2);
  z2 = (-c.^2 + sqrt(c.^4 + 4*K.^2.*(2*pi*loss(2,1))^2))./(2*K.^2);
  sig0 = 6*log(10)*(-z2/loss(1,2) + z1/loss(2,2))./(z1-z2);
  sig1 = 6*log(10)*(1/loss(1,2) - 1/loss(2,2))./( z1-z2);
end

%%%%% grid

k = 1/SR;                   % time step
hmin = sqrt(0.5* (c.^2*k^2+sqrt(c.^4*k^4+16*K.^2.*k.^2)) );    % minimal grid spacing

N = floor(L./hmin);          % number of segments (N+1 is number of grid points)
h = L./N;                    % adjusted grid spacing

lambda = c*k./h;             % Courant number
mu = K*k./h.^2;               % numerical stiffness constant

N = N + 1;

%%%%% I/O

Nf = floor(SR*Tf);          % number of time steps

li = floor(xi*N);         % grid index of excitation
lo = floor(xo*N);         % grid index of output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check that xo and xi are greater than h distance away from boundaries
if any([(xi*L)<h (L-(xi*L))<h])
  warning([sprintf('xo is too close to string boundary\n'),...
  sprintf('changing xo to equal h\n')]);
  xi = abs(round(xi)-h);
  if any(li<=2)
    li(li<=2) = 3;
  end
  if any(li>=N)
    li(li>=N) = N(li>=N)-1;
  end
end

if any([(xo*L)<h (L-(xo*L))<h])
  warning([sprintf('xo is too close to string boundary\n'),...
  sprintf('changing xo to equal h\n')]);
  xo = abs(round(xo)-h);
  if any(lo<=2)
    lo(lo<=2) = 3;
  end
  if any(lo>=N)
    lo(lo>=N) = N(lo>=N)-1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Force Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = sparse(Nf,strNum);                   % input Force signal
durint = floor(dur*SR);                 % duration of Force signal, in samples
exc_st_int = floor(exc_st*SR);          % start time index for excitation
exc_end_int = exc_st_int+durint-1;
d0 = (k^2./(h.*rho.*A*(1+ k*sig0)));               % Input coefficient

for en = 1:strNum
  n = (exc_st_int(en):exc_end_int(en))+1;
  f(n,en)=famp * 0.5 * ...
  (1-cos((2/itype)*pi.*(n'-exc_st_int(en))/durint))*d0(en);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% initialise output
y = zeros(Nf,2);

%%%% Initialise Master Matrices
AA=cell(1,strNum);
BB=AA;
CC=AA;

%%%% Grid space count
cN = zeros(1,strNum+1);

tic

%%%%% Scheme Matrices
Dxx = fidimat(N,'xx');
Dxxxx = fidimat(N,'xxxx',bctype);
% sI = fidimat(N,'I');
sI = speye(N);

%%%%% Coefficient Matrices
AA = 1/(1+k*sig0);
BB = ((lambda^2 * Dxx) + (2*sI) - (mu^2*Dxxxx) + ((2*sig1*k/h^2)*Dxx)) *AA;
CC = -((1-k*sig0)*sI + ((2*sig1*k/h^2)*Dxx)) * AA;
BB([1 end],:)=0; CC([1 end],:)=0;




%%%%% initialise scheme variables
u2 = zeros(N,1);              % state
u1 = u2;                        % state
u = u2;                         % state

%%%% input
fvect = u2;                     % vector for including f(n,s)
sumstr = ones(1,strNum);        % sum strings together


%%% redefine read out points from offset count cN
li = li;
lo = lo;

%%%%% Main loop
for n=1:Nf

  % input update
  fvect(li) = f(n,:);

  u = BB*u1 + CC*u2 + fvect;

  if plot_string
    plot(u)
    ylim([-2e-3 2e-3])
    drawnow
  end

  if (outtype==1)
    y(n,:) = panLR*u(lo);
  elseif (outtype==2)
    y(n,:) = panLR*(SR*(u(lo)-u1(lo)));
  end

  % shift state
  u2 = u1; u1 = u;
end

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
winL = floor(window_dur*SR); % window length in samples
fade = 0.5*(1-cos(pi.*([1:winL]/winL)+ pi))';
y(end-winL+1:end,:) = y(end-winL+1:end,:).*repmat(fade,1,2); %Fade out audio at end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% play sound
if play
  y = y/max(abs(y(:)));
  soundsc(y,SR);
end

%%%%% plot spectrum
if plot_waves
  figure(1)
  y = (y/abs(max(y(:))));
  y1 = (y1/abs(max(y1(:))));

  yfft = 10*log10(abs(fft(y(:,2))));
  plot(([0:Nf-1]'/Nf)*SR, yfft, 'k')
  fft_ax = gca;
  % xlim([f0*.9 f0*1.1])
  xlim([0 f0*10])
  hold on
  line([f0, f0], fft_ax.YLim,'Color',[0 0 0],'LineWidth',2)
  plot([0:length(y1)-1]*SR1/length(y1),10*log10(abs(fft(y1(:,1)))),'g')
  legend('FFT Output_{dB}','F0_{Hz}')
  set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen
end

if benchtest
  profile viewer
  profile off
end
