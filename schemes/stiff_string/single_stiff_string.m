clearvars, 
close all

%%%%% flags

plot_on = 0;                % in-loop plotting on (1) or off (0)
itype = 1;                  % type of input: 1: struck, 2: plucked
bctype = 1;                 % boundary condition type: 1: simply supported, 2: clamped
outtype = 1;                % output type: 1: string displacement, 2: string velocity

%%%%% parameters

% physical string parameters (see Q7a)
L(1:6) = 0.648;                 % guitar scale length(m)
rho(1:6)= 7750;                 % Density of steel (kg/m^3)
% Strings:     E        B        G       D       A     E
% gauges:   .046     .036     .026    .017    .013   .01
r   = [1.1684e-3 9.144e-4 6.604e-4 4.3e-4 3.3e-4 2.5e-4];   % string radius (m)

% Young's modulus (Pa) taken from
% http://srjcstaff.santarosa.edu/~yataiiya/E45/PROJECTS/guitar%20strings.pdf
% considering the lowest three strings to be wound and the top three
% single-wire
E   = [6.4e10 6.4e10 6.4e10 1.88e11 1.88e11 1.88e11] ;
freq= [82.41 110 146.83 196 246.94 329.63]; % Note frequencies



% T60 (s) was left to 5 s. Although longer T60 were found by analysing
% decays of plucked guitar strings, 5 s produced nicer sound. 
T60(1:6) = 5;

% I/O
SR = 44100;                 % sample rate (Hz)
Tf = 5;                     % duration of simulation (s)

xi(1:6) = 0.3;              % coordinate of excitation (normalised, 0-1)
famp = [1 1 1 1 1 1];                        % peak amplitude of excitation (N)
dur = [0.001 0.001 0.001 0.001 0.001 0.001]; % duration of excitation (s)
exc_st = [0 0.5 1 1.5 2 2.5];                % start time of excitation (s)
xo(1:6) = 0.2;              % coordinate of output (normalised, 0-1)
window_dur = 0.01;          % duration of fade-out window (s). See Q6.

%%%%% Error checking

if any([r, rho, T60, L, SR, Tf, famp, window_dur]<=0)
    error('Parameters must be positive');
end
if any(E<0)
    error('Young''s modulus must be non-negative');
end
if any(exc_st<0 | exc_st>Tf)
    error('Start time of excitation are beyond the simulation time')
end
if any([xi, xo]>=1 | [xi, xo]<=0)
    error('Excitation and output coordinates must be between 0 and 1')
end
if any([all(plot_on~=[0 1]) all(itype~=[1 2])...
        all(bctype~=[1 2]) all(outtype~=[1 2])])
    error('Flag values are incorrect')
end
if any(exc_st+dur>Tf)
    error('Excitation duration exceeds simulation duration')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% derived parameters
A = pi*r.^2;                    % string cross-sectional area
I = 0.25*pi*r.^4;               % string moment of intertia
K = sqrt(E.*I./(rho.*A));       % stiffness constant

% Tension (N) using stiff string frequency equation (7.21) from NSS
T = (4*L.^2.*freq.^2+K.^2*pi^2).*rho.*A;

c = sqrt(T./(rho.*A));        % wave speed

sig = 6*log(10)./T60;         % loss parameter

%%%%% grid

k = 1/SR;                   % time step
hmin = sqrt(0.5*(c.^2.*k^2+sqrt(c.^4.*k^4+16*K.^2.*k^2)));    % minimal grid spacing

N = floor(L./hmin);          % number of segments (N+1 is number of grid points)

%%%%% rest of error checking
if any(N>10000)
    error('Too many grid points.')
end

h = L./N;                    % adjusted grid spacing

lambda = c*k./h;             % Courant number
mu = K*k./h.^2;               % numerical stiffness constant

%%%%% I/O

Nf = floor(SR*Tf);          % number of time steps

li = floor(xi.*N);         % grid index of excitation
lo = floor(xo.*N);         % grid index of output

Nstr = length(T);            % number of strings


%%%% Adjusting li and lo in case they are out of range and calculating
%%%% indeces for I/O in the case of concatenated u

mli=zeros(1,Nstr);      % Input for matrix form
mlo=zeros(1,Nstr);      % Output for matrix form

for s=1:Nstr
    if any([li(s),lo(s)]<2 | [li(s), lo(s)]>N(s))
        warning('Param:outRange',...
            ['Excitation and/or output is outside the range of h to L-h.' ...
            '\n Forcing it to the range.'])
    end
    if li(s)<2
        li(s)=2;
    elseif li(s)>N(s)
        li(s)=N(s);
    end
    if lo(s)<2
        lo(s)=2;
    elseif lo(s)>N(s)
        lo(s)=N(s);
    end
    
    % Adding the summed number of grid points of the previous strings and
    % subtracting s number of points, since the 1st point of each string is
    % missing
    mli(s)=li(s)+sum(N(1:s-1))-s;
    mlo(s)=lo(s)+sum(N(1:s-1))-s;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create force signal

f = zeros(Nf,Nstr);                     % input force signal
durint = floor(dur*SR);                 % duration of force signal, in samples
exc_st_int = floor(exc_st*SR)+1;        % start time index for excitation


% Computing forcing terms for each string
for s=1:Nstr
    f(exc_st_int(s):exc_st_int(s)+durint(s)-1,s)=famp(s)*0.5.*...
        (1-cos(2/itype*pi*(0:durint(s)-1)/durint(s)));  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% scheme coefficients
out = zeros(Nf, 1);                % output
tic
% loop over the strings
AA=cell(1,Nstr);
BB=cell(1,Nstr);
CC=cell(1,Nstr);

% Finding matrices for each string
for s=1:Nstr
    % interior
    I = speye(N(s) - 1);
    Dxx = fidimat(N(s)-1,'xx', 1);
    Dxxxx = fidimat(N(s)-1,'xxxx', 1);
        
    a = 1/(1+sig(s)*k);
    BB{s} = a*(2*I + (lambda(s)^2*Dxx) + (-mu(s)^2) * Dxxxx);
    CC{s} = a*(-(1-sig(s)*k)*I);         
end

% Concatenating them into block diagonal matrices
AA=blkdiag(AA{:});
BB=blkdiag(BB{:});
CC=blkdiag(CC{:});


% input
d0 = 1./(1+sig*k).*(k^2./(h.*rho.*A));

%%%%% initialise scheme variables

u2 = zeros(sum(N)-Nstr,1);      % state
u1 = u2;                        % state
u = u2;                         % state

% plot
if(plot_on==1)
    % draw current state
    fig=figure;
    for s=1:Nstr
        subplot(Nstr,1,s)
        pl(s)=plot((0:N(s)-2)'*h(s), u(1-s+1+sum(N(1:s-1)):sum(N(1:s))-s), 'k');
        axis([0 L(s)-h(s) -0.0005 0.0005])
        % set data source (to be later used for refreshing data)
        set(pl(s),'YDataSource',['u(' num2str(1-s+1+sum(N(1:s-1)))...
            ':' num2str(sum(N(1:s))-s) ')'])
        ylabel(['String ' num2str(s)])
    end
    drawnow
end
    
%%%%% main loop
for n=1:Nf
    
    u=(BB*u1+CC*u2);
    
    % send in input
    u(mli) = u(mli)+(d0.*f(n,:))';
    
    % read output
    if(outtype==1)
        out(n)= sum(u(mlo));
    end
    if(outtype==2)
        out(n) = sum((u(mlo)-u1(mlo))/k);
    end
    
    % plot
    if(plot_on==1)
        refreshdata(fig)
        drawnow
    end
    
    % shift state
    u2 = u1;
    u1 = u;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% windowing the end
nWin=floor(SR*window_dur);
out(end-nWin:end)=0.5*(1+cos(pi*(0:nWin)/nWin))'.*out(end-nWin:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc
%%%%% play sound
soundsc(out,SR);

%%%%% plot spectrum

figure(1)
outfft = 10*log10(abs(fft(out)));
plot((0:SR/Nf:SR*(1-1/Nf))', outfft, 'k')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Frequency Spectrum of Output Signal')
set(gca, 'FontSize', 12)
xlim([0 SR/2])

