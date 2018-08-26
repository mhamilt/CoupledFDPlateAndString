%%%%% FDTD 2D Wave Model
%%%%% Matthew Hamilton s0674653
%%%%% Description:
%%%%%
%%%%% A Kirchhoff Thin Plate FDTD Model

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Instrument File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inst % file containing all relevent parameters for this set of schemes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Motion Coefficients
D = (E*(H)^3)/(12*(1-(nu^2)));
kappa = sqrt(D / (rho*  H) );

%%%%% Scheme Spacing
k = 1/SR;                   % time step
hmin = 2*sqrt(k*kappa);      % NSS (equation 12.5)
N = floor(L./hmin);          % number of segments
h = L./(N);                    % adjusted grid spacing
mu = (kappa * k)/(h^2);      % scheme parameter

N = N+1;                     % for includng only ALL grid points %%NOTE%%
ss = N*N;                    % total grid size.
%%%%% I/O
Nf = Tf*SR;          % number of time steps


% force = (sin(2*pi*(0:Nf-1)*k*fmode))*2;
force = [(sin(2*pi*(0:22e3-1)*k*fmode))*2 zeros(1,Nf-22e3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % Loss coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if losstype ==0
  % frequency independant loss
  sigma0 = 0;
  sigma1 = 0;
end

if losstype ==1
  % frequency independant loss
  sigma0 = 6*log(10)/loss(1,2);
  sigma1 = 0;
end

if losstype == 2
  %% this is simply loss from NSS with 'c' removed as tension is not present.
  z1 = 2*kappa*(2*pi*loss(1,1))/(2*kappa.^2); %from 1D
  z2 = 2*kappa*(2*pi*loss(2,1))/(2*kappa.^2);

  sigma0 = 6*log(10)*(-z2/loss(1,2) + z1/loss(2,2))./(z1-z2);
  sigma1 = 6*log(10)*(1/loss(1,2) - 1/loss(2,2))./(z1-z2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % Read In/Out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lo = rp*N;
lo = floor(sub2ind([N N],lo(1), lo(2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Force Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d0 = (k^2)/(rho*H*(h^2))*(1/(1+k*sigma0)); % input force coefficient
% force = [force;zeros(zpad*SR,1)]*d0;       % input oscillator
force = force*d0;       % input oscillator
fIn = sub2ind([N N],floor(ctr(1)*N),floor(ctr(2)*N));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficient Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BH = biharm2(N,N,bctype); % biharmonic matrix
LA = laplace2(N,N,bctype); % Laplacian matrix

A = (1/(1+k*sigma0))*speye(ss);   %NOTE% Currently inverted
B = (-(mu^2)*BH + (2*sigma1*k/(h^2))*LA + 2*speye(ss)) * A;
C = (-(2*sigma1*k/(h^2))*LA - (1-sigma0*k)*speye(ss))  * A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialise I/O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% initialise output
u2 = zeros(N*N,1);
u1 = u2;
u  = u2;
Fv = u2;              % force vector
%%% initialise scheme variables


%%%% input


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_on
  fig = figure(1);
  fig.Color = [0 0 0];
  ax = axes('XLim',[0 N],'YLim',[0 N],'ZLim',[-pTr pTr]);
  % colormap copper
  % shading interp
  material SHINY
  ax.Projection = 'perspective';

  set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen
  view(15,15)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Energy Analysis
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

%%%%%
%%%%% Main loop
tic
if run
  for n = 1:Nf

    % update input forces
    Fv(fIn) = force(n);

    % main operation
    u = B*u1 + C*u2 + Fv;

    % plotting
    if plot_on
      if ~mod(n,FrameDrop)

        %%% plate
        surf(ax,reshape(u,N,N),'FaceLighting','gouraud','FaceColor','interp',...
        'AmbientStrength',0.5)
        light('Position',[1 0 1],'Style','local')
        colormap copper
        ax.Projection = 'perspective';
        % zlim([-pTr pTr])
        % caxis([-1e-7 3e-6]);
        % shading interp
        % view(45+spin*n/(Nf*FrameDrop),30)
        % view(2)
        view(45,45)
        ax.XLim =[0 N]; ax.YLim =[0 N]; ax.ZLim =[-pTr pTr];
        axis off
        F(n/FrameDrop) = getframe(gca);
        drawnow

      end
    end

    % read output
    % if (outtype==1)
    %   y(n,:) = u(lo);
    %
    % elseif (outtype==2)
    %   y(n,:) = (SR*(u(lo)-u1(lo)));
    %
    % end

    % shift state
    u2 = u1; u1 = u;

    if energyAn % shift energy state
      KE(n) = coE(1)*sum((u1-u2).^2); PE(n) = coE(2)*((LA*u1)' * (LA*u2));
      Energy(n) = coE*[sum((u1-u2).^2); (LA*u1)' * (LA*u2)];
    end

    if any(isnan(u))
      error('UNSTABLE script terminated')
    end


  end
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (Housekeeping)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Modal Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculating plate modes

p = [1:5];
q = [1:5]';

[P,Q] = meshgrid(p,q);

mfreqs = (pi*kappa/(2*L^2))*((P.^2) + (Q.^2));     % mode frequencies


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
font_s = 14; % font point size
set(0,'DefaultAxesFontSize',font_s);
YF = abs(fft(y));

% Output FFT
fy = linspace(0, SR, length(y)); % frequency axis

figure(1);
plot(fy, 20*log10(YF));
ax1 = gca;
line([mfreqs(:) mfreqs(:)],ax1.YLim,'Color',[1 0 0])
ax1.XLim = [0, 1500];
title('Scheme Output FFT', 'FontSize', font_s)
xlabel('Freq._{Hz}');
ylabel('magnitude');
legend('FFT Output_{dB}','Mode Frequencies_{Hz}')
set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if benchtest
  profile viewer
  profile off
end

% EOF
