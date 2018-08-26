%%%%% FDTD 2D Wave Model
%%%%% Matthew Hamilton s0674653
%%%%% Description:
%%%%%
%%%%% A Kirchhoff Thin Plate FDTD Model with loss

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Instrument File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inst % file containing all relevent parameters for this set of schemes


%%%%%%%%% %%%%%%%%%% %%%%%%%


%%%%%%%%% DONT EDIT THESE %%%%%%%
SR = SR*OSR;                        % redefine SR by OSR

%%%%%%%%% %%%%%%%%%%%%%%% %%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Subsection Or Code Extension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Error Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
skpErr = true; % skip error checking for debugging

% Array of variable names and conditions they must meet.
paraLs = {'SR'};
paraCons = [(SR<1e3),];

if ~skpErr

  if any(paraCons)

    for paraError = paraLs(paraCons)

      switch char(paraError)

        %%% WARNING ERRORS %%%
      case 'string'
        % case

        %%% FATAL ERRORS %%%

      otherwise
        error([sprintf('The following variables are invalid\n\n')...
        sprintf('''%s'' \n',paraLs{paraCons}) ...
        sprintf('\nCheck their values and try again\n')])

      end

    end

  end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Motion Coefficients
% kappa = sqrt( (E*(H)^2) / (12*rho*(L^4)*(1-(nu^2))) );
Dx = (Ex*(H)^3)/(12*(1-((nu_xy*nu_yx))));
Dy = (Ey*(H)^3)/( 12* (1- (nu_xy*nu_yx) ));
Dxy = nu_yx*Dx + nu_xy*Dy + (G_xy*H^3)/3;

kappa_x = sqrt(Dx / (rho *  H ) );
kappa_y = sqrt(Dy / (rho *  H ) );
kappa_xy = sqrt(Dxy/(rho *  H ) );



%%%%% Scheme Spacing
k = 1/SR;                   % time step

hminx = 2*sqrt(k*kappa_x);      % NSS (equation 12.5)
Nx = floor(L./hminx);          % number of segments
hx = L./(Nx);                    % adjusted grid spacing (only internal points)

hminy = 2*sqrt(k*kappa_y);      % NSS (equation 12.5)
Ny = floor(L./hminy);          % number of segments
hy = L./(Ny);                    % adjusted grid spacing (only internal points)

% scheme parameter
mu_x = (kappa_x * k)/(hx^2);
mu_y = (kappa_y * k)/(hy^2);
mu_xy = (kappa_xy * k)/(hy*hx);

Nf = floor(SR*Tf);          % number of time steps

Nx = Nx+1;
Ny = Ny+1;
ss = Nx*Ny;                    % total grid size.
%%%%% I/O
Nf = floor(SR*Tf);          % number of time steps

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
  % sigma0 = 0
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
lo = rp*Nx;
lo = floor(sub2ind([Ny Nx],lo(1), lo(2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Force Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y] = meshgrid([1:Nx]*hx, [1:Ny]*hy);         % Grid of point in value of meters
dist = sqrt((X-(ctr(1)*L)).^2 + (Y-(ctr(2)*L)).^2); % distance of points from excitation
ind = sign(max(-dist+(wid*0.5),0));         % displaced grid points (logical)
rc = .5*ind.*(1+cos(2*pi*dist/wid));        % displacement
rc = rc(:);                                 % 2D plane as vector


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficient Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

XXXX = fidimat(Ny,Nx,'xxxx', 1); % biharmonic matrix
YYYY = fidimat(Ny,Nx,'yyyy', 1); % biharmonic matrix
XXYY = fidimat(Ny,Nx,'xxyy', 1);  % gradient matrix
LA = fidimat(Ny,Nx,'laplace', 1);  % gradient matrix
spI  = fidimat(Ny,Nx,'I', 1);

A = (1/(1+k*sigma0))*speye(ss);   %NOTE% Currently inverted
B = (-(mu_x^2)*XXXX - (mu_y^2)*YYYY - (mu_xy^2)*XXYY + (2*sigma1*k/(hx*hy))*LA + 2*spI) * A;
C = (-(2*sigma1*k/(hx*hy))*LA - (1-sigma0*k)*speye(ss))  * A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialise I/O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% initialise output
u2 = u0*rc;
u1 = (u0+(k*v0))*rc;
u  = u2;


%%% initialise scheme variables


%%%% input


%%%% Energy %%%NOTE%%% Following NSS, still not sure about Potential energy
if energyAn

  % energy coefficients = [Kinetic, Potential];
  coE = [.5*(hx*hy)/k^2 , .5*(kappa^2)*(1/hx*hy)];

  coL = [(hx*hy)*sigma0/(2*k^2), sigma1/(4*k^2)]; % Loss energy coefficients

  Energy = zeros(Nf,1); % total
  KE = Energy;          % kineteic
  PE = Energy;          % potential
  Loss = Energy;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_on
  view1 = [0 Nx 0 Ny -pTr pTr];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Main loop
tic
if run
  for n = 1:Nf

    % update input forces

    % main operation
    u = B*u1 + C*u2;

    % plotting
    if plot_on
      if ~mod(n,FrameDrop)
        surf(reshape(u,Ny,Nx));
        axis(view1);
        view(2)
        shading interp
        % F(n) = getframe(gca);
        drawnow
      end
    end

    % read output
    if (outtype==1)
      y(n,:) = u(lo);

    elseif (outtype==2)
      y(n,:) = (SR*(u(lo)-u1(lo)));

    end

    if energyAn % shift energy state
      KE(n) = coE(1)*sum((u-u1).^2);
      PE(n) = coE(2)*((LA*u)' * (LA*u1));

      Loss(n) = coL*[(u-u2)'*(u-u2); (GR*(u-u1))' * (GR*(u1-u2))];

      Energy(n) = coE*[(u-u1)'*(u-u1); (LA*u)' * (LA*u1)];

    end
    % shift state
    u2 = u1; u1 = u;


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

p = [1:10];
q = [1:10]';

[P,Q] = meshgrid(p,q);

mfreqs = (pi*kappa/(2*L^2))*((P.^2) + (Q.^2));     % mode frequencies


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
font_s = 14; % font point size
set(0,'DefaultAxesFontSize',font_s);

% Output FFT

if modeAn
  YF = abs(fft(y));
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
end

if energyAn
  %change in energy

  dEdt = (Energy(2:end-1)-Energy(1:end-2))/k;
  figure(2)
  subplot(2,1,1)
  plot(Energy);
  hold on
  plot(dEdt)
  plot(KE)
  plot(PE)
  plot(Loss)
  ax2 = gca;
  title('Change in Energy at Time Step k', 'FontSize', font_s)
  xlabel('Sample_{k}');
  ylabel('Power_{Not Watts}');
  set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen
  legend('Total Energy','Power_{dE/dt}','Kinetic','Potential', 'Loss')


  subplot(2,1,2)
  plot(dEdt,'.')
  hold on
  plot(Loss(2:end-1), '.')
  plot(dEdt+Loss(2:end-1),'.')
  title('dE/dt rounding error', 'FontSize', font_s)
  xlabel('Sample_{k}');
  ylabel('dE/dt');
end


% if plot_on
%   v = VideoWriter('2DWave','Motion JPEG AVI');
%   % v.CompressionRatio = 3;
%   v.FrameRate = 24;
%   open(v)
%   writeVideo(v,F(1:4000))
%   close(v)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if benchtest
  profile viewer
  profile off
end

% EOF
