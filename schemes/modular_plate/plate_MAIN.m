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

plate_scheme_vars % file containing all relevent parameters for this set of schemes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Error Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Array of variable names and conditions they must meet.

if ~skpErr
  plate_error
  disp('Error checking complete')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plate_coefs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % Read In/Out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lo = rp*N;
lo = floor(sub2ind([N N],lo(1), lo(2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Force Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y] = meshgrid([1:N]*h, [1:N]*h);         % Grid of point in value of meters
dist = sqrt((X-(ctr(1)*L)).^2 + (Y-(ctr(2)*L)).^2); % distance of points from excitation
ind = sign(max(-dist+(wid*0.5),0));         % displaced grid points (logical)
rc = .5*ind.*(1+cos(2*pi*dist/wid));        % displacement
rc = rc(:);                                 % 2D plane as vector


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficient Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LA = laplace2(N,N,bctype); % Laplacian Matrix
BH = biharm2(N,N,bctype); % biharmonic matrix

B = -(mu^2)*BH + 2*speye(ss);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialise I/O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% initialise output
u2 = u0*rc;
u1 = (u0+(k*v0))*rc;
u  = u2;
y = zeros(Nf,1);

%%%% input



%%%% Energy %%%NOTE%%% Following NSS, still not sure about Potential energy
if energyAn

  % energy coefficients = [Kinetic, Potential];
  coE = [.5*(h^2)/k^2 , .5*(kappa^2)*(1/h^2)];

  Energy = zeros(Nf,1); % total
  KE = Energy;          % kineteic
  PE = Energy;          % potential

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_on
  view = [0 N 0 N -0.00015 0.00015];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('All variables initialised')

%%%%% Main loop
tic
if run
  for n = 1:Nf

    % update input forces


    % main operation
    u = B*u1 - u2;


    % plotting
    if plot_on
      mesh(reshape(u,N,N));
      axis(view);
      F(n) = getframe(gca);
      drawnow
    end

    % read output
    if (outtype==1)
      y(n,:) = u(lo);

    elseif (outtype==2)
      y(n,:) = (SR*(u(lo)-u1(lo)));

    end

    % shift state
    u2 = u1; u1 = u;

    if energyAn % shift energy state
      KE(n) = coE(1)*sum((u1-u2).^2); PE(n) = coE(2)*((LA*u1)' * (LA*u2));
      Energy(n) = coE*[sum((u1-u2).^2); (LA*u1)' * (LA*u2)];
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

if modeAn
  p = [1:5];
  q = [1:5]';

  [P,Q] = meshgrid(p,q);

  mfreqs = (pi*kappa/(2*L^2))*((P.^2) + (Q.^2));     % mode frequencies
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output
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
  dEdt = (Energy(2:end)-Energy(1:end-1));
  figure(2)
  subplot(2,1,1)
  plot(Energy);
  hold on
  plot(dEdt)
  plot(KE)
  plot(PE)
  ax2 = gca;
  title('Change in Energy at Time Step k', 'FontSize', font_s)
  xlabel('Sample_{k}');
  ylabel('Power_{Not Watts}');
  set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen
  legend('Total Energy','Power_{dE/dt}','Kinetic','Potential')


  subplot(2,1,2)
  plot(dEdt,'.')
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
