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

D = (E*(H)^3)/(12*(1-(nu^2)));

kappa = sqrt(D / (rho*  H) );

%%%%% Scheme Spacing
k = 1/SR;                   % time step
hmin = 2*sqrt(k*kappa);      % NSS (equation 12.5)
N = floor(L./hmin);          % number of segments
h = L./(N);                    % adjusted grid spacing (only internal points)
mu = (kappa * k)/(h^2);      % scheme parameter
Nf = floor(SR*Tf);          % number of time steps

N = N+1;
ss = N*N;                    % total grid size.
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
lo = rp*N;
lo = floor(sub2ind([N N],lo(1), lo(2)));
li = ctr*N;
li = floor(sub2ind([N N],li(1), li(2)))
%1+(N*(lo(2)-1)) + lo(1);
% 1+(N*( ( (rp(2)*N)-1)) ) + (rp(1)*N)-1;
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

BH = biharm2(N,N,bctype); % biharmonic matrix
LA = laplace2(N,N,bctype); % Laplacian matrix
GR = fidimat(N,N,'grad');  % gradient matrix

A = (1/(1+k*sigma0))*speye(ss);   %NOTE% Currently inverted
B = (-(mu^2)*BH + (2*sigma1*k/(h^2))*LA + 2*speye(ss)) * A;
C = (-(2*sigma1*k/(h^2))*LA - (1-sigma0*k)*speye(ss))  * A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialise I/O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% initialise output
u2 = u0*rc;
u1 = (u0+(k*v0))*rc;
u  = u2;


%%% initialise scheme variables


%%%% input
y = zeros(Nf,1);

%%%% Energy %%%NOTE%%% Following NSS, still not sure about Potential energy
if energyAn

  % energy coefficients = [Kinetic, Potential];
  coE = [.5*(h^2)/k^2 , .5*(kappa^2)*(1/h^2)];

  coL = [(h^2)*sigma0/(2*k^2), sigma1/(4*k^2)]; % Loss energy coefficients

  Energy = zeros(Nf,1); % total
  KE = Energy;          % kineteic
  PE = Energy;          % potential
  Loss = Energy;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_on
  fig = figure(1);
  fig.Color = [0 0 0];
  ax = axes('XLim',[0 N],'YLim',[0 N],'ZLim',[-pTr pTr]);
  colormap copper
  shading interp
  material SHINY

  set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen

end
%%%%% Main loop
tic
if run
  for n = 1:Nf

    % main operation
    u = B*u1 + C*u2;

    % plotting
    if plot_on
      if ~mod(n,FrameDrop)

        %%% plate
        surf(ax,reshape(u,N,N),'FaceLighting','gouraud','FaceColor','interp',...
        'AmbientStrength',0.5)
        light('Position',[1 0 1],'Style','local')
        shading interp
        % view(45+spin*n/(Nf*FrameDrop),30)

        ax.XLim =[0 N]; ax.YLim =[0 N]; ax.ZLim =[-pTr pTr];
        axis off
        F(n/FrameDrop) = getframe(gca);
        drawnow

      end
    end

    % read output
    if (outtype==1)
      y(n,:) = u(lo);

    elseif (outtype==2)
      y(n,:) = (SR*(u(lo)-u1(lo)));

    end

    % if energyAn % shift energy state
    %   KE(n) = coE(1)*sum((u-u1).^2);
    %   PE(n) = coE(2)*((LA*u)' * (LA*u1));
    %
    %   Loss(n) = coL*[(u-u2)'*(u-u2); (GR*(u-u1))' * (GR*(u1-u2))];
    %
    %   Energy(n) = coE*[(u-u1)'*(u-u1); (LA*u)' * (LA*u1)];
    %
    % end
    % shift state
    u2 = u1; u1 = u;

    % if any(isnan(u))
    %   error('UNSTABLE script terminated')
    % end

  end
end
toc

%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Modal Analysis
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% calculating plate modes
%
% p = [1:10];
% q = [1:10]';
%
% [P,Q] = meshgrid(p,q);
%
% mfreqs = (pi*kappa/(2*L^2))*((P.^2) + (Q.^2));     % mode frequencies
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% % Plotting
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% font_s = 14; % font point size
% set(0,'DefaultAxesFontSize',font_s);
%
% % Output FFT
%
% if modeAn
%   YF = abs(fft(y));
%   fy = linspace(0, SR, length(y)); % frequency axis
%   figure(1);
%   plot(fy, 20*log10(YF),'k');
%   ax1 = gca;
%   line([mfreqs(:) mfreqs(:)],ax1.YLim,'Color',[1 0 0])
%   ax1.XLim = [0, 500];
%   title(sprintf('Modes: Sample Rate %.0f_{Hz}',SR), 'FontSize', font_s)
%   xlabel('Freq._{Hz}');
%   ylabel('magnitude');
%   legend('FFT Output_{dB}','Mode Frequencies_{Hz}')
%   set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen
% end
%
% if energyAn
%   %change in energy
%
%   dEdt = (Energy(2:end-1)-Energy(1:end-2))/k;
%   figure(2)
%   subplot(2,1,1)
%   plot(Energy);
%   hold on
%   plot(dEdt)
%   plot(KE)
%   plot(PE)
%   plot(Loss)
%   ax2 = gca;
%   title('Change in Energy at Time Step k', 'FontSize', font_s)
%   xlabel('Sample_{k}');
%   ylabel('Power_{Not Watts}');
%   set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen
%   legend('Total Energy','Power_{dE/dt}','Kinetic','Potential', 'Loss')
%
%
%   subplot(2,1,2)
%   plot(dEdt,'.','Color','k')  %% Lossless
%   hold on
%   % plot(Loss(2:end-1), '.')
%   % plot(dEdt+Loss(2:end-1),'.','Color','k')
%   title('dE/dt Rounding Error', 'FontSize', font_s)
%   xlabel('Sample_{k}');
%   ylabel('dE/dt');
% end
%
%
% if plot_on
%   v = VideoWriter('2DWave','Motion JPEG AVI');
%   % v.CompressionRatio = 3;
%   v.FrameRate = 24;
%   open(v)
%   file_open = true;
%   frameBlock = 0
%   while file_open
%
%     try
%       writeVideo(v,F((1:100)+100*frameBlock));
%       frameBlock = frameBlock+1;
%     catch
%       writeVideo(v,F((1+100*(frameBlock-1)):end) ) ;
%       file_open = false;
%     end
%
%   end
%   close(v)
% end
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% if benchtest
%   profile viewer
%   profile off
% end
%
% EOF
