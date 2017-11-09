%%%%% FDTD 2D Wave Model
%%%%% Matthew Hamilton s0674653
%%%%% Description:
%%%%%
%%%%% A Kirchhoff Thin Plate FDTD Model

clear all
close all hidden

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Instrument File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plateVerbVars % file containing all relevent parameters for this set of schemes

%%%%%%%%% DONT EDIT THESE %%%%%%%
% SR = SR*OSR;                        % redefine SR by OSR

filetypes = {'*.wav';'*.flac';'*.ogg';'*.mp3'; '*.m4a'; '*.mp4'};

[sfile, spath] =  uigetfile(filetypes,'Choose a Sound File'); %% No Message in OS X
if ~sfile; error('Load file cancelled'); end

[force, SR] = audioread([spath sfile]);

%%% Normalise audio
force = force/max(abs(force));
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
disp('Error checking complete')
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
Nx = floor((Lx)/hmin);          % number of segments
Ny = floor((Ly)/hmin);          % number of segments
h = sqrt(Lx*Ly)/sqrt(Nx*Ny);    % adjusted grid spacing
mu = (kappa * k)/(h^2);         % scheme parameter
Nx = Nx+1;                      % for includng only internal grid points %%NOTE%%
Ny = Ny+1;
ss = Nx*Ny;                     % total grid size.
%%%%% I/O
Nf = length(force)+ zpad*SR;          % number of time steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % Loss coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if losstype ==1
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
lo = rp.*[Nx Ny;Nx Ny];
lo = [floor(sub2ind([Nx Ny],lo(1), lo(3))), floor(sub2ind([Nx Ny],lo(2), lo(4)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Force Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d0 = (k^2)/(rho*H*(h^2))*(1/(1+k*sigma0)); % input force coefficient
force = [force;zeros(zpad*SR,1)]*d0;       % input oscillator
fIn = floor(sub2ind([Ny Nx],ctr(1)*Ny, ctr(2)*Nx));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficient Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BH = biharm2(Ny,Nx,bctype); % biharmonic matrix
LA = laplace2(Ny,Nx,bctype); % Laplacian matrix

A = (1/(1+k*sigma0))*speye(ss);   %NOTE% Currently inverted
B = (-(mu^2)*BH + (2*sigma1*k/(h^2))*LA + 2*speye(ss)) * A;
C = (-(2*sigma1*k/(h^2))*LA - (1-sigma0*k)*speye(ss))  * A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialise I/O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% initialise output
u2 = zeros(ss,1);
u1 = u2;
u  = u2;
Fv = sparse(u2);              % force vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_on
  fig = figure(1);
  fig.Color = [0 0 0];
  ax = axes('XLim',[0 Nx],'YLim',[0 Ny],'ZLim',[-pTr pTr]);
  colormap copper
  shading interp
  material SHINY

  set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen

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
disp('All variables initialised')

%%%%%
%%%%% Main loop
tic
if run

  %%%%% Load Bar
  loadBar = waitbar(0, 'Processing...');
  loadDrop = floor(Nf*.01);
  oNf = 1/Nf; %%% 1 over Nf for loading bar

  for n = 1:Nf

    % update input forces
    Fv(fIn) = force(n);

    % main operation
    u = B*u1 + C*u2 + Fv;

    %%% plotting
    % if plot_on
    %   if ~mod(n,FrameDrop)
    %
    %     %%% plate
    %     surf(ax,reshape(u,Ny,Nx),'FaceLighting','gouraud','FaceColor','interp',...
    %     'AmbientStrength',0.5)
    %     % light('Position',[1 0 1],'Style','local')
    %     % shading interp
    %     % view(45+spin*n/(Nf*FrameDrop),30)
    %
    %     ax.XLim =[0 max([Nx Ny])]; ax.YLim =[0 max([Nx Ny])]; ax.ZLim =[-pTr pTr];
    %     axis off
    %     % F(n/FrameDrop) = getframe(gca);
    %     drawnow
    %
    %   end
    % end

    % read output
    if (outtype==1)
      y(n,:) = u(lo);

    elseif (outtype==2)
      y(n,:) = (SR*(u(lo)-u1(lo)));

    end

    % shift state
    u2 = u1; u1 = u;
    %
    % if energyAn % shift energy state
    %   KE(n) = coE(1)*sum((u1-u2).^2); PE(n) = coE(2)*((LA*u1)' * (LA*u2));
    %   Energy(n) = coE*[sum((u1-u2).^2); (LA*u1)' * (LA*u2)];
    % end

    % if any(isnan(u))
    %   error('UNSTABLE script terminated')
    % end
    if ~mod(n,loadDrop)
      waitbar(n*oNf, loadBar)
    end


  end


    %% Save Output
    y = decimate(y,OSR);                                  %% downsample
    save_file = sprintf('%s_plate_reverb.wav',sfile);        %% save file name
    audiowrite(save_file, (y/abs(max(y(:)))), SR/OSR);    %% Normalise output

    %% Close Loading Bar
    close(loadBar)
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

freqs = (pi*kappa/(2*L^2))*((P.^2) + (Q.^2));     % mode frequencies


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
line([freqs(:) freqs(:)],ax1.YLim,'Color',[1 0 0])
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
