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
li = floor(sub2ind([N N],li(1), li(2)));
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

% BH = biharm2(N,N,bctype); % biharmonic matrix
% LA = laplace2(N,N,bctype); % Laplacian matrix
% GR = fidimat(N,N,'grad');  % gradient matrix
%
% A = (1/(1+k*sigma0))*speye(ss);   %NOTE% Currently inverted
% B = (-(mu^2)*BH + (2*sigma1*k/(h^2))*LA + 2*speye(ss)) * A;
% C = (-(2*sigma1*k/(h^2))*LA - (1-sigma0*k)*speye(ss))  * A;

% coefficients are named based on position on the x and y axes.
A00 = 1/(1+k*sigma0); % Central Loss Coeeffient (INVERTED)

%% Current time step (B) coeffients
% There are six unique coefficients for B coefs
B00 = (-(mu^2)*20 + (2*sigma1*k/(h^2))*-4 + 2) * A00;	% center
B01 = (-(mu^2)*-8 + (2*sigma1*k/(h^2))) * A00;			% 1-off
B11 = (-(mu^2)*2) * A00;									% diag
B02 = (-(mu^2)*1) * A00;									% 2-off

% Simply Supported Boundary Coefficients
BC1 = (-(mu^2)*19 + (2*sigma1*k/(h^2))*-4 + 2) * A00; % Side
BC2 = (-(mu^2)*18 + (2*sigma1*k/(h^2))*-4 + 2) * A00; % Corner

%% Previous time step (C) coeffients
C00 = (-(2*sigma1*k/(h^2))*-4 - (1-sigma0*k))  * A00;
C01 = -(2*sigma1*k/(h^2))  * A00;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialise I/O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% initialise output
% u2 = u0*rc;
% u1 = (u0+(k*v0))*rc;
u2=zeros(size(rc)); u1 = u2;

u2(li) = u0;
u1(li) = (u0+(k*v0));
u  = u2;


%%% initialise scheme variables


%%%% input
y = zeros(Nf,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cp = li;
%
% acTest = B00*u1(cp) + ...
% B01*( u1(cp-1) + u1(cp+1) + u1(cp-N) + u1(cp+N) ) + ...
% B02*( u1(cp-2) + u1(cp+2) +u1(cp-(2*N)) + u1(cp+(2*N)) ) + ...
% B11*( u1(cp-1-N) + u1(cp+1-N) +u1(cp+1+N) + u1(cp-1+N) ) + ...
% C00*u2(cp) + ...
% C01*( u2(cp-1) + u2(cp+1) + u2(cp-N) + u2(cp+N) );

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
  % Nf = 2
  for n = 1:Nf


    for xi = 2:N-3
      for yi = 3:N-2
        cp = yi+(xi * N);

        u(cp) = B00*u1(cp) + ...
        B01*( u1(cp-1) + u1(cp+1) + u1(cp-N) + u1(cp+N) ) + ...
        B02*( u1(cp-2) + u1(cp+2) +u1(cp-(2*N)) + u1(cp+(2*N)) ) + ...
        B11*( u1(cp-1-N) + u1(cp+1-N) +u1(cp+1+N) + u1(cp-1+N) ) + ...
        C00*u2(cp) + ...
        C01*( u2(cp-1) + u2(cp+1) + u2(cp-N) + u2(cp+N) );
      end
    end


    for xi=2:N-3
      %North
      cp = 2 +(xi * N); % current point
      u(cp)  = BC1*u1(cp) + ...
      B01*( u1(cp+1) + u1(cp-N) + u1(cp+N) ) + ...
      B02*( u1(cp-2) + u1(cp+2) +u1(cp-(2*N)) + u1(cp+(2*N)) ) + ...
      B11*( u1(cp+1-N) + u1(cp+1+N) ) + ...
      C00*u2(cp) + ...
      C01*( u2(cp+1) + u2(cp-N) + u2(cp+N) );

      %South
      cp = N-1 +(xi * N); % current point
      u(cp)  = BC1*u1(cp) + ...
      B01*( u1(cp-1) + u1(cp-N) + u1(cp+N) ) + ...
      B02*( u1(cp-2) + u1(cp-(2*N)) + u1(cp+(2*N)) ) + ...
      B11*( u1(cp-1-N) + u1(cp-1+N) ) + ...
      C00*u2(cp) + ...
      C01*( u2(cp-1) + u2(cp-N) + u2(cp+N) );

    end

    % y -axis sides
    for yi=3:N-2

      %West
      cp = yi+N; % current point
      u(cp)  = BC1*u1(cp) + ...
      B01*( u1(cp-1) + u1(cp+1) + u1(cp+N) ) + ...
      B02*( u1(cp-2) + u1(cp+2) + u1(cp+(2*N)) ) + ...
      B11*( u1(cp+1+N) + u1(cp-1+N) ) + ...
      C00*u2(cp) + ...
      C01*( u2(cp-1) + u2(cp+1) + u2(cp+N) );

      %East
      cp = yi + N*(N-2); % current point
      u(cp)  = BC1*u1(cp) + ...
      B01*( u1(cp-1) + u1(cp+1) + u1(cp-N) ) + ...
      B02*( u1(cp-2) + u1(cp+2) +u1(cp-(2*N)) ) + ...
      B11*( u1(cp-1-N) + u1(cp+1-N) ) + ...
      C00*u2(cp) + ...
      C01*( u2(cp-1) + u2(cp+1) + u2(cp-N) );

    end

    % %		printf("CORNER POINTS");
    % % corner points

    cp = N+2;
    u(cp) = BC2*u1(cp) + ...
    B01*( u1(cp-1) + u1(cp+1) + u1(cp-N) + u1(cp+N) ) + ...
    B02*( u1(cp+2) + u1(cp+(2*N)) ) + ...
    B11*( u1(cp-1-N) + u1(cp+1-N) +u1(cp+1+N) + u1(cp-1+N) ) + ...
    C00*u2(cp) + ...
    C01*( u2(cp-1) + u2(cp+1) + u2(cp-N) + u2(cp+N) );

    cp = 2*N-1;
    u(cp) = BC2*u1(cp) + ...
    B01*( u1(cp-1) + u1(cp+1) + u1(cp-N) + u1(cp+N) ) + ...
    B02*( u1(cp-2) + u1(cp+(2*N)) ) + ...
    B11*( u1(cp-1-N) + u1(cp+1-N) +u1(cp+1+N) + u1(cp-1+N) ) + ...
    C00*u2(cp) + ...
    C01*( u2(cp-1) + u2(cp+1) + u2(cp-N) + u2(cp+N) );

    cp = N*(N-2)+2;
    u(cp) = BC2*u1(cp) + ...
    B01*( u1(cp-1) + u1(cp+1) + u1(cp-N) + u1(cp+N) ) + ...
    B02*( u1(cp+2) + u1(cp-(2*N)) ) + ...
    B11*( u1(cp-1-N) + u1(cp+1-N) +u1(cp+1+N) + u1(cp-1+N) ) + ...
    C00*u2(cp) + ...
    C01*( u2(cp-1) + u2(cp+1) + u2(cp-N) + u2(cp+N) );

    cp = N*(N-1) - 1;
    u(cp) = BC2*u1(cp) + ...
    B01*( u1(cp-1) + u1(cp+1) + u1(cp-N) + u1(cp+N) ) + ...
    B02*( u1(cp-2) + u1(cp-(2*N)) ) + ...
    B11*( u1(cp-1-N) + u1(cp+1-N) +u1(cp+1+N) + u1(cp-1+N) ) + ...
    C00*u2(cp) + ...
    C01*( u2(cp-1) + u2(cp+1) + u2(cp-N) + u2(cp+N) );

    % read output
    if (outtype==1)
      y(n,:) = u(lo);

    elseif (outtype==2)
      y(n,:) = (SR*(u(lo)-u1(lo)));

    end


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


    u(li)

    % shift state
    u2 = u1; u1 = u;

    if any(isnan(u))
      error('UNSTABLE script terminated')
    end

  end
  toc
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% EOF
