% coefficient Module for plate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Motion Coefficients
kappa = sqrt((E*(H)^2)/(12*rho*(L^4)*(1-nu)));

%%%%% Scheme Spacing
k = 1/SR;                   % time step
hmin = 2*sqrt(k*kappa);      % NSS (equation 12.5)
N = floor(L./hmin)-1;          % number of segments
h = L./(N+1);                    % adjusted grid spacing (only internal points)
mu = (kappa * k)/(h^2);      % scheme parameter
Nf = floor(SR*Tf);          % number of time steps

ss = N*N;                    % total grid size.
%%%%% I/O
Nf = floor(SR*Tf);          % number of time steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % Loss coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
