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

instPiano % file containing all relevent parameters for this set of schemes
SR = OSR*SR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plate Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Motion Coefficients
pl_D = (pl_E*(pl_H)^3)/(12*(1-(nu^2)));
pl_kappa = sqrt(pl_D / (pl_rho*pl_H) );

%%%%% Scheme Spacing
k = 1/SR;                           % time step
pl_hmin = 2*sqrt(k*pl_kappa);       % NSS (equation 12.5)
pl_Nx = floor((pl_Lx)/pl_hmin);          % number of segments
pl_Ny = floor((pl_Ly)/pl_hmin);          % number of segments
pl_h = sqrt(pl_Lx*pl_Ly)/sqrt(pl_Nx*pl_Ny);    % adjusted grid spacing
pl_mu = (pl_kappa * k)/(pl_h^2);    % scheme parameter
Nf = floor(SR*Tf);                  % number of time steps

pl_Nx = pl_Nx +1;
pl_Ny = pl_Ny +1;
ss = pl_Nx*pl_Ny;                     % total grid size.


%%%%% I/O
Nf = floor(SR*Tf);          % number of time steps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% String Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% derived parameters

for t = 1:3
  st_A{t} = pi*st_r{t}.^2;                             % string cross-sectional area
  st_I{t} = 0.25*pi*st_r{t}.^4;                        % string moment of intertia
  st_c{t} = sqrt(st_T{t}./(st_rho*st_A{t}));         % wave speed
  st_K{t} = sqrt(st_E{t}.*st_I{t}./(st_rho.*st_A{t}));   % stiffness constant

  %%%%% grid
  st_hmin{t} = sqrt(0.5* (st_c{t}.^2*k^2+sqrt(st_c{t}.^4*k^4+16*st_K{t}.^2.*k.^2)) );
  st_N{t} = floor(st_L{2}./st_hmin{t});          % number of segments (N+1 is number of grid points)
  st_h{t} = st_L{t}./st_N{t};                    % adjusted grid spacing
  st_lambda{t} = st_c{t}*k./st_h{t};             % Courant number
  st_mu{t} = st_K{t}*k./st_h{t}.^2;               % numerical stiffness constant
  st_N{t} = st_N{t}+1;                       % change st_N to be number of grid points
  %%%%% I/O
  li{t} = floor(st_xi*st_N{t}(notes));         % grid index of excitation

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Error Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Array of variable names and conditions they must meet.

if ~skpErr
  paraLs = {'SR' 'strNum', 'MIDI'};
  paraCons = [(SR<1e3), strNum~=length(unique(floor(pl_ctr*pl_Ny))),...
  (any(MIDI(:,1)<1) || any(MIDI(:,1)>88))];

  if any(paraCons)

    for paraError = paraLs(paraCons)

      switch char(paraError)

        %%% WARNING ERRORS %%%
      case 'strNum'
        error([sprintf('two string share the same grid point, make pl_N bigger')])

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
%%%% % String Loss coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if losstype ==0
  % frequency independant loss
  st_sigma0 = zeros(1,88);
  st_sigma1 = zeros(1,88);

  pl_sigma0 = 0;
  pl_sigma1 = 0;
end

if losstype ==1
  % frequency independant loss
  st_sigma0 = (6*log(10)/st_loss(1,2))*ones(1,88);
  st_sigma1 = zeros(1,88);

  pl_sigma0 = 6*log(10)/pl_loss(1,2);
  pl_sigma1 = 0;
end

if losstype == 2

  for t = 1
    st_z1 = (-st_c{t}.^2 + sqrt(st_c{t}.^4 + 4*st_K{t}.^2.*(2*pi*st_loss(1,1:Ts{t})).^2))./(2*st_K{t}.^2);
    st_z2 = (-st_c{t}.^2 + sqrt(st_c{t}.^4 + 4*st_K{t}.^2.*(2*pi*st_loss(3,1:Ts{t})).^2))./(2*st_K{t}.^2);

    st_sigma0{t} = 6*log(10)*(-st_z2./st_loss(2,1:Ts{t}) + st_z1./st_loss(4,1:Ts{t}))./(st_z1-st_z2);
    st_sigma1{t} = 6*log(10)*(1./st_loss(2,1:Ts{t}) - 1./st_loss(4,1:Ts{t}))./(st_z1-st_z2);
  end

  %% this is simply loss from NSS with 'c' removed as tension is not present.
  %%% plate loss coefficients
  pl_z1 = 2*pl_kappa*(2*pi*pl_loss(1,1))/(2*pl_kappa.^2); %from 1D
  pl_z2 = 2*pl_kappa*(2*pi*pl_loss(2,1))/(2*pl_kappa.^2);

  pl_sigma0 = 6*log(10)*(-pl_z2/pl_loss(1,2) + pl_z1/pl_loss(2,2))./(pl_z1-pl_z2);
  pl_sigma1 = 6*log(10)*(1/pl_loss(1,2) - 1/pl_loss(2,2))./(pl_z1-pl_z2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % Read In/Out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lo = rp.*[pl_Nx pl_Ny];
lo = [floor(sub2ind([pl_Nx pl_Ny],lo(1), lo(3))), floor(sub2ind([pl_Nx pl_Ny],lo(2), lo(4)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IV Read In/Out Error Check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 1:3
  % Check that xo and xi are greater than h distance away from boundaries
  if any([(st_xi*st_L{t})<st_h{t} (st_L{t}-(st_xi*st_L{t}))<st_h{t}])
    warning([sprintf('xo is too close to string boundary\n'),...
    sprintf('changing xo to equal h\n')]);
    % st_xi = abs(round(st_xi)-st_h);
    if any(li{t}<=2)
      li{t}(li<=2) = 3;
    elseif any(li>=st_N{t}(notes))
      li{t}(li>=st_N{t}(notes)) = st_N{t}(li{t}>=st_N{t}(notes))-1;
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Force Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 1
  st_d0 = (k^2./(st_h.*st_rho.*st_A.*(1+ k.*st_sigma0{t})));          % force foefficient
  % st_d0 = (k^2./(st_h.*st_rho.*st_A));          % force foefficient

  st_f = sparse(Nf,length(notes));      % input force signal
  hitN = length(exc_st(:,1));            % number of note excitations
  durint = floor(st_dur*SR);            % duration of force signal, in samples

  exc_st_int = floor(exc_st(:,1)*SR) + randi(jiggle,hitN,1);          % start time index for excitation


  for hit = 1:hitN

    % string index of forced string
    iStr = exc_st(hit,2);

    %% sample duration of force
    int_dur = durint(notes(iStr));

    % excitation duration in samples depending (string depndant)
    dur = (exc_st_int(hit):exc_st_int(hit)+int_dur-1)+1;

    %%% loop through and change sample values in force matrix
    st_f(dur,iStr)=st_famp * 0.5 * ...
    (1-cos((2/itype)*pi.*(dur'-exc_st_int(hit))/int_dur))*st_d0(notes(iStr)) + ...
    st_f(dur,iStr); % this accomodates for rogue MIDI notes that are very close together.

  end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plate Coefficient Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LA = laplace2(pl_Ny,pl_Nx,bctype); % Laplacian Matrix
BH = biharm2(pl_Ny,pl_Nx,bctype); % biharmonic matrix
pl_mI = fidimat(pl_Ny,pl_Nx,'I');

pl_mA = (1/(1+k*pl_sigma0))*speye(ss);
pl_mB = (-(pl_mu^2)*BH + 2*pl_mI + (2*k*pl_sigma1/pl_h^2)*LA) * pl_mA;
pl_mC = (-(1-pl_sigma0*k)*pl_mI - (2*k*pl_sigma1/pl_h^2)*LA) * pl_mA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% String Coefficients Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% countN = zeros(1,strNum+1);


%%% with this model of generating strings the lowest index will now refer
%%% to the highest string which should make some generation more simplistic

for t = 1:3 %%another loop for each tier of strings
  for s = 1:Ts{t}

    %%% current matrix cell
    cmc = s+To{t};

    %%%%% scheme coefficients
    Dn = fidimat(st_N{t}(s),'x+',bctype);
    Dnn = fidimat(st_N{t}(s),'xx',1);
    Dnnnn = fidimat(st_N{t}(s),'xxxx',st_bctype);
    Dnnnn(2,1) = -2;
    st_mI = speye(st_N{t}(s));


    %%% String sparsity matrix, uncoupled.
    st_mA = (1/(1+k*st_sigma0(s)));
    st_mB{cmc} = ( (st_lambda{t}(s)^2)*Dnn - (st_mu{t}(s)^2)*Dnnnn + (2*st_mI) +...
    (2*st_sigma1(s)*k/st_h{t}(s)^2)*Dn ) * st_mA ;
    st_mC{cmc} = -((1-st_sigma0(s)*k)*st_mI + (2*st_sigma1(s)*k/st_h{t}(s)^2)*Dnn) * st_mA;

    %%%%% Alter matrices for coupling string boundary
    st_mB{cmc}}([1 end],:) = 0;st_mC{cmc}([1 end],:) = 0;
    st_mB{cmc}(:,end) = 0; st_mC{cmc}(:,end) = 0;
    st_mB{cmc}(1,:) = ((st_lambda(s)^2)*Dn(1,:) - (st_mu{t}(s)^2)*Dnn(2,:) + (2*st_mI(1,:))) * st_mA;
    st_mC{cmc}(1,:) = -((1-st_sigma0(s)*k)*st_mI(1,:));

    % countN(s+1) = countN(s)+st_N{t}(s); % increment count of grid size

  end
  countN{t} = accrue(st_N{t},1);
end

st_mB = blkdiag(st_mB{:});
st_mC = blkdiag(st_mC{:});
st_mI = speye(sum([st_N{1},st_N{2},st_N{3}]));

%% total number of string grid points
st_ss = sum(st_N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Couple Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% coupling points
for t = 1:3
  pl_w0{t} = sub2ind([pl_Ny, pl_Nx],floor(pl_ctr(1:Ts{t})*pl_Ny),floor(pl_rx*pl_Nx)*ones(1,Ts{t})+ (t-1)*pl_Nx );
  % pl_wL = sub2ind([pl_N, pl_N],floor(pl_ctr(2,1)*pl_N),floor(pl_ctr(2,2)*pl_N));
end
%%% NOTE Very long winded an in the process of being tidied up


%%% Mass Ratio Coefficients
st_Mr = 1./(st_rho.*st_A.*st_h.*(1+k*st_sigma0));
pl_Mr = 1/((pl_rho*pl_H*pl_h^2)*(1+k*pl_sigma0));

if interpJ
  %%% Linear interpolation indeces
  lii = [pl_w0; pl_w0+1; pl_w0+pl_Nx;pl_w0+pl_Nx+1];
  %%% Linear interpolation coeffcients
  lic = [(1-pl_ax)*(1-pl_ay);(1-pl_ax)*pl_ay;pl_ax*(1-pl_ay);pl_ax*pl_ay];

  %%% Spreading Operators
  J = sparse(Nt,strNum); J(lii+(0:strNum-1)*Nt) = pl_Mr*lic/(1+k*pl_sigma0);
  J(st_w0+ss+(0:strNum-1)*Nt) = -st_Mr./(1+k*st_sigma0);
  pl_J = sparse(ss,strNum); pl_J(lii+(0:strNum-1)*ss) = lic;
  st_J = sparse(st_ss,strNum); st_J(st_w0+(0:strNum-1)*st_ss) = 1;
else
  %%% Spreading Operators
  J = sparse(Nt,strNum); J(pl_w0+(0:strNum-1)*Nt) = pl_Mr/(1+k*pl_sigma0);
  J(st_w0+ss+(0:strNum-1)*Nt) = -st_Mr./(1+k*st_sigma0);
  pl_J = sparse(ss,strNum); pl_J(pl_w0+(0:strNum-1)*ss) = 1;
  st_J = sparse(st_ss,strNum); st_J(st_w0+(0:strNum-1)*st_ss) = 1;
end

%%% Mass Ratio Matrix
Mr = inv( spdiags(st_Mr', 0, strNum,strNum) +...
(1/(pl_rho*pl_H*(pl_h^2)*(1+k*pl_sigma0)))*(pl_J'*pl_J));

%%% Coupling Vector
F = Mr*[-pl_J'*pl_mB,st_J'*st_mB];        % coupling vector for B matrix
F1 = Mr*[pl_J'*(((1+pl_sigma0*k)/(1+k*pl_sigma0))*pl_mI - pl_mC),...
-2*inv(spdiags([1+k*st_sigma0]',0,strNum,strNum))*(st_J'*st_mI)];   % coupling vector for C matrix

%%% Main Sparisty Matrices
B = blkdiag(pl_mB,st_mB) + J*F;
C = blkdiag(pl_mC,st_mC) + J*F1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialise I/O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% initialise output
u = zeros(ss,1);
w = zeros(sum(st_N),1);

%%% Joined vectors
uw = [u;w];
uw1 = uw;
uw2 = uw;

%%%% input
fvect = uw;                     % force vector

%%% for the input, off set by the number of points on the plate, ss
%%% then add the relevant number of points to offset each note
li = ss + li + countN(notes);

%%%% output
y = zeros(Nf,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('All variables initialised')


%%%%% Main loop

if run
  %%%%% Load Bar
  loadBar = waitbar(0, 'Processing...');
  loadDrop = floor(Nf*.01);
  oNf = 1/Nf; %%% 1 over Nf for loading bar

  %%% for plotting
  st_w0x = floor(pl_rx*pl_Nx)*ones(1,88);
  st_w0y = floor(pl_ctr*pl_Ny)';


  tic
  for n = 1:Nf
    % for n = 1:SR

    %%% update input forces
    fvect(li) = st_f(n,:);

    %%% main operation
    uw = B*uw1 + C*uw2 + fvect;


    % read output
    if (outtype==1)
      y(n,:) = uw(lo);

    elseif (outtype==2)
      y(n,:) = (SR*(uw(lo)-uw1(lo)));

    end


    if plot_on
      if ~mod(n,FrameDrop)
        figure(1)
        %%% plate
        % subplot(2,1,1)
        surf(reshape(uw(1:ss) ,pl_Ny,pl_Nx));
        axis([0 max([pl_Nx pl_Ny]) 0 max([pl_Nx pl_Ny]) -pTr pTr]);
        % shading interp
        % view(2)
        % axis square
        % view(45+10000*n/Nf,45)
        axis off
        hold on

        %%%%plot connecting points
        plot3(st_w0x,st_w0y, uw(pl_w0),...
        '^','MarkerSize',5,'LineWidth',1,...
        'MarkerEdgeColor','k','MarkerFaceColor','r')
        set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen
        hold off


        drawnow
      end
    end

    %%% shift state
    uw2 = uw1; uw1 = uw;

    if ~mod(n,loadDrop)
      waitbar(n*oNf, loadBar)
    end


  end

  %% Save Output
  y = decimate(y,OSR);                                  %% downsample
  save_file = sprintf('%s_pianoModel.wav',file);        %% save file name
  audiowrite(save_file, (y/abs(max(y(:)))), SR/OSR);    %% Normalise output
  close(loadBar)
end
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EOF
