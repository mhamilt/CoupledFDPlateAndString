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

%%% redefine variables
pl_Nx = pl_Nx +1;
pl_Ny = pl_Ny +1;
ss = pl_Nx*pl_Ny;                     % total grid size.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% String Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% derived parameters

st_A = pi*st_r.^2;                             % string cross-sectional area
st_I = 0.25*pi*st_r.^4;                        % string moment of intertia
st_T = st_T .*st_TT;  %% compensate tuning

st_c = sqrt(st_T./(st_rho.*st_A));         % wave speed
st_K = sqrt(st_E.*st_I./(st_rho.*st_A));   % stiffness constant

%%%%% grid

st_hmin = sqrt(0.5* (st_c.^2*k^2+sqrt(st_c.^4*k^4+16*st_K.^2.*k.^2)) );
st_N = floor(st_L./st_hmin);          % number of segments (N+1 is number of grid points)
st_h = st_L./st_N;                    % adjusted grid spacing
st_lambda = st_c*k./st_h;             % Courant number
st_mu = st_K*k./st_h.^2;               % numerical stiffness constant
st_N = st_N+1;                       % change st_N to be number of grid points
st_ss = sum(st_N);                     % total number of string grid points
%%%%% I/O
li = floor(st_xi*st_N(notes));         % grid index of excitation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Error Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Array of variable names and conditions they must meet.

if (~skpErr && ~interpJ)

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

  %%% String Loss
  st_z1 = (-st_c.^2 + sqrt(st_c.^4 + 4*st_K.^2.*(2*pi*st_loss(1,:)).^2))./(2*st_K.^2);
  st_z2 = (-st_c.^2 + sqrt(st_c.^4 + 4*st_K.^2.*(2*pi*st_loss(3,:)).^2))./(2*st_K.^2);

  st_sigma0 = 6*log(10)*(-st_z2./st_loss(2,:) + st_z1./st_loss(4,:))./(st_z1-st_z2);
  st_sigma1 = 6*log(10)*(1./st_loss(2,:) - 1./st_loss(4,:))./(st_z1-st_z2);

  %%% Plate Loss
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

% Check that xo and xi are greater than h distance away from boundaries
if any([(st_xi*st_L)<st_h (st_L-(st_xi*st_L))<st_h])
  warning([sprintf('xo is too close to string boundary\n'),...
  sprintf('changing xo to equal h\n')]);
  % st_xi = abs(round(st_xi)-st_h);
  if any(li<=2)
    li(li<=2) = 3;
  elseif any(li>=st_N(notes))
    li(li>=st_N(notes)) = st_N(li>=st_N(notes))-1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Force Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st_d0 = (k^2./(st_h.*st_rho.*st_A.*(1+ k.*st_sigma0))); % force foefficient

st_f = sparse(Nf,length(notes));      % input force signal

hitN = length(exc_st(:,1));            % number of note excitations

durint = floor(st_dur*SR);            % duration of force signal, in samples

exc_st_int = floor(exc_st(:,1)*SR) + randi(jiggle,hitN,1);          % start time index for excitation

if itype == 1;

  %%% scaling force [min max]
  tr = [1.15e-3, 2.5e-3]; % milliseconds, time range
  fr = [1 40];          % force range, Newtons.
  vs = [1 127];          % range of velocities

end

if itype == 2;

  %%% scaling force [min max]
  tr = [4.15e-3, 5.5e-3]; % milliseconds, time range
  fr = [1 40];          % force range, Newtons.
  vs = [1 127];          % range of velocities

end


for hit = 1:hitN

  % string index of forced string
  ist = exc_st(hit,2);

  %%current note number
  cn = notes(ist);

  %% force velocity
  fv = exc_st(hit,3);

  %%% Force in Netons
  fn = (fr(2)-fr(1))*(fv-1)/(vs(2)-vs(1)) + fr(1);

  %% sample duration of force
  int_dur = floor((((tr(2)-tr(1))*(fn - fr(1))/(fr(2) - fr(1)))+tr(1))*SR);

  % excitation interval in samples depending (string depndant)
  exc_int = 0:int_dur;

  % excitation on set time
  exc_t = exc_st_int(hit)+exc_int;

  %% raised cosine
  rc = (1-cos((2/itype)*pi*(exc_int)/int_dur))';

  %%% loop through and change sample values in force matrix
  %%% + st_f... accomodates for rogue MIDI notes that
  %%% are very close together.
  st_f(exc_t,ist) = fn * 0.5 * rc *st_d0(cn) + st_f(exc_t,ist);

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

countN = zeros(1,strNum+1);

for s = 1:strNum
  %%%%% scheme coefficients
  Dn = fidimat(st_N(s),'x+',bctype);
  Dnn = fidimat(st_N(s),'xx',1);
  Dnnnn = fidimat(st_N(s),'xxxx',st_bctype);
  Dnnnn(2,1) = -2;
  st_mI = speye(st_N(s));

  %%% String sparsity matrix, uncoupled.
  st_mA = (1/(1+k*st_sigma0(s)));
  st_mB{s} = ( (st_lambda(s)^2)*Dnn - (st_mu(s)^2)*Dnnnn + (2*st_mI) +...
  (2*st_sigma1(s)*k/st_h(s)^2)*Dn ) * st_mA ;
  st_mC{s} = -((1-st_sigma0(s)*k)*st_mI + (2*st_sigma1(s)*k/st_h(s)^2)*Dnn) * st_mA;

  %%%%% Alter matrices for coupling string boundary
  st_mB{s}([1 end],:) = 0;st_mC{s}([1 end],:) = 0;
  st_mB{s}(:,end) = 0; st_mC{s}(:,end) = 0;
  st_mB{s}(1,:) = ((st_lambda(s)^2)*Dn(1,:) - (st_mu(s)^2)*Dnn(2,:) + (2*st_mI(1,:))) * st_mA;
  st_mC{s}(1,:) = -((1-st_sigma0(s)*k)*st_mI(1,:));

  countN(s+1) = countN(s)+st_N(s); % increment count of grid size

end

st_mB = blkdiag(st_mB{:});
st_mC = blkdiag(st_mC{:});
st_mI = speye(sum(st_N));

%%% Matrix index for debugging
mi = ss + countN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Couple Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% coupling points
pl_w0 = sub2ind([pl_Ny, pl_Nx],floor(pl_ctr*pl_Ny),floor(pl_rx*pl_Nx)*ones(1,strNum));
% pl_wL = sub2ind([pl_N, pl_N],floor(pl_ctr(2,1)*pl_N),floor(pl_ctr(2,2)*pl_N));
st_w0 = countN(1:end-1) + 1;

%%% Interpolation off set
pl_ax = (pl_rx*pl_Nx) - floor(pl_rx*pl_Nx);
pl_ay = (pl_ctr*pl_Ny) - floor(pl_ctr*pl_Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Spreading vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Total Number of Grid points
Nt = ss+ st_ss;

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

Mr = inv( spdiags(st_Mr', 0, strNum,strNum) +...
(1/(pl_rho*pl_H*(pl_h^2)*(1+k*pl_sigma0)))*(pl_J'*pl_J));

F = Mr*[-pl_J'*pl_mB,st_J'*st_mB];        % coupling vector for B matrix
% F1 = Mr*[pl_J'*(2*pl_mI + (2*pl_sigma1*k/pl_h^2)*LA)/(1+k*pl_sigma0),...
% -2*inv(spdiags([1+k*st_sigma0]',0,strNum,strNum))*(st_J'*st_mI)];   % coupling vector for C matrix
F1 = Mr*[pl_J'*(((1+pl_sigma0*k)/(1+k*pl_sigma0))*pl_mI - pl_mC),...
-2*inv(spdiags([1+k*st_sigma0]',0,strNum,strNum))*(st_J'*st_mI)];   % coupling vector for C matrix


%%% an alternative to F1
% F1 = pl_J'*-(pl_mC+((1-pl_sigma0*k)/(1+k*pl_sigma0))*pl_mI - (2/(1+k*pl_sigma0))*pl_mI)
% F1 = pl_J'*(((1+pl_sigma0*k)/(1+k*pl_sigma0))*pl_mI - pl_mC);

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

%%% String boundary coupling
Js = sparse(st_ss,strNum);
% Js(1+countN(1:end-1)+(0:strNum-1)*st_ss) = 1;
Js(1+countN(1:end-1)+(0:strNum-1)*st_ss) = 1./(st_rho.*st_A.*st_h);

%%%%% Spreading operator
Jp = J/(pl_rho*pl_H*pl_h^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Visualisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_on
  fig = figure(1);
  fig.Color = [0 0 0];
  ax = axes('XLim',[0 pl_Nx],'YLim',[0 pl_Ny],'ZLim',[-pTr pTr+sTr]);
  colormap copper
  material METAL
  lighting gouraud
  set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen
  view(15,15)

  for s = 1:strNum

    % st_x0{s} = linspace(pl_ctr(s)*pl_Nx,pl_ctr(s)*pl_Nx+st_L(s),st_N(s));
    st_x0{s} = linspace(pl_rx*pl_Nx,pl_rx*pl_Nx+st_L(s)/pl_h,st_N(s));

    % st_y0{s} = pl_rx*pl_Nx*ones(1,st_N(s))+s;
    st_y0{s} = pl_ctr(s)*pl_Ny*ones(1,st_N(s));

    % st_y0{s} = linspace(floor(pl_ctr(1,2)*pl_N),pl_N,st_N);
    % st_x0{s} = pl_ctr(1,1)*pl_N*ones(1,st_N);
  end
  pl_ylim = max(pl_rx*pl_Nx+st_L/pl_h);
end

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

  tic
  for n = 1:Nf


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

        %%% plate
        % colormap copper
        surf(reshape(uw(1:ss) ,pl_Ny,pl_Nx),'FaceLighting','gouraud','FaceColor','interp',...
        'AmbientStrength',0.5)
        light('Position',[1 0 1],'Style','local')
        colormap copper
        ax.Projection = 'perspective';
        shading interp
        material SHINY
        % caxis([-3e-6,3e-6])

        %%% string
        hold on
        for s=1:strNum
          plot3(st_x0{s},st_y0{s},uw(ss+1+countN(s):ss+countN(s+1))*sTr/3e-2 +...
          sTr+2e-6,'r','LineWidth',.5)
        end
        hold off

        % view(15,15)
        view(45+spin*n/(Nf*FrameDrop),15)
        axis off
        ax.XLim =[0 pl_Nx]; ax.YLim =[0 pl_Ny]; ax.ZLim =[-pTr pTr+sTr];
        drawnow
        % F(n) = getframe(gca);
      end
    end

    %%% shift state
    uw2 = uw1; uw1 = uw;

    if ~mod(n,loadDrop)
      waitbar(n*oNf, loadBar)
    end

    if any(isnan(uw))
      error('UNSTABLE script terminated')
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
