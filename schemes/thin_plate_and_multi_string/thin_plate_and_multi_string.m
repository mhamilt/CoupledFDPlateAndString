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

plate_and_multi_string_vars % file containing all relevent parameters for this set of schemes
SR = OSR*SR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Error checking complete')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plate Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Motion Coefficients
pl_D = (pl_E*(pl_H)^3)/(12*(1-(nu^2)));
pl_kappa = sqrt(pl_D / (pl_rho*pl_H) );

%%%%% Scheme Spacing
k = 1/SR;                           % time step
pl_hmin = 2*sqrt(k*pl_kappa);       % NSS (equation 12.5)
pl_N = floor(pl_L/pl_hmin);         % number of segments
pl_h = pl_L/(pl_N);                 % adjusted grid spacing
pl_mu = (pl_kappa * k)/(pl_h^2);    % scheme parameter
Nf = floor(SR*Tf);                  % number of time steps

pl_N = pl_N+1;                      % for includng only internal grid points
ss = pl_N*pl_N;                     % total grid size.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% String Derived Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% derived parameters

st_A = pi*st_r^2;                       % string cross-sectional area
st_I = 0.25*pi*st_r^4;                  % string moment of intertia

st_c = sqrt(st_T/(st_rho*st_A));        % wave speed
if stiff
    st_K = sqrt(st_E*st_I/(st_rho*st_A));   % stiffness constant
else
    st_K = zeros(strNum,1);
end
%%%%% grid

st_hmin = st_c*k;
st_N = floor(st_L./st_hmin);          % number of segments (N+1 is number of grid points)
st_h = st_L./st_N;                    % adjusted grid spacing
st_lambda = st_c*k./st_h;             % Courant number

if stiff
    st_mu = st_K*k./st_h.^2;               % numerical stiffness constant
else
    st_mu = zeros(strNum,1);
end


st_N = st_N+1;                       % change st_N to be number of grid point
st_ss = sum(st_N);
%%%%% I/O
li = floor(st_xi*st_N);         % grid index of excitation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Error Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Array of variable names and conditions they must meet.

if ~skpErr
    paraLs = {'SR' 'strNum'};
    paraCons = [(SR<1e3), strNum~=length(unique(floor(pl_ctr*pl_N)))];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % String Loss coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if losstype ==0
    % frequency independant loss
    st_sigma0 = 0;
    st_sigma1 = 0;

    pl_sigma0 = 0;
    pl_sigma1 = 0;
end

if losstype ==1
    % frequency independant loss
    st_sigma0 = (6*log(10)/st_loss(1,2));
    st_sigma1 = zeros(strNum,1);

    pl_sigma0 = 6*log(10)/pl_loss(1,2);
    pl_sigma1 = 0;
end

if losstype == 2

    %%% string loss coefficients
    st_sigma0 = (6*log(10)/st_loss(1,2));
    st_sigma1 = zeros(strNum,1);


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

lo = rp*pl_N;
lo = [floor(sub2ind([pl_N pl_N],lo(1), lo(3))), floor(sub2ind([pl_N pl_N],lo(2), lo(4)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IV Read In/Out Error Check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check that xo and xi are greater than h distance away from boundaries
if any([(st_xi*st_L)<st_h (st_L-(st_xi*st_L))<st_h])
    warning([sprintf('xo is too close to string boundary\n'),...
        sprintf('changing xo to equal h\n')]);
    st_xi = abs(round(st_xi)-st_h);
    if li<=2
        li = 3;
    elseif li>=st_N
        li = st_N-1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Force Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st_d0 = (k^2./(st_h.*st_rho*st_A));          % force foefficient

% create force signal
st_f = sparse(zeros(Nf,strNum));                        % input force signal
durint = floor(st_dur*SR);                 % duration of force signal, in samples
exc_st_int = floor(exc_st*SR);
exc_end_int = exc_st_int+durint-1;
for hit = 1:strNum
    n = (exc_st_int(hit):exc_end_int(hit))+1;
    st_f(n,hit)=st_famp * 0.5 * ...
        (1-cos((2/itype)*pi.*(n'-exc_st_int(hit))/durint))*st_d0(hit);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plate Coefficient Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LA = laplace2(pl_N,pl_N,bctype); % Laplacian Matrix
BH = biharm2(pl_N,pl_N,bctype); % biharmonic matrix
pl_mI = fidimat(pl_N,pl_N,'I');

pl_mA = (1/(1+k*pl_sigma0))*speye(ss);
pl_mB = (-(pl_mu^2)*BH + 2*pl_mI + (2*k*pl_sigma1/pl_h^2)*LA) * pl_mA;
pl_mC = (-(1-pl_sigma0*k)*pl_mI - (2*k*pl_sigma1/pl_h^2)*LA) * pl_mA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% String Coefficients Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

countN = zeros(1,strNum+1);

for s = 1:strNum

    %%%% scheme coefficients
    Dn = fidimat(st_N(s),'x+');
    Dnn = fidimat(st_N(s),'xx');

    Dnn(1,:) = Dn(1,:);

    %%% String sparsity matrix
    st_mA = (1/(1+k*st_sigma0))*speye(st_N(s));
    st_mB{s} = (2*speye(st_N(s))+(st_lambda(s)^2)*Dnn) * st_mA ;
    st_mC{s} = -(1-st_sigma0*k)*speye(st_N(s)) * st_mA;

    countN(s+1) = countN(s)+st_N(s); % increment count of grid size

end

st_mB = blkdiag(st_mB{:});
st_mC = blkdiag(st_mC{:});
st_mI = speye(sum(st_N));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Couple Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% coupling points
pl_w0 = sub2ind([pl_N, pl_N],floor(pl_ctr*pl_N),floor(pl_ctrx*pl_N)*ones(1,strNum));
% pl_wL = sub2ind([pl_N, pl_N],floor(pl_ctr(2,1)*pl_N),floor(pl_ctr(2,2)*pl_N));

%%% Interpolation off set
pl_ax = (pl_ctrx*pl_N) - floor(pl_ctrx*pl_N);
pl_ay = (pl_ctr*pl_N) - floor(pl_ctr*pl_N);

%%% Mass Ratio Coefficients
st_Mr = (pl_rho*pl_H*pl_h^2)*(1+k*pl_sigma0)./((pl_rho*pl_H*(pl_h^2)*(1+k*pl_sigma0))+(st_rho*st_A*st_h*(1+k*st_sigma0)));
pl_Mr = (st_rho*st_A*st_h)*(1+k*st_sigma0)./((pl_rho*pl_H*(pl_h^2)*(1+k*pl_sigma0))+(st_rho*st_A*st_h*(1+k*st_sigma0)));

%%% Spreading vector
J = sparse(ss,strNum);
J([pl_w0; pl_w0+1; pl_w0+pl_N;pl_w0+pl_N+1]) = ...
    [(1-pl_ax)*(1-pl_ay);(1-pl_ax)*pl_ay;pl_ax*(1-pl_ay);pl_ax*pl_ay];

% Main Sparisty Matrices
B = blkdiag(pl_mB,st_mB);
C = blkdiag(pl_mC,st_mC);
MS = size(B); MS = MS(1);
%%% Coupling Vector
% for s = 1:strNum
J = sparse(ss,strNum); J(pl_w0+(0:strNum-1)*ss) = 1;

%%% Coupling Vector
%   cpl_v = [-pl_mB(pl_w0(s),:),st_mB(1+countN(s),:)];        % coupling vector for B matrix
%   cpl_vc = [(2*pl_mI(pl_w0(s),:) + (2*pl_sigma1*k/pl_h^2)*LA(pl_w0(s),:))/(1+k*pl_sigma0), -2*st_mI(1+countN(s),:)/(1+k*st_sigma0)];   % coupling vector for C matrix
%     cpl_v = [-J0'*pl_mB,st_mB(1+countN(s),:)];        % coupling vector for B matrix
%     cpl_vc = [J0'*(2*pl_mI + (2*pl_sigma1*k/pl_h^2)*LA)/(1+k*pl_sigma0), -2*st_mI(1+countN(s),:)/(1+k*st_sigma0)];   % coupling vector for C matrix

F = [-J'*pl_mB,st_mB(1+countN(1:end-1),:)];        % coupling vector for B matrix
Fc = [J'*(2*pl_mI + (2*pl_sigma1*k/pl_h^2)*LA)/(1+k*pl_sigma0), -2*st_mI(1+countN(1:end-1),:)/(1+k*st_sigma0)];   % coupling vector for C matrix


%   % %%%% Add coupling to Matrices
%   B(pl_w0(s),:) = B(pl_w0(s),:) + pl_Mr(s)*cpl_v;
%   B(ss+1+countN(s),:) = B(ss+1+countN(s),: ) - st_Mr(s)*cpl_v;
% B(pl_w0,:) = B(pl_w0,:) + repmat(pl_Mr,1,MS).*F;
% B(ss+1+countN(1:end-1),:) = B(ss+1+countN(1:end-1),: ) - repmat(st_Mr,1,MS).*F;
%
%   % %%%careful of signs, exploding systems tend to be due to backward time step
%   C(pl_w0(s),:) = C(pl_w0(s),:) + pl_Mr(s)*cpl_vc;
%   C(ss+1+countN(s),: ) = C(ss+1+countN(s),: ) - st_Mr(s)*cpl_vc;
%
% C(pl_w0,:) = C(pl_w0,:) + repmat(pl_Mr,1,MS).*Fc;
% C(ss+1+countN(1:end-1),: ) = C(ss+1+countN(1:end-1),: ) - repmat(st_Mr,1,MS).*Fc;
pl_F  = repmat(pl_Mr,1,MS).*F;
pl_Fc = repmat(pl_Mr,1,MS).*Fc;
st_F  = -repmat(st_Mr,1,MS).*F;
st_Fc = -repmat(st_Mr,1,MS).*Fc;


% end
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
fvect = [u;w];                     % force vector
li = li' + countN(1:end-1);
%%%% output
y = zeros(Nf,2);
Js = sparse(st_ss,strNum);
Js(1+countN(1:end-1)+(0:strNum-1)*st_ss) = 1;


%%%% Energy %%%NOTE%%% Following NSS, still not sure about Potential energy
if energyAn
    % energy coefficients = [Kinetic, Potential];
    pl_coE = [.5*(pl_h^2)/k^2 , .5*(pl_kappa^2)/(pl_h^2)];
    st_coE = [.5*(st_h)/k^2 , .5*(st_c^2)/st_h];

    pEnergy = zeros(Nf,1); % total
    pKE = pEnergy;          % kineteic
    pPE = pEnergy;          % potential

    sEnergy = zeros(Nf,1); % total
    sKE = sEnergy;          % kineteic
    sPE = sEnergy;          % potential

    fEnergy = sEnergy;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('All variables initialised')

%%%%% Load Bar
loadBar = waitbar(0, 'Processing...');
loadDrop = floor(Nf*.01);
oNf = 1/Nf;

%%%%% Main loop
tic
if run
    for n = 1:Nf

        %%% update input forces
        fvect(ss + li) = st_f(n,:);

        %%% main operation
        %         uw = B*uw1 + C*uw2 + fvect;

        uw = B*uw1 + C*uw2 + fvect + [J*pl_F*uw1;Js*st_F*uw1] + [J*pl_Fc*uw2;Js*st_Fc*uw2];
%         uw = B*uw1 + C*uw2 + fvect + [J;Js]*F*uw1 + [J;Js]*Fc*uw2;


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
                subplot(2,1,1)
                surf(reshape(uw(1:ss) ,pl_N,pl_N));
                zlim([-pTr pTr])
                % shading interp
                % view(2)
                %         axis square

                %%% string
                subplot(2,1,2)
                % plot([1:sum(st_N)]'*st_h, uw(ss+1:end), 'k');
                % axis([0 sum(st_L) -sTr sTr])
                plot(uw(ss+1:end), 'k');
                ylim([-sTr sTr])
                line([countN(:), countN(:)], .25*[-sTr sTr],'Color',[0 0 0],'LineWidth',2)
                % surf(sZ+3,sY+3,sX+uw(ss+1:end))
                % axis([0 st_L -sTr sTr -sTr sTr])
                % axis square
                set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen

                % F(n) = getframe(gca);
                if any(isnan(uw))
                    break
                end
                drawnow
            end
        end

        if energyAn % shift energy state
            u = uw(1:ss); u1 = uw1(1:ss);
            w = uw(ss+1:end); w1 = uw1(ss+1:end);

            pKE(n) = pl_coE(1)*((u-u1)'*(u-u1)); pPE(n) = pl_coE(2)*((LA*u)' * (LA*u1));
            pEnergy(n) = pl_coE*[sum(((u-u1).^2)); (LA*u)' * (LA*u1)];

            sKE(n) = st_coE(1)*((w-w1)'*(w-w1)); sPE(n) = st_coE(2)*((Dn*w)' * (Dn*w1));
            sEnergy(n) = st_coE*[sum(((w-w1).^2)); (Dn*w)' * (Dn*w1)];
            fEnergy(n) = (k*st_c^2/(2*st_h^2))*(w(2)-w(1)) + SR*(w(1)-w1(1)) + (pl_kappa^2*k*.5)*pl_mB(pl_w0,:)*u - SR*(u(pl_w0) - u1(pl_w0));
        end

        %%% shift state
        uw2 = uw1; uw1 = uw;

        if ~mod(n,loadDrop)
            waitbar(n*oNf, loadBar);
        end

    end
end
toc
close(loadBar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis Module for Thin Plate and string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plate_string_analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_file = sprintf('%s.wav',file);
audiowrite(save_file, (y/abs(max(y(:)))), SR);

if benchtest
    profile viewer
    profile off
end

% EOF
