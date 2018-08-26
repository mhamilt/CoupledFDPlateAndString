%% Calculate energy at a time step for a 2D FTDT scheme

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialise variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E1 = 0;   % energy at current time step
E2 = 0;   % energy at previous time step
Energy = zeros(NF,1);


% LA = laplacian matrix
E1 =  (sum(((u1 - u2)^2)/(2*(k^2)) + ((kappa^2)/(2*(h^4)*k)*LA)*(u1^2)))/k
E2 =  (sum(((u2)^2)/(2*(k^2)) + ((kappa^2)/(2*(h^4)*k)*LA)*(u2^2)))/k
dEdt = E1 - E2

%% OR: a variation on the above to reduce computation
coE = [1/(2*(k^2)), (kappa^2)/(2*(h^4)*k)];
E1 = coE*[sum((u-u1).^2); sum(BH*(u.^2))];
E2 = coE*[sum((u1-u2).^2); sum(BH*(u1.^2))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OR

coE = [.5*(h^2)/k^2 , .5*(kappa^2)*(1/h^2)];
Energy = zeros(Nf,1);
Energy(n) = coE*[sum((u1-u2).^2); (LA*u1)' * (LA*u2)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
