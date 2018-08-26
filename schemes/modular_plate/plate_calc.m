%% Main calculation for FD thin plate

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
