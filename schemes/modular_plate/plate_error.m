% Error Checking Module for plate FD Scheme
paraLs = {'SR'};
paraCons = [(SR<1e3),];

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
