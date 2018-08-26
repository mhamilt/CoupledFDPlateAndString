% Analysis Module for Thin Plate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Modal Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculating plate modes

if modeAn
  p = [1:10];
  q = [1:10]';

  [P,Q] = meshgrid(p,q);

    mfreqs = (pi*kappa/(2))*((P.^2)*Lx/Ly + (Q.^2)*Ly/Lx);     % mode frequencies
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
font_s = 14; % font point size
set(0,'DefaultAxesFontSize',font_s);

% Output FFT

if modeAn
  YF = abs(fft(y));
  fy = linspace(0, SR, length(y)); % frequency axis
  figure(1);
  plot(fy, 20*log10(YF));
  ax1 = gca;
  line([mfreqs(:) mfreqs(:)],ax1.YLim,'Color',[1 0 0])
  ax1.XLim = [0, 1500];
  title('Scheme Output FFT', 'FontSize', font_s)
  xlabel('Freq._{Hz}');
  ylabel('magnitude');
  legend('FFT Output_{dB}','Mode Frequencies_{Hz}')
  set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen
end

if energyAn
  %change in energy
  dEdt = (Energy(2:end)-Energy(1:end-1));
  figure(2)
  subplot(2,1,1)
  plot(Energy);
  hold on
  plot(dEdt)
  plot(KE)
  plot(PE)
  ax2 = gca;
  title('Change in Energy at Time Step k', 'FontSize', font_s)
  xlabel('Sample_{k}');
  ylabel('Power_{Not Watts}');
  set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen
  legend('Total Energy','Power_{dE/dt}','Kinetic','Potential')


  subplot(2,1,2)
  plot(dEdt,'.')
  title('dE/dt rounding error', 'FontSize', font_s)
  xlabel('Sample_{k}');
  ylabel('dE/dt');
end

% if plot_on
%   v = VideoWriter('2DWave','Motion JPEG AVI');
%   % v.CompressionRatio = 3;
%   v.FrameRate = 24;
%   open(v)
%   writeVideo(v,F(1:4000))
%   close(v)
% end
