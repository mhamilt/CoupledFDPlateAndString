% Analysis Module for Thin Plate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Modal Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculating plate modes

if modeAn
  p = [1:5];
  q = [1:5]';

  [P,Q] = meshgrid(p,q);
  %%NOTE double check calculation
  mfreqs = (pi*kappa/(2))*((P.^2/Ly^2) + (Q.^2/Lx^2));     % mode frequencies
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

  if bctype == 1
    line([mfreqs(:) mfreqs(:)],ax1.YLim,'Color',[1 0 0])
  end
  if bctype == 2
    mfreqs = kappa*[35.10, 72.90, 107.47, 131.63, 164.39, 210.35, 219.32, 242.20, 295.69, 370.66]/(2*pi);
    line([mfreqs(:) mfreqs(:)],ax1.YLim,'Color',[1 0 0])
  end
  ax1.XLim = [0, 1500];

  title('Scheme Output FFT', 'FontSize', font_s)
  xlabel('Freq._{Hz}');
  ylabel('magnitude');
  legend('FFT Output_{dB}','Mode Frequencies_{Hz}')
  set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen
end

if energyAn
  %change in energy

  dEdt = (Energy(2:end-1)-Energy(1:end-2))/k;
  figure(2)
  subplot(2,1,1)
  plot(Energy);
  hold on
  plot(dEdt)
  plot(KE)
  plot(PE)
  % plot(Loss)
  ax2 = gca;
  title('Change in Energy at Time Step k', 'FontSize', font_s)
  xlabel('Sample_{k}');
  ylabel('Power_{Not Watts}');
  set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); % fullscreen
  legend('Total Energy','Power_{dE/dt}','Kinetic','Potential')


  subplot(2,1,2)
  plot(dEdt,'.')
  hold on
  % plot(Loss(2:end-1), '.')
  title('dE/dt rounding error', 'FontSize', font_s)
  xlabel('Sample_{k}');
  ylabel('dE/dt');
  % plot(dEdt+Loss(2:end-1),'.')
end

%%% generate spy graph and demarcate matrix points
figure(3)
spy(B,'k')
spy_ax = gca;
line([ss+.5 ss+.5],spy_ax.YLim,'Color',[1 0 0],'LineWidth',1)
line(spy_ax.XLim,[ss+.5 ss+.5],'Color',[1 0 0],'LineWidth',1)
% if plot_on
%   v = VideoWriter('2DWave','Motion JPEG AVI');
%   % v.CompressionRatio = 3;
%   v.FrameRate = 24;
%   open(v)
%   writeVideo(v,F(1:4000))
%   close(v)
% end
