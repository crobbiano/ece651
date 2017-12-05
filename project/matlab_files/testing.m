
f0 = .1;
n=1:2500;
% signal_snr = 20;
vars = .1;
signal_snr = -10*log10(vars)

% true_signal = cos(f0.*n);
true_signal = gen_ux0(length(n), 8);
noisy_true_signal = awgn(true_signal, signal_snr);
% truevar = var(noisy_true_signal-true_signal);
truevar = vars;

var_offset = .05;
mean_offset = .5;
mod_snr = -10*log10(truevar + var_offset)
mod_signal = true_signal;
mod_signal(1:99) = awgn(true_signal(1:99), signal_snr);
mod_signal(100:499) = awgn(true_signal(100:499), mod_snr);
mod_signal(500:end) = awgn(true_signal(500:end), mod_snr) + mean_offset;
% mod_signal(200:end) = awgn(true_signal(200:end), mod_snr) + linspace(0,mean_offset,length(n) - 200 + 1);
x = mod_signal;
change_points = [100 500];

h = 5.5;
[kk,cc]=CUSUM(mod_signal',true_signal',sqrt(truevar),h);
shew3=Shewhart(mod_signal',true_signal',sqrt(truevar)); % Shewhart.m
[SRp,Wp]=SRpoi(mod_signal',true_signal',sqrt(truevar)); % SRpoi.m
[SRn,Wn]=SRnorm(mod_signal',true_signal',sqrt(truevar)); % SRnorm.m


xlen=1:length(mod_signal);
% x=x.*100./60.*v;
disp(['The Shewart method responded after ',num2str(shew3),' measurements.']);
disp(['The CUSUM method responded after ',num2str(kk),' measurements.']);
disp(['The Poisson SR method responded after ',num2str(SRp),' measurements.']);
disp(['The normal SR method responded after ',num2str(SRn),' measurements.']);

if (kk~=length(n))
    display(['Found change at t=' num2str(kk)])
else
    display(['!! No change detected !!'])
end



% Subplot
xlimit = max([kk, SRn, shew3]);

fig=figure(99);clf;
set(0, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex')


% Plot signal with standard deviation lines
fig4 = subplot(4,1,1); hold on
plot(fig4, n,x,'-',n,true_signal,'--', change_points, x(change_points), 'ko')
plot(fig4, n,true_signal+sqrt(truevar),'--k',n,true_signal-sqrt(truevar),'--k')
line([kk kk],ylim,'Color',[1,0,0])
line([SRn SRn],ylim,'Color',[0,1,0])
line([shew3 shew3],ylim,'Color',[0,0,1])
for i=1:length(change_points)
    line([change_points(i) change_points(i)],ylim,'Color',[0,0,0],'LineStyle','--')
end
xlim([0 xlimit+10])
grid minor
xlabel('Samples')
ylabel('Signal Amplitude')
title(['SNR=' num2str(signal_snr) '; 1st Change (n=100): SNR=' num2str(mod_snr) '; 2nd Change (n=500): DC mean shift= ' num2str(mean_offset)])

% Plot CUSUM statistic
fig2 = subplot(4,1,3);
plot(fig2, xlen,cc,'black')
hline=refline(0,h);
set(hline,'LineStyle',':','Color','black');
title('CUSUM')
xlabel('Samples');
ylabel('CUSUM Statistic');
line([kk kk],ylim,'Color',[1,0,0])
ylim([0, 10])
xlim([0 xlimit+10])
for i=1:length(change_points)
    line([change_points(i) change_points(i)],ylim,'Color',[0,0,0],'LineStyle','--')
end

% Plot SR Gauss statistic
fig3 = subplot(4,1,2);
plot(fig3, xlen,Wn,'black')
% hline=refline(0,h);
% set(hline,'LineStyle',':','Color','black');
title('SR Gauss')
xlabel('Samples');
ylabel('SR Gauss Statistic');
line([SRn SRn],ylim,'Color',[0,1,0])
ylim([0, Wn(SRn)+.2*Wn(SRn)])
xlim([0 xlimit+10])
for i=1:length(change_points)
    line([change_points(i) change_points(i)],ylim,'Color',[0,0,0],'LineStyle','--')
end

% Plot Shewart statistic
fig5 = subplot(4,1,3);
plot(fig5, xlen,x,'b',xlen,true_signal+3*sqrt(truevar),'k--')
% hline=refline(0,h);
% set(hline,'LineStyle',':','Color','black');
title('Shewart 3-sigma')
xlabel('Samples');
ylabel('Signal Amplitude');
line([shew3 shew3],ylim,'Color',[0,0,1])
xlim([0 xlimit+10])
for i=1:length(change_points)
    line([change_points(i) change_points(i)],ylim,'Color',[0,0,0],'LineStyle','--')
end

% Plot CUSUM statistic
fig2 = subplot(4,1,4);
plot(fig2, xlen,cc,'black')
hline=refline(0,h);
set(hline,'LineStyle',':','Color','black');
title('CUSUM')
xlabel('Samples');
ylabel('CUSUM Statistic');
line([kk kk],ylim,'Color',[1,0,0])
ylim([0, 10])
xlim([0 xlimit+10])
for i=1:length(change_points)
    line([change_points(i) change_points(i)],ylim,'Color',[0,0,0],'LineStyle','--')
end
