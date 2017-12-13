
clear all
clc
%%

f0 = 1;
n=1:1000;
% signal_snr = 20;
truevar = .2;
signal_snr = -10*log10(truevar)

% true_signal = cos(f0.*n);
true_signal = gen_ux0(length(n), 30)+5;
% noisy_true_signal = awgn(true_signal, signal_snr);
% truevar = var(noisy_true_signal-true_signal);

var_offset = .0;
mean_offset = -2;
mod_snr = -10*log10(truevar + var_offset)
mod_signal = true_signal;
mod_signal(1:199) = awgn(true_signal(1:199), signal_snr);
mod_signal(200:499) = awgn(true_signal(200:499), mod_snr);
mod_signal(500:end) = awgn(true_signal(500:end), mod_snr) + mean_offset;
% mod_signal(200:end) = awgn(true_signal(200:end), mod_snr) + linspace(0,mean_offset,length(n) - 200 + 1);
x = mod_signal;
change_points = [200 500];

h = 5.5;
[kk,cc]=CUSUM(mod_signal',true_signal',sqrt(truevar),h);
cu_thresh = -sqrt(truevar)*5;
cu_win_sz = 20;
[cu_var_time,cu_var_cc]=CUSUM_var(mod_signal',true_signal',sqrt(truevar),cu_thresh,cu_win_sz);
shew3=Shewhart(mod_signal',true_signal',sqrt(truevar)); % Shewhart.m
Wstop = 700;
[SRn,Wn]=SRnorm(mod_signal',true_signal',sqrt(truevar),Wstop); % SRnorm.m

%% Need to put this into its own function
% cusum adopted for both var and mean
data = mod_signal;
m=100;
s=truevar;
a=true_signal;
t=length(n);
for i=1:t
    if (i<=m)
        
        sig_est = (1/i)*sum((data(1:i)-a(1:i)).^2);
        mu_est = (1/i)*sum(data(1:i)-a(1:i));
        c(i) = .5*log(s/sig_est) + .5*(1/sig_est - 1/s)*sum((data(1:i)-a(1:i)).^2) + (1/sig_est)*sum(((data(1:i)-a(1:i))*mu_est) - mu_est^2);
    else
        sig_est = (1/m)*sum((data(i-m:i)-a(i-m:i)).^2);
        mu_est = (1/m)*sum(data(i-m:i)-a(i-m:i));
        c(i) = .5*log(s/sig_est) + .5*(1/sig_est - 1/s)*sum((data(i-m:i)-a(i-m:i)).^2) + (1/sig_est)*sum(((data(i-m:i)-a(i-m:i))*mu_est) - mu_est^2);
    end
%     cmin(i) = min([0,c(i)]);
    cmin(i) = c(i);
end

if (sum(cmin<-22)>0)
    num1 = find(cmin < -22);
    num1 = num1(1);
else
    num1=length(n)
end

%%

xlen=1:length(mod_signal);
% x=x.*100./60.*v;
disp(['The Shewart method responded after ',num2str(shew3),' measurements.']);
disp(['The CUSUM method responded after ',num2str(kk),' measurements.']);
disp(['The CUSUM adopted method responded after ',num2str(cu_var_time),' measurements.']);
disp(['The CUSUM adopted2 method responded after ',num2str(num1),' measurements.']);
disp(['The normal SR method responded after ',num2str(SRn),' measurements.']);

if (kk~=length(n) )
    display(['Found change at t=' num2str(kk)])
else
    display(['!! No change detected !!'])
end



% Subplot
xlimit = max([kk, SRn, shew3, cu_var_time, num1]);

fig=figure(99);clf;
set(0, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex')


% Plot signal with standard deviation lines
fig4 = subplot(5,1,1); hold on
plot(fig4, n,x,'b-',n,true_signal,'r--', change_points, x(change_points), 'ko')
% plot(fig4, n,true_signal+sqrt(truevar),'--k',n,true_signal-sqrt(truevar),'--k')

xlim([0 xlimit+10])
ylim([2.5 7.5])
line([kk kk],ylim,'Color',[1,0,0])
% line([cu_var_time cu_var_time],ylim,'Color',[.7,.5,.6])
line([SRn SRn],ylim,'Color',[0,1,0])
line([shew3 shew3],ylim,'Color',[0,1,1])
line([num1 num1],ylim,'Color',[.2,.7,.7])
for i=1:length(change_points)
    line([change_points(i) change_points(i)],ylim,'Color',[0,0,0],'LineStyle','--')
end
grid minor
% xlabel('Samples')
ylabel('Signal Amplitude')
title(['SNR=' num2str(signal_snr) '; 1st Change (n=100): SNR=' num2str(mod_snr) '; 2nd Change (n=500): DC mean shift= ' num2str(mean_offset)])

% Plot SR Gauss statistic
fig3 = subplot(5,1,2);
plot(fig3, xlen,Wn,'black')
% hline=refline(0,h);
% set(hline,'LineStyle',':','Color','black');
title('SR Gauss')
% xlabel('Samples');
ylabel('SR Gauss Statistic');
line([SRn SRn],ylim,'Color',[0,1,0])
Wstopline=refline(0,Wstop);
set(Wstopline,'LineStyle',':','Color','black');
ylim([0, Wn(SRn)+.2*Wn(SRn)])
xlim([0 xlimit+10])
for i=1:length(change_points)
    line([change_points(i) change_points(i)],ylim,'Color',[0,0,0],'LineStyle','--')
end

% Plot Shewart statistic
fig5 = subplot(5,1,3);
plot(fig5, xlen,x,'b',xlen,true_signal+3*sqrt(truevar),'k--')
% hline=refline(0,h);
% set(hline,'LineStyle',':','Color','black');
title('Shewart 3-sigma')
% xlabel('Samples');
ylabel('Signal Amplitude');
line([shew3 shew3],ylim,'Color',[0,1,1])
xlim([0 xlimit+10])
for i=1:length(change_points)
    line([change_points(i) change_points(i)],ylim,'Color',[0,0,0],'LineStyle','--')
end

% Plot CUSUM statistic
fig2 = subplot(5,1,4);
plot(fig2, xlen,cc,'black')
hline=refline(0,h);
set(hline,'LineStyle',':','Color','black');
title('CUSUM')
% xlabel('Samples');
ylabel('CUSUM Statistic');
line([kk kk],ylim,'Color',[1,0,0])
ylim([-.5, 10])
xlim([0 xlimit+10])
for i=1:length(change_points)
    line([change_points(i) change_points(i)],ylim,'Color',[0,0,0],'LineStyle','--')
end

%
% subplot(5,1,5)
% plot(cu_var_cc,'b')
% title('CUSUM adopted for change in variance')
% xlabel('Samples');
% ylabel('CUSUM Statistic');
% ylim([cu_thresh*2 - 1, .5])
% xlim([0 xlimit+10])
% for i=1:length(change_points)
%     line([change_points(i) change_points(i)],ylim,'Color',[0,0,0],'LineStyle','--')
% end
% hline=refline(0,cu_thresh);
% set(hline,'LineStyle',':','Color','black');
% line([cu_var_time cu_var_time],ylim,'Color',[.7,.5,.6])


subplot(5,1,5)
plot(cmin,'k')
title('SMAGLR')
xlabel('Samples');
ylabel('SMAGLR Statistic');
for i=1:length(change_points)
    line([change_points(i) change_points(i)],ylim,'Color',[0,0,0],'LineStyle','--')
end
xlim([0 xlimit+10])
ylim([-50 0])
hline=refline(0,-20);
set(hline,'LineStyle',':','Color','black');
line([num1 num1],ylim,'Color',[.2,.7,.7])

%% more testing of new code
