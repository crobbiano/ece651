%% more testing
clear all
clc

f0 = .02;
n=1:500;
vars = .1;
signal_snr = -10*log10(vars)

true_signal = gen_ux0(length(n));
noisy_true_signal = awgn(true_signal, signal_snr);
truevar = var(noisy_true_signal-true_signal);

mean_offset = 1;
mod_signal = true_signal;
mod_signal(180:260) = true_signal(180:260) + linspace(0,mean_offset,260-180 + 1);
mod_signal(261:end) = true_signal(261:end) + mean_offset;
x = awgn(mod_signal, signal_snr);
change_points = [180 ];

alpha = [-10:.01:10]';
mean_thresh = 1.5*sqrt(truevar);

alarm_window_sz = 5; % size of window to count alarm detections in
alarm_time = nan;
alarm_window_idx = 1;
alarms = [];

%% Testing for deviation form the known mean function
xlikelihood = nan(length(alpha), 1);
% xlikelihood = nan(length(alpha), length(n));


for t=n % This is the time index
    % At each step in time, find the maximal likelihood of the mea
    % offset parameter, alpha.  If the difference between the observed
    % signal and the true mean function is > the threshold then tag that
    % time index as the alarm time,t_a.  Repeat this process, and if there are
    % alarm_window_sz discrepencies in a row then notify that a change in
    % the distribution has occured
    
    if (alarm_window_idx >= alarm_window_sz)
        display(['Change in mean of the distribution at time t = ' num2str(alarm_time) '; delay time = ' num2str(alarm_time - change_points(1))])
        break;
    end
    
    % calculate the likelihood of x(n) with different mean offsets
    for a=1:length(alpha);
        xlikelihood(a,1) = (1/sqrt(2*pi*truevar))*exp(-(1/(2*truevar))*(x(t) + alpha(a) - true_signal(t))^2);
    end
    
    % find the mean offset that produced the greatest likelihood
    [rowmax rowmaxidx] = max(xlikelihood);
    [colmax colmaxidx] = max(max(xlikelihood));
    detected_mu = -alpha(rowmaxidx(colmaxidx));
    
    % compare the mean offset to the mean threshold
    if (abs(detected_mu) > mean_thresh )
        if (t ~= alarm_time + alarm_window_idx)
            alarm_window_idx = 1;
            alarm_time = t;
        else
            alarm_window_idx = alarm_window_idx + 1;
        end
        alarms(end+1) = t;
    end
end



%
% [rowmax rowmaxidx] = max(xmodlikelihood);
% [colmax colmaxidx] = max(max(xmodlikelihood));
% fake = -alpha(rowmaxidx(colmaxidx))

% Compare the recovered mean to the true mean and report if we are outside
% of the bounds

%% plot things
figure(10);clf;hold on
plot(n,x,'-',n,true_signal,'--', change_points, x(change_points), 'ko', alarms, zeros(length(alarms)),'xr')
line([t t],ylim,'Color',[1,0,0])
for i=1:length(change_points)
    line([change_points(i) change_points(i)],ylim,'Color',[0,0,0],'LineStyle','--')
end