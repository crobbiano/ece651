%% more testing
clear all
clc

f0 = .02;
n=1:1500;
% signal_snr = 20;
vars = 1;
signal_snr = -10*log10(vars)

true_signal = cos(f0.*n);
true_signal = gen_ux0(length(n),2);
noisy_true_signal = awgn(true_signal, signal_snr);
truevar = var(noisy_true_signal-true_signal);

var_offset = .5;
mod_snr = -10*log10(vars + var_offset)
mod_signal = true_signal;
mod_signal(1:250) = awgn(true_signal(1:250), signal_snr);
% mod_signal(40:89) = awgn(true_signal(40:89), signal_snr)+3;
mod_signal(251:end) = awgn(true_signal(251:end), mod_snr);
x = mod_signal;
modvar = var(true_signal(251:350) - mod_signal(251:350))
change_points = [251 ];

beta = [-truevar+.000001:.001:2*truevar]';
var_thresh = (truevar)/1;
var_thresh2 = truevar*8;

alarm_window_sz = 20; % size of window to count alarm detections in
alarm_window_sz2 = 20; % size of window to count alarm detections in
alarm_time = nan;
alarm_window_idx = 1;
alarms = [];

found_var = 0;

%% Testing for deviation form the known mean function
xlikelihood = nan(length(beta), 1);
% xlikelihood = nan(length(alpha), length(n));


for t=n % This is the time index
    % At each step in time, find the maximal likelihood of the mea
    % offset parameter, alpha.  If the difference between the observed
    % signal and the true mean function is > the threshold then tag that
    % time index as the alarm time,t_a.  Repeat this process, and if there are
    % alarm_window_sz discrepencies in a row then notify that a change in
    % the distribution has occured
    
    if (alarm_window_idx >= alarm_window_sz || found_var)
        display(['Change in variance of the distribution at time t = ' num2str(alarm_time) '; delay time = ' num2str(alarm_time - change_points(1))])
        break;
    end
    
    % calculate the likelihood of x(n) with different mean offsets
    for a=1:length(beta);
        xlikelihood(a,1) = (1/sqrt(2*pi*(truevar + beta(a))))*exp(-(1/(2*(truevar + beta(a))))*(x(t) - true_signal(t))^2);
    end
    
    % find the mean offset that produced the greatest likelihood
    [rowmax rowmaxidx] = max(xlikelihood);
    detected_var = beta(rowmaxidx);
    
    % compare the mean offset to the mean threshold
    if ((abs(detected_var)) > var_thresh )
        %     if ((abs(detected_var)) > var_thresh  )
        if (t ~= alarm_time + alarm_window_idx)
            alarm_window_idx = 1;
            alarm_time = t;
        else
            alarm_window_idx = alarm_window_idx + 1;
        end
        alarms(end+1) = t;
        
        % count the number of alarms in the past alarm_window_sz and flag
        % if more than the threshold
        alarm_thresh = alarm_window_sz/2;
        num_alarms = sum(alarms > t - alarm_window_sz);
        if (num_alarms >= alarm_thresh)
            found_var = 1;
            display(['Found thresh'])
        end
    end
    
    if (t < alarm_time + alarm_window_sz2 && alarm_time - alarm_window_sz2 > 1)
        % track the distrbution of the alarms from the current alarm window
        offset = floor(alarm_window_sz2); % look into the past if the info is available
        windowvar = var(x(alarm_time-offset:t)-true_signal(alarm_time-offset:t));
        if (windowvar > var_thresh2)
            found_var = 1;
            display(['Found window'])
        end
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