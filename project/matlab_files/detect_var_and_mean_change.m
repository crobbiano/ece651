%% more testing
clear all
clc

f0 = .1;
n=1:1000;
% signal_snr = 20;
vars = .1;
signal_snr = -10*log10(vars)

% true_signal = cos(f0.*n);
true_signal = gen_ux0(length(n), 2);
noisy_true_signal = awgn(true_signal, signal_snr);
% truevar = var(noisy_true_signal-true_signal);
truevar = vars;

var_offset = .0;
mean_offset = .3;
mod_snr = -10*log10(truevar + var_offset)
mod_signal = true_signal;
mod_signal(1:99) = awgn(true_signal(1:99), signal_snr);
mod_signal(100:199) = awgn(true_signal(100:199), mod_snr);
mod_signal(200:end) = awgn(true_signal(200:end), mod_snr) + mean_offset;
% mod_signal(200:end) = awgn(true_signal(200:end), mod_snr) + linspace(0,mean_offset,length(n) - 200 + 1);
x = mod_signal;
modvar = truevar + var_offset
change_points = [100 200];

alpha = [-10:.01:10]';
mean_thresh = 3*sqrt(truevar);

% beta = [-truevar+.0001:.001:2*truevar]';
beta = [0:.001:2*truevar]';
var_thresh = (truevar)*1.5;
% var_thresh2 = truevar*8;
var_thresh2 = sqrt(truevar)*2;

alarm_window_sz = 30; % size of window to count alarm detections in
alarm_window_sz2 = 20; % size of window to count alarm detections in
alarm_time = nan;
alarm_time_m = nan;
alarm_window_idx = 1;
alarm_window_idx_m = 1;
alarms = [];
alarms_m = [];
final_time = -1;

found_var = 0;
found_mu = 0;

%% Testing for deviation form the known mean function
% xlikelihood = nan(length(beta), 1);
xlikelihood = nan( length(beta), 1);
xlikelihood_m = nan( length(alpha), 1);


for t=n % This is the time index
    % At each step in time, find the maximal likelihood of the mea
    % offset parameter, alpha.  If the difference between the observed
    % signal and the true mean function is > the threshold then tag that
    % time index as the alarm time,t_a.  Repeat this process, and if there are
    % alarm_window_sz discrepencies in a row then notify that a change in
    % the distribution has occured
    
    if ( found_var || found_mu)
        display(['Change in the distribution at time t = ' num2str(max([alarm_time alarm_time_m])) '; delay time = ' num2str(final_time)])
        break;
    end
     
    
        % calculate the likelihood of x(n) with different mean offsets
    for a=1:length(alpha);
%         xlikelihood_m(a) = (1/sqrt(2*pi*truevar))*exp(-(1/(2*truevar))*(x(t) + alpha(a) - true_signal(t))^2);
        xlikelihood_m(a) = (1/sqrt(2*pi*truevar))*exp(-(1/(2*truevar))*(x(t) - alpha(a))^2);
    end
    
    % find the mean offset that produced the greatest likelihood
    [rowmax rowmaxidx] = max(xlikelihood_m);
    [colmax colmaxidx] = max(max(xlikelihood_m));
    detected_mu = -alpha(rowmaxidx);
    
    % compare the mean offset to the mean threshold
    if (abs(detected_mu) > mean_thresh )
%     if (abs(detected_mu) > 1.5*abs(true_signal(t)) )
        if (t ~= alarm_time_m + alarm_window_idx_m)
            alarm_window_idx_m = 1;
            alarm_time_m = t;
        else
            alarm_window_idx_m = alarm_window_idx_m + 1;
        end
        alarms_m(end+1) = t;
    end
    % count the number of alarms in the past alarm_window_sz and flag
    % if more than the threshold
    alarm_thresh_m = alarm_window_sz/2;
    num_alarms_m = sum(alarms_m > t - alarm_window_sz);
    if (num_alarms_m >= alarm_thresh_m)
        final_time = alarm_time_m - change_points(2);
        found_mu = 1;
        display(['Found mean of ' num2str(detected_mu) '; true = ' num2str(true_signal(t))])
    end
    
    % calculate the likelihood of x(n) with different mean offsets
    for b=1:length(beta);
%         xlikelihood(b) = (1/sqrt(2*pi*(truevar + beta(b))))*exp(-(1/(2*(truevar + beta(b))))*(x(t) - true_signal(t))^2);
        xlikelihood(b) = (1/sqrt(2*pi*(beta(b))))*exp(-(1/(2*(beta(b))))*(x(t) - true_signal(t))^2);
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
            final_time = alarm_time - change_points(1);
            display(['Found thresh'])
        end
    end
    
    if (t < alarm_time + alarm_window_sz2 && alarm_time - alarm_window_sz2 > 1)
        % track the distrbution of the alarms from the current alarm window
        offset = floor(alarm_window_sz2); % look into the past if the info is available
        windowvar = var(x(alarm_time-offset:t)-true_signal(alarm_time-offset:t));
        if (windowvar > var_thresh2)
            found_var = 1;
            final_time = alarm_time - change_points(1);
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
plot(n,x,'-',n,true_signal,'--', change_points, x(change_points), 'ko', alarms, zeros(length(alarms)),'xr',alarms_m, zeros(length(alarms_m)),'ob')
plot(n,true_signal+sqrt(var_thresh),'--.k',n,true_signal-sqrt(var_thresh),'--.k')
line([t t],ylim,'Color',[1,0,0])
for i=1:length(change_points)
    line([change_points(i) change_points(i)],ylim,'Color',[0,0,0],'LineStyle','--')
end
xlim([0 t+10])
grid minor