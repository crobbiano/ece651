%% more testing
clear all
clc

%% Testing for deviation form the known mean function
% mean_thresh_mult = 0:.1:3;
mean_thresh_mult = 2;
for thresh_idx=1:length(mean_thresh_mult)
    thresh_idx
    found_mu = 0;
    for test_idx = 1:1
        f0 = .02;
        n=1:1000;
        vars = .1;
        signal_snr = -10*log10(vars);
        
        true_signal = gen_ux0(length(n),3);
        noisy_true_signal = awgn(true_signal, signal_snr);
        truevar = var(noisy_true_signal-true_signal);
        
        mean_offset = .4*max(true_signal);
        mod_signal = true_signal;
        mod_signal(760:1000) = true_signal(760:1000) + linspace(0,mean_offset,1000 -760 + 1);
        mod_signal(760:1000) = true_signal(760:1000) + mean_offset;
        mod_signal(1001:end) = true_signal(1001:end) + mean_offset;
        x = awgn(mod_signal, signal_snr);
        change_points = [760 ];
        
        alpha = [-10:.01:10]';
        mean_thresh = mean_thresh_mult(thresh_idx)*sqrt(truevar);
        
        
        alarm_window_sz = 20; % size of window to count alarm detections in
        alarm_time = nan;
        alarm_window_idx = 1;
        alarms = [];
        final_time = -1;
        
        xlikelihood = nan(length(alpha), 1);
        % xlikelihood = nan(length(alpha), length(n));
        for t=n % This is the time index
            % At each step in time, find the maximal likelihood of the mea
            % offset parameter, alpha.  If the difference between the observed
            % signal and the true mean function is > the threshold then tag that
            % time index as the alarm time,t_a.  Repeat this process, and if there are
            % alarm_window_sz discrepencies in a row then notify that a change in
            % the distribution has occured
            
            if (found_mu)
                display(['Change in mean of the distribution at time t = ' num2str(alarm_time) '; delay time = ' num2str(final_time)])
                break;
            end
            
            % calculate the likelihood of x(n) with different mean offsets
            for a=1:length(alpha);
                %                 xlikelihood(a,1) = (1/sqrt(2*pi*truevar))*exp(-(1/(2*truevar))*(x(t) + alpha(a) - true_signal(t))^2);
                xlikelihood(a,1) = (1/sqrt(2*pi*truevar))*exp(-(1/(2*truevar))*(x(t) - alpha(a))^2);
            end
            
            % find the mean offset that produced the greatest likelihood
            [rowmax rowmaxidx] = max(xlikelihood);
            [colmax colmaxidx] = max(max(xlikelihood));
            detected_mu = alpha(rowmaxidx(colmaxidx))
            
            llr = xlikelihood(rowmaxidx(colmaxidx))/((1/sqrt(2*pi*truevar))*exp(-(1/(2*truevar))*(x(t) - true_signal(t))^2))
            
            
            delta_mu = true_signal(t) - detected_mu
            % compare the mean offset to the mean threshold
            %             if (abs(delta_mu) > 2*sqrt(truevar))
            %             if (abs(detected_mu) > 6*abs(true_signal(t)) )
            if (llr > 10 )
                if (t ~= alarm_time + alarm_window_idx)
                    alarm_window_idx = 1;
                    alarm_time = t;
                else
                    alarm_window_idx = alarm_window_idx + 1;
                end
                alarms(end+1) = t;
            end
            alarm_thresh_m = alarm_window_sz/4;
            num_alarms_m = sum(alarms > t - alarm_window_sz);
            if (num_alarms_m >= alarm_thresh_m)
                final_time = alarm_time - change_points(1);
                found_mu = 1;
                display(['Found mean of ' num2str(detected_mu) '; true = ' num2str(true_signal(t))])
            end
        end
        
    end
    alarm_times(test_idx) = final_time;
    
    %% Calculate some statistics
    Pfa = sum(alarm_times<0)/length(alarm_times)
    Pd = (length(alarm_times)-sum(alarm_times==-1))/length(alarm_times)
    alarm_avg = mean(alarm_times(alarm_times>0))
    alarm_var = var(alarm_times(alarm_times>0))
    
    data(thresh_idx).Pfa = Pfa;
    data(thresh_idx).Pd = Pd;
    data(thresh_idx).alarm_avg = alarm_avg;
    data(thresh_idx).alarm_var = alarm_var;
    data(thresh_idx).mean_thresh = mean_thresh;
    data(thresh_idx).signal_var = vars;
    data(thresh_idx).signal_snr = signal_snr;
    data(thresh_idx).mean_offset = mean_offset;
    
end
% save('mean_data.mat','data');
%% plot things
figure(10);clf;hold on
plot(n,x,'-',n,true_signal,'--', change_points, x(change_points), 'ko', alarms, zeros(length(alarms)),'xr')
plot(n,true_signal+sqrt(truevar),'--.k',n,true_signal-sqrt(truevar),'--.k')
line([t t],ylim,'Color',[1,0,0])
for i=1:length(change_points)
    line([change_points(i) change_points(i)],ylim,'Color',[0,0,0],'LineStyle','--')
end
%% plot density
% x = 0:200;
% figure(398);clf;
% plot(x, normpdf(x,alarm_avg, sqrt(alarm_var)))