%% gen_results
% generate results based on the comparison_of_methods


clear all
clc
%%

f0 = .1;
n=1:1000;

var_list = [0.001 0.01 0.1 0.2 0.3 0.5 0.7 0.9];
var_offset_list = [0 0.001 0.01 0.1 0.2 0.3 0.5 0.7 0.9];
mean_list = -3:.25:3;
numInnerLoops = 100;

% data = (truevar, var_offset, mean_offset, signal_snr, mod_snr, kk, shew3, SRn, num1)
data = zeros(length(var_list)*length(var_offset_list)*length(mean_list)*numInnerLoops, 9);
data_idx = 0;

gen_data = 0;
if (gen_data)
    for truevar = var_list;
        for var_offset = var_offset_list;
            for mean_offset = mean_list;
                for innerLoopNum = 1:numInnerLoops;
                    data_idx = data_idx + 1;
                    if (mod(data_idx,100) == 1)
                        display(['Entry ' num2str(data_idx)])
                    end
                    true_signal = gen_ux0(length(n), 8)+5;
                    
                    signal_snr = -10*log10(truevar);
                    mod_snr = -10*log10(truevar + var_offset);
                    change_points = [200 500];
                    
                    mod_signal = zeros(1, length(n));
                    mod_signal(1:change_points(1)-1) = awgn(true_signal(1:change_points(1)-1), signal_snr);
                    mod_signal(change_points(1):change_points(2)-1) = awgn(true_signal(change_points(1):change_points(2)-1), mod_snr);
                    mod_signal(change_points(2):end) = awgn(true_signal(change_points(2):end), mod_snr) + mean_offset;
                    
                    h = 5.5;
                    [kk,cc]=CUSUM(mod_signal',true_signal',sqrt(truevar),h);
                    
                    shew3=Shewhart(mod_signal',true_signal',sqrt(truevar)); % Shewhart.m
                    
                    Wstop = 700;
                    [SRn,Wn]=SRnorm(mod_signal',true_signal',sqrt(truevar),Wstop); % SRnorm.m
                    
                    %% Need to put this into its own function
                    % cusum adopted for both var and mean
                    m=100;
                    c = zeros(1,length(n));
                    cmin = zeros(1,length(n));
                    c_thresh = -20;
                    
                    for i=1:length(n);
                        if (i<=m)
                            c(i) = 0;
                            cmin(i) = 0;
                        else
                            sig_est = (1/m)*sum((mod_signal(i-m:i)-true_signal(i-m:i)).^2);
                            mu_est = (1/m)*sum(mod_signal(i-m:i)-true_signal(i-m:i));
                            c(i) = .5*log(truevar) - .5 * log(sig_est) + .5*(1/sig_est - 1/truevar)*sum((mod_signal(i-m:i)-true_signal(i-m:i)).^2);
                            cmin(i) = min([0,.5*log(truevar) - .5 * log(sig_est) + .5*(1/sig_est - 1/truevar)*sum((mod_signal(i-m:i)-true_signal(i-m:i)).^2)]);
                        end
                    end
                    
                    if (sum(cmin < c_thresh)>0)
                        num1 = find(cmin < c_thresh);
                        num1 = num1(1);
                    else
                        num1=length(n);
                    end
                    
                    % store some data
                    data(data_idx,:) = [truevar, var_offset, mean_offset, signal_snr, mod_snr, kk, shew3, SRn, num1];
                end
            end
        end
    end
end
%% save the data
% save('full_data.mat','data');
load('full_data.mat')
% load('full_set.mat')
%% process all the data
% data = (truevar, var_offset, mean_offset, signal_snr, mod_snr, kk, shew3, SRn, num1)

% easy false alarm probs
zerorows = find(and(any(data(:,[2])==0, 2), any(data(:,[3]) == 0, 2))); zerodata = data(zerorows,:);
zerovar = find(and(any(data(:,[2])==0, 2), any(data(:,[3]) ~= 0, 2))); zerovardata = data(zerovar,:);
zeromean = find(and(any(data(:,[2])~=0, 2), any(data(:,[3]) == 0, 2))); zeromeandata = data(zeromean,:);
for i=1:length(zeromeandata)
    zeromeandata(i,10) = zeromeandata(i,5) - (zeromeandata(i,4));
end

p_total = 0;
p_mean_total = 0;
p_var_total = 0;
pfa_mean_sr = 0; pfa_mean_cusum = 0; pfa_mean_shewart = 0; pfa_mean_modcu = 0;
pd_mean_sr = 0; pd_mean_cusum = 0; pd_mean_shewart = 0; pd_mean_modcu = 0;
pfa_var_sr = 0; pfa_var_cusum = 0; pfa_var_shewart = 0; pfa_var_modcu = 0;
pd_var_sr = 0; pd_var_cusum = 0; pd_var_shewart = 0; pd_var_modcu = 0;

add_sr = [];
add_cusum = [];
add_shewart = [];
add_modcu = [];

for idx = 1:size(zerodata,1)
    if (zerodata(idx,6)<length(n))
        pfa_var_cusum = pfa_var_cusum + 1;
    end
    if (zerodata(idx,7)<length(n))
        pfa_var_shewart = pfa_var_shewart + 1;
    end
    if (zerodata(idx,8)<length(n))
        pfa_var_sr = pfa_var_sr + 1;
    end
    if (zerodata(idx,9)<length(n))
        pfa_var_modcu = pfa_var_modcu + 1;
    end
    
    p_var_total = p_var_total + 1; % FIXME - this amount should only be added to total when calculating Pfa
end

% zer var data - delta mean
for idx = 1:size(zerovardata,1)
    if (zerovardata(idx,6)<length(n))
        pd_var_cusum = pd_var_cusum + 1;
        if ( zerovardata(idx,6)>=500)
            add_cusum(end + 1) = zerovardata(idx,6)-500;
        end
    end
    if (zerovardata(idx,6)<500)
        pfa_var_cusum = pfa_var_cusum + 1;
    end
    if ( zerovardata(idx,7)<length(n))
        pd_var_shewart = pd_var_shewart + 1;
        
        if ( zerovardata(idx,7)>=500)
            add_shewart(end + 1) = zerovardata(idx,7)-500;
        end
    end
    if (zerovardata(idx,7)<500)
        pfa_var_shewart = pfa_var_shewart + 1;
    end
    if ( zerovardata(idx,8)<length(n))
        pd_var_sr = pd_var_sr + 1;
        if ( zerovardata(idx,8)>=500)
            add_sr(end + 1) = zerovardata(idx,8)-500;
        end
    end
    if (zerovardata(idx,8)<500)
        pfa_var_sr = pfa_var_sr + 1;
    end
    if ( zerovardata(idx,9)<length(n))
        pd_var_modcu = pd_var_modcu + 1;
        
        if ( zerovardata(idx,9)>=500)
            add_modcu(end + 1) = zerovardata(idx,9)-500;
        end
    end
    if (zerovardata(idx,9)<500)
        pfa_var_modcu = pfa_var_modcu + 1;
    end
    
    p_var_total = p_var_total + 1;
end

add_cusum_var = sum(add_cusum)/numel(add_cusum);
add_shewart_var = sum(add_shewart)/numel(add_shewart);
add_sr_var = sum(add_sr)/numel(add_sr);
add_modcu_var = sum(add_modcu)/numel(add_modcu);

% zero mean data - delta var
for idx = 1:size(zeromeandata,1)
    if ( zeromeandata(idx,6)<length(n))
        pd_mean_cusum = pd_mean_cusum + 1;
        if ( zerovardata(idx,6)>=200)
            add_cusum(end + 1) = zerovardata(idx,6)-200;
        end
    end
    if (zeromeandata(idx,6)<200)
        pfa_mean_cusum = pfa_mean_cusum + 1;
    end
    if ( zeromeandata(idx,7)<length(n))
        
        pd_mean_shewart = pd_mean_shewart + 1;
        if ( zerovardata(idx,7)>=200)
            add_shewart(end + 1) = zerovardata(idx,7)-200;
        end
    end
    if (zeromeandata(idx,7)<200)
        pfa_mean_shewart = pfa_mean_shewart + 1;
    end
    if (zeromeandata(idx,8)<length(n))
        pd_mean_sr = pd_mean_sr + 1;
        if ( zerovardata(idx,8)>=200)
            add_sr(end + 1) = zerovardata(idx,8)-200;
        end
    end
    if (zeromeandata(idx,8)<200)
        pfa_mean_sr = pfa_mean_sr + 1;
    end
    if ( zeromeandata(idx,9)<length(n))
        pd_mean_modcu = pd_mean_modcu + 1;
        if ( zerovardata(idx,9)>=200)
            add_modcu(end + 1) = zerovardata(idx,9)-200;
        end
    end
    if (zeromeandata(idx,9)<200)
        pfa_mean_modcu = pfa_mean_modcu + 1;
    end
    
    p_mean_total = p_mean_total + 1;
end
add_cusum_mean = sum(add_cusum)/numel(add_cusum);
add_shewart_mean = sum(add_shewart)/numel(add_shewart);
add_sr_mean = sum(add_sr)/numel(add_sr);
add_modcu_mean = sum(add_modcu)/numel(add_modcu);

add_cusum = (add_cusum_mean + add_cusum_var)/2
add_shewart = (add_shewart_mean + add_shewart_var)/2
add_sr = (add_sr_mean + add_sr_var)/2
add_modcu = (add_modcu_mean + add_modcu_var)/2


pfa_sr = (pfa_mean_sr/p_mean_total + pfa_var_sr/p_var_total)/2
pfa_cusum = (pfa_mean_cusum/p_mean_total + pfa_var_cusum/p_var_total)/2
pfa_shewart = (pfa_mean_shewart/p_mean_total + pfa_var_shewart/p_var_total)/2
pfa_modcu = (pfa_mean_modcu/p_mean_total + pfa_var_modcu/p_var_total)/2


pd_sr = (pd_mean_sr/p_mean_total + pd_var_sr/p_var_total)/2
pd_cusum = (pd_mean_cusum/p_mean_total + pd_var_cusum/p_var_total)/2
pd_shewart = (pd_mean_shewart/p_mean_total + pd_var_shewart/p_var_total)/2
pd_modcu = (pd_mean_modcu/p_mean_total + pd_var_modcu/p_var_total)/2

%% more data processing - find above stats per delta SNR
% easy false alarm probs
zerorows = find(and(any(data(:,[2])==0, 2), any(data(:,[3]) == 0, 2))); zerodata = data(zerorows,:);
zerovar = find(and(any(data(:,[2])==0, 2), any(data(:,[3]) ~= 0, 2))); zerovardata = data(zerovar,:);
zeromean = find(and(any(data(:,[2])~=0, 2), any(data(:,[3]) == 0, 2))); zeromeandata = data(zeromean,:);
for i=1:length(zeromeandata)
    zeromeandata(i,10) = zeromeandata(i,5) - (zeromeandata(i,4));
end

delta_snrs = sort(unique(zeromeandata(:,10)));

for snrs_idx = 1:length(delta_snrs)    
    snrsdata = find(zeromeandata(:,10)==delta_snrs(snrs_idx));
    snrsdata = zeromeandata(snrsdata,:);
    
    p_total = 0;
    p_var_total = 0;
    pfa_var_sr = 0; pfa_var_cusum = 0; pfa_var_shewart = 0; pfa_var_modcu = 0;
    pd_var_sr = 0; pd_var_cusum = 0; pd_var_shewart = 0; pd_var_modcu = 0;
    
    add_sr = [length(n)];    add_cusum = [length(n)];    add_shewart = [length(n)];    add_modcu = [length(n)];
    

    % zero mean data - delta var
    for idx = 1:size(snrsdata,1)
        if (snrsdata(idx,6)<length(n))
            pd_var_cusum = pd_var_cusum + 1;
            if ( snrsdata(idx,6)>=200)
                add_cusum(end + 1) = snrsdata(idx,6)-200;
            end
        end
        if (snrsdata(idx,6)<200)
            pfa_var_cusum = pfa_var_cusum + 1;
        end
        if ( snrsdata(idx,7)<length(n))
            pd_var_shewart = pd_var_shewart + 1;
            
            if ( snrsdata(idx,7)>=200)
                add_shewart(end + 1) = snrsdata(idx,7)-200;
            end
        end
        if (snrsdata(idx,7)<200)
            pfa_var_shewart = pfa_var_shewart + 1;
        end
        if ( snrsdata(idx,8)<length(n))
            pd_var_sr = pd_var_sr + 1;
            if ( snrsdata(idx,8)>=200)
                add_sr(end + 1) = snrsdata(idx,8)-200;
            end
        end
        if (snrsdata(idx,8)<200)
            pfa_var_sr = pfa_var_sr + 1;
        end
        if ( snrsdata(idx,9)<length(n))
            pd_var_modcu = pd_var_modcu + 1;
            
            if ( snrsdata(idx,9)>=200)
                add_modcu(end + 1) = snrsdata(idx,9)-200;
            end
        end
        if (snrsdata(idx,9)<200)
            pfa_var_modcu = pfa_var_modcu + 1;
        end
        
        p_var_total = p_var_total + 1;
    end
    
    add_cusum_var = sum(add_cusum)/numel(add_cusum);
    add_shewart_var = sum(add_shewart)/numel(add_shewart);
    add_sr_var = sum(add_sr)/numel(add_sr);
    add_modcu_var = sum(add_modcu)/numel(add_modcu);

    % save the stuff into global arrays
    var_data(snrs_idx,:) = [...
        add_cusum_var, ...
        add_shewart_var,...
        add_sr_var,...
        add_modcu_var,...
        pfa_var_cusum/p_var_total,...
        pfa_var_shewart/p_var_total,...
        pfa_var_sr/p_var_total,...
        pfa_var_modcu/p_var_total,...
        pd_var_cusum/p_var_total,...
        pd_var_shewart/p_var_total,...
        pd_var_sr/p_var_total,...
        pd_var_modcu/p_var_total,...
        delta_snrs(snrs_idx),...
        ];
end

%% plot some stuff
figure(89);clf; 
subplot(4,1,1)
grid minor; hold on
plot(-var_data(:,13), var_data(:,12),'--')
plot(-var_data(:,13), var_data(:,8),'-.')
title('Pd and Pfa for different changes in SNR')
xlabel('Change in SNR (dB) - Modified CUSUM')
legend('P_D','P_F_A')

subplot(4,1,2)
grid minor; hold on
plot(-var_data(:,13), var_data(:,9),'--')
plot(-var_data(:,13), var_data(:,5),'-.')
% title('P_D and P_F_A for different changes in SNR')
xlabel('Change in SNR (dB) -  CUSUM')
legend('P_D','P_F_A')

subplot(4,1,3)
grid minor; hold on
plot(-var_data(:,13), var_data(:,10),'--')
plot(-var_data(:,13), var_data(:,6),'-.')
% title('P_D and P_F_A for different changes in SNR')
xlabel('Change in SNR (dB) - Shewart')
legend('P_D','P_F_A')

subplot(4,1,4)
grid minor; hold on
plot(-var_data(:,13), var_data(:,11),'--')
plot(-var_data(:,13), var_data(:,7),'-.')
% title('P_D and P_F_A for different changes in SNR')
xlabel('Change in SNR (dB) - SR Gauss')
legend('P_D','P_F_A')


f = @(b,x) b(1).*exp(b(2).*x) + b(3);
nrmrsd = @(b) norm( var_data(:,1) - f(b,-var_data(:,13)))                          % Residual Norm Cost Function
B0 = rand(3,1);                                                     % Choose Appropriate Initial Estimates
[B,rnrm] = fminsearch(nrmrsd, B0);    
f1 = f(B,-var_data(:,13))
nrmrsd = @(b) norm( var_data(:,2) - f(b,-var_data(:,13)))                          % Residual Norm Cost Function
B0 = rand(3,1);                                                     % Choose Appropriate Initial Estimates
[B,rnrm] = fminsearch(nrmrsd, B0);    
f2 = f(B,-var_data(:,13))
nrmrsd = @(b) norm( var_data(:,3) - f(b,-var_data(:,13)))                          % Residual Norm Cost Function
B0 = rand(3,1);                                                     % Choose Appropriate Initial Estimates
[B,rnrm] = fminsearch(nrmrsd, B0);    
f3 = f(B,-var_data(:,13))
nrmrsd = @(b) norm( var_data(:,4) - f(b,-var_data(:,13)))                          % Residual Norm Cost Function
B0 = rand(3,1);                                                     % Choose Appropriate Initial Estimates
[B,rnrm] = fminsearch(nrmrsd, B0);    
f4 = f(B,-var_data(:,13))

figure (90);clf;
grid minor; hold on;
h1=plot(-var_data(:,13),f1,'b-.',-var_data(:,13),var_data(:,1),'b.')
h2=plot(-var_data(:,13),f2,'r-.',-var_data(:,13),var_data(:,2),'r.')
h3=plot(-var_data(:,13),f3,'m-.',-var_data(:,13),var_data(:,3),'m.')
h4=plot(-var_data(:,13),f4,'g-.',-var_data(:,13),var_data(:,4),'g.')
legend([h1(2) h2(2) h3(2) h4(2)],'CUSUM','Shewart','SR Gauss','Modified CUSUM')
xlabel('Change in SNR (dB)')
ylabel('Average number of steps before detection')
%% more data processing - find above stats per delta SNR
% easy false alarm probs
zerorows = find(and(any(data(:,[2])==0, 2), any(data(:,[3]) == 0, 2))); zerodata = data(zerorows,:);
zerovar = find(and(any(data(:,[2])==0, 2), any(data(:,[3]) ~= 0, 2))); zerovardata = data(zerovar,:);
zeromean = find(and(any(data(:,[2])~=0, 2), any(data(:,[3]) == 0, 2))); zeromeandata = data(zeromean,:);
for i=1:length(zeromeandata)
    zerovardata(i,10) = zerovardata(i,3);
end

delta_snrs = sort(unique(zerovardata(:,10)));
var_data = [];
for snrs_idx = 1:length(delta_snrs)    
    snrsdata = find(zerovardata(:,10)==delta_snrs(snrs_idx));
    snrsdata = zerovardata(snrsdata,:);
    
    p_total = 0;
    p_var_total = 0;
    pfa_var_sr = 0; pfa_var_cusum = 0; pfa_var_shewart = 0; pfa_var_modcu = 0;
    pd_var_sr = 0; pd_var_cusum = 0; pd_var_shewart = 0; pd_var_modcu = 0;
    
    add_sr = [length(n)];    add_cusum = [length(n)];    add_shewart = [length(n)];    add_modcu = [length(n)];
    

    % zero mean data - delta var
    for idx = 1:size(snrsdata,1)
        if (snrsdata(idx,6)<length(n))
            pd_var_cusum = pd_var_cusum + 1;
            if ( snrsdata(idx,6)>=500)
                add_cusum(end + 1) = snrsdata(idx,6)-500;
            end
        end
        if (snrsdata(idx,6)<500)
            pfa_var_cusum = pfa_var_cusum + 1;
        end
        if ( snrsdata(idx,7)<length(n))
            pd_var_shewart = pd_var_shewart + 1;
            
            if ( snrsdata(idx,7)>=500)
                add_shewart(end + 1) = snrsdata(idx,7)-500;
            end
        end
        if (snrsdata(idx,7)<500)
            pfa_var_shewart = pfa_var_shewart + 1;
        end
        if ( snrsdata(idx,8)<length(n))
            pd_var_sr = pd_var_sr + 1;
            if ( snrsdata(idx,8)>=500)
                add_sr(end + 1) = snrsdata(idx,8)-500;
            end
        end
        if (snrsdata(idx,8)<500)
            pfa_var_sr = pfa_var_sr + 1;
        end
        if ( snrsdata(idx,9)<length(n))
            pd_var_modcu = pd_var_modcu + 1;
            
            if ( snrsdata(idx,9)>=500)
                add_modcu(end + 1) = snrsdata(idx,9)-500;
            end
        end
        if (snrsdata(idx,9)<500)
            pfa_var_modcu = pfa_var_modcu + 1;
        end
        
        p_var_total = p_var_total + 1;
    end
    
    add_cusum_var = sum(add_cusum)/numel(add_cusum);
    add_shewart_var = sum(add_shewart)/numel(add_shewart);
    add_sr_var = sum(add_sr)/numel(add_sr);
    add_modcu_var = sum(add_modcu)/numel(add_modcu);

    % save the stuff into global arrays
    var_data(snrs_idx,:) = [...
        add_cusum_var, ...
        add_shewart_var,...
        add_sr_var,...
        add_modcu_var,...
        pfa_var_cusum/p_var_total,...
        pfa_var_shewart/p_var_total,...
        pfa_var_sr/p_var_total,...
        pfa_var_modcu/p_var_total,...
        pd_var_cusum/p_var_total,...
        pd_var_shewart/p_var_total,...
        pd_var_sr/p_var_total,...
        pd_var_modcu/p_var_total,...
        delta_snrs(snrs_idx),...
        ];
end

%% plot some stuff
figure(891);clf; 
subplot(4,1,1)
grid minor; hold on
plot(var_data(:,13), var_data(:,12),'--')
plot(var_data(:,13), var_data(:,8),'-.')
title('Pd and Pfa for different changes in mean')
xlabel('Change in mean - Modified CUSUM')
legend('P_D','P_F_A')

subplot(4,1,2)
grid minor; hold on
plot(var_data(:,13), var_data(:,9),'--')
plot(var_data(:,13), var_data(:,5),'-.')
% title('Pd and Pfa for different changes in mean')
xlabel('Change in mean -  CUSUM')
legend('P_D','P_F_A')

subplot(4,1,3)
grid minor; hold on
plot(var_data(:,13), var_data(:,10),'--')
plot(var_data(:,13), var_data(:,6),'-.')
% title('Pd and Pfa for different changes in mean')
xlabel('Change in mean - Shewart')
legend('P_D','P_F_A')

subplot(4,1,4)
grid minor; hold on
plot(var_data(:,13), var_data(:,11),'--')
plot(var_data(:,13), var_data(:,7),'-.')
% title('Pd and Pfa for different changes in mean')
xlabel('Change in mean - SR Gauss')
legend('P_D','P_F_A')


f = @(b,x) b(1).*exp(b(2).*x) + b(3);
nrmrsd = @(b) norm( var_data(:,1) - f(b,-var_data(:,13)))                          % Residual Norm Cost Function
B0 = rand(3,1);                                                     % Choose Appropriate Initial Estimates
[B,rnrm] = fminsearch(nrmrsd, B0);    
f1 = f(B,-var_data(:,13))
nrmrsd = @(b) norm( var_data(:,2) - f(b,-var_data(:,13)))                          % Residual Norm Cost Function
B0 = rand(3,1);                                                     % Choose Appropriate Initial Estimates
[B,rnrm] = fminsearch(nrmrsd, B0);    
f2 = f(B,-var_data(:,13))
nrmrsd = @(b) norm( var_data(:,3) - f(b,-var_data(:,13)))                          % Residual Norm Cost Function
B0 = rand(3,1);                                                     % Choose Appropriate Initial Estimates
[B,rnrm] = fminsearch(nrmrsd, B0);    
f3 = f(B,-var_data(:,13))
nrmrsd = @(b) norm( var_data(:,4) - f(b,-var_data(:,13)))                          % Residual Norm Cost Function
B0 = rand(3,1);                                                     % Choose Appropriate Initial Estimates
[B,rnrm] = fminsearch(nrmrsd, B0);    
f4 = f(B,-var_data(:,13))

figure (901);clf;
grid minor; hold on;
h5=plot(var_data(:,13),var_data(:,1),'b.')
h6=plot(var_data(:,13),var_data(:,2),'r.')
h7=plot(var_data(:,13),var_data(:,3),'m.')
h8=plot(var_data(:,13),var_data(:,4),'g.')
legend([h5(1) h6(1) h7(1) h8(1)],'CUSUM','Shewart','SR Gauss','Modified CUSUM')
xlabel('Change in mean')
ylabel('Average number of steps before detection')
ylim([0 50])