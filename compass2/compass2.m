%% Save the boobies!
%  ECE651 computer assignment 2
clear all
clc
%% Read in test images
cd('images\croppedOrig');
% ts = dir('test*.png');
ts = dir('Image*.png');
for i=1:length(ts)
    %     images(i) = load(ts(i).name);
    images{i} = imread(ts(i).name);
end
cd('..\..')
%% Read in sample cov and mean
load('images\stats.mat')

%% Look at just the first image for now and scan windows for tumors
%  for each window, check the test statistics and determine what hypothesis
%  should be accepted

% Cheat and test against a signal
% currImage = reshape(sig(6,:),45,45);

plotThings = 0;

coords=[];
% C1inv = inv(Cs + Cn);
% C0inv = inv(Cn);
% Cest = Cs*C1inv;

slideLen = 1;

% for picIdx=1:length(images)
detectMap = [];
Pds = [];
Pfas = [];
picNum = 0;
for picIdx=[1 3 7 8 9 ]
% for picIdx=[1 2]
% for picIdx=[2]
    picNum = picNum +1 ;
    tic
    currImage = images{picIdx};
    %     currImage = imread('images_p2\tumors\tumor3.png');
    %     currImage = imread('test2.png');
    
    detectMap = [];
    
    testVals = zeros(size(currImage,1)-windowSize+1,size(currImage,2)-windowSize+1);
    display(['Image: ' num2str(picIdx)])
    iIdx=0;jIdx=0;
    for i=1:slideLen:size(currImage,1)-windowSize+1
        iIdx = iIdx+1;jIdx=0;
        if (mod(i,20)==0)
            display(['i: ' num2str(i)])
        end
        for j=1:slideLen:size(currImage,2)-windowSize+1
            jIdx=jIdx+1;
            
            %             window = (reshape(double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),2025,1)-mu');
            ROI = currImage(i:windowSize+(i-1),j:windowSize+(j-1));
            % Take current window and split it into sub windows and process
            i2Idx=0;j2Idx=0;
            t=[];
            for i2=1:windowSize
                i2Idx = i2Idx + 1;
                j2Idx=0;
                for j2=1:stripSize:windowSize-stripSize
                    j2Idx = j2Idx + 1;
                    %                     subWindow(:,i2Idx) = ROI(j2:j2+stripSize-1,i2);
                    subWindow = ROI(j2:j2+stripSize-1,i2);
                    t(i2Idx,j2Idx) = norm(double(subWindow)'*H)^2;
                end
            end
            [mm,nn]=size(t);
            %             bigMap((iIdx-1)*mm+1:(iIdx-1)*mm+mm,(jIdx-1)*nn+1:(jIdx-1)*nn+nn) = t;
            detects = sum(sum(t))/numel(t);
            detectMap(iIdx,jIdx) = detects;
        end
    end
    toc
    %     filename = sprintf('images_p2/results/result%02d.mat',picIdx)
    %     imwrite(testVals ,filename)
    %     save(filename, 'testVals');
    
    %% processes the detectMap
    % look at single pixel, then if majority of neighbors passes threshold then
    % decide there is a detection
    
    % loop over different values of the thresh and do the following:
    %   find center of tumor
    %   determine radius from windowSize
    %   count # of 1's in detections - if this # == # of pixels in radius then we have 100% Pd
    %   count # of 1's everywhere else - this is Pfa
    centers = [100,128;...
        190,160;...
        207,147;...
        55,65;...
        125,135;...
        105,145;...
        195,260;...
        55,50;...
        195,145;...
        50,176;];
    radiuses = [30,33,40,29,44,20,20,29,26,45];
    
    for cent=1:size(centers,1)
        %         centers(cent,:) = centers(cent,:) - [ceil(windowSize/4), ceil(windowSize/4)];
        centers(cent,:) = centers(cent,:) - [floor((windowSize-10)/4), floor((windowSize-10)/4)];
    end
    
    [m,n]= size(detectMap);
%     foldername = sprintf('images/results/image%02d',picIdx);
%     mkdir(foldername)
%     filename = sprintf('images/results/image%02d/result%02d_orig.png',picIdx,picIdx);
%     imwrite(currImage ,filename)        
%     filename = sprintf('images/results/image%02d/result%02d_detMap.png',picIdx,picIdx)
%     imwrite(mat2gray(detectMap),filename)
    
    
    pIdx = 0;
        for threshInc=0:.01:1
%     for threshInc=.87
        pIdx = pIdx + 1;
        detections = zeros(m,n);
        thresh =threshInc*max(max(detectMap));
        %         thresh = 42480115.2;
        for i=1:m
            for j=1:n
                if (detectMap(i,j)>thresh)
                    detections(i,j) = 1;
                end
            end
        end
        
        % Get center and radius
        radius = ceil(45/2);
        center = centers(picIdx,:) - [radius, radius];
        radius = floor(radiuses(picIdx)/3);
        % fix center for tumors by edges
        if (center(1)<=radius)
            center(1) = radius+1;
        end
        if (center(2)<=radius)
            center(2) = radius+1;
        end
        %         radius = ceil(windowSize/4);
        %         radius = 3;
        % Count # of 1's ball around radius
        sigWindow = detections(center(2)-radius:center(2)+radius-1, center(1)-radius:center(1)+radius-1);
        numInSigWindow = sum(sum(sigWindow));
        Pds(picNum,pIdx) = numInSigWindow / (2*radius*2*radius);
        % break image into 4 out of signal sections
        numNotInSigWindow = sum(sum( detections(1:center(2)-radius-1, :) ));
        numNotInSigWindow = numNotInSigWindow + sum(sum( detections(center(2)+radius:end, :) ));
        numNotInSigWindow = numNotInSigWindow + sum(sum( detections(center(2)-radius:center(2)+radius-1, 1:center(1)-radius-1) ));
        numNotInSigWindow = numNotInSigWindow + sum(sum( detections(center(2)-radius:center(2)+radius-1, center(1)+radius:end) ));
        Pfas(picNum,pIdx) = numNotInSigWindow / (m*n - 2*2*radius*radius);
        
        
% 
%         filename = sprintf('images/results/image%02d/result%02d_t%.2f_detections.png',picIdx,picIdx,threshInc);
%         imwrite(detections ,filename)
        detectMapID = detectMap;
        detectMapID(center(2)-radius:center(2)+radius-1, center(1)-radius:center(1)+radius-1) = 0;
        detectionsID = detections;
        detectionsID(center(2)-radius:center(2)+radius-1, center(1)-radius:center(1)+radius-1) = checkerboard(1,length(center(2)-radius:center(2)+radius-1)/2);
        if (plotThings)
            figure(2);clf;
            subplot(1,3,1)
            % imagesc(detectMap)
            imshow(detectMapID,[])
            subplot(1,3,2)
            imshow(currImage,[])
            subplot(1,3,3)
            imshow(detectionsID,[])
            pause(.001)
        end
    end
    
end

%% Reinterpolate the pds and pfas
% Pfas(1,lastIdx:end),Pds(1,lastIdx:end),
for pfaIdx=1:size(Pfas,1)
    % pfaIdx = 3;
    lastIdx = find (Pfas(pfaIdx,:)==1);lastIdx = lastIdx(end);
    [C,ia,idx] = unique(Pfas(pfaIdx,lastIdx:end),'stable');
    val = accumarray(idx,Pds(pfaIdx,lastIdx:end),[],@mean); %You could take something other than the mean.
    newPds(pfaIdx,:)= interp1(C,val,linspace(1,0,1000),'linear'); %see interp1() for other interpolation methods. Extrapolation is dange
end
newPdmean=mean(newPds);
% interp1(Pfas(6,36:end),Pds(6,36:end),linspace(1,0,1000))
%% At this point we should have all of the Pd and Pfa info, now just take
% the mean of each then plot!
Pdmean = mean(Pds);
Pfamean = mean(Pfas);
figure(13);clf; hold on
% plot(Pfamean,Pdmean,'--o')
plot(linspace(1,0,1000),newPdmean,'--')
% plot(Pfas,Pds)
for i=1:size(Pfas,1)
    plot(Pfas(i,:),Pds(i,:))
end
xlabel('P_F_A')
ylabel('P_D')
xlim([-.05 1.05])
% xlim([-.05 .2])
ylim([-.05 1.05])
% ylim([.8 1.05])
legend('mean','1','3','7','8','9')
grid minor
%% calc theoretical ROC
theta=H\double(signalImages');
theta=mean(theta,2);
lam = theta'*theta;
lam = lam/varN;
pfa_ideal=chi2cdf(0:.1:100,3,'upper');
% pd_ideal=ncx2cdf(0:.1:100,1,lam/24,'upper');
pd_ideal=ncx2cdf(0:.1:100,3,lam/1,'upper');

figure(292);clf;plot(pfa_ideal,pd_ideal);
xlabel('P_F_A')
ylabel('P_D')
xlim([-.05 1.05])
ylim([-.05 1.05])
grid minor
%% Plot both theoretical and actual
figure(666);clf;
hold on
plot(linspace(1,0,1000),newPdmean)
plot(pfa_ideal,pd_ideal,'--');
legend('Empirical', 'Theoretical')
xlabel('P_F_A')
ylabel('P_D')
xlim([-.05 1.05])
ylim([-.05 1.05])
grid minor