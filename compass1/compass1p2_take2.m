%% Save the boobies!
%  ECE651 computer assignment 1 part 2
clear all
clc
%% Read in test images
cd('images_p2\croppedOrig');
% ts = dir('test*.png');
ts = dir('Image*.png');
for i=1:length(ts)
    %     images(i) = load(ts(i).name);
    images{i} = imread(ts(i).name);
end
cd('..\..')
%% Read in sample cov and mean
load('images_p2\noise\stats2.mat')

%% Look at just the first image for now and scan windows for tumors
%  for each window, check the test statistics and determine what hypothesis
%  should be accepted

% Cheat and test against a signal
% currImage = reshape(sig(6,:),45,45);

plotThings = 0;
if (plotThings)
    figure(234);clf;
end
coords=[];
C1inv = inv(Cs + Cn);
C0inv = inv(Cn);
Cest = Cs*C1inv;

slideLen =1;

% for picIdx=1:length(images)
detectMap = [];
Pds = [];
Pfas = [];
picNum = 0;
for picIdx=[2 3 4 5 6 8 9]
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
            window = currImage(i:windowSize+(i-1),j:windowSize+(j-1));
            % Take current window and split it into sub windows and process
            i2Idx=0;j2Idx=0;
            t=[];
            for i2=1:subWindowSize:size(window,1)-subWindowSize
                i2Idx = i2Idx + 1;
                j2Idx=0;
                for j2=1:subWindowSize:size(window,2)-subWindowSize
                    j2Idx = j2Idx + 1;
                    tempSub = window(i2:subWindowSize+(i2-1),j2:subWindowSize+(j2-1));
                    subWindow(:,i2Idx) = reshape(double(tempSub), subWindowSize*subWindowSize,[]);
                    subWindow(:,i2Idx) = subWindow(:,i2Idx) - muN';
                    %                     subWindow(:,i2Idx) = subWindow(:,i2Idx);
                    %                     subWindow(:,i2Idx) = subWindow(:,i2Idx) - muS' - muN';
                    t(i2Idx,j2Idx) = subWindow(:,i2Idx)'*C1inv*muS' + .5*subWindow(:,i2Idx)'*C0inv*Cest*subWindow(:,i2Idx);
                    %                     t(i2Idx,j2Idx) = .5*subWindow(:,i2Idx)'*C0inv*(Cest*subWindow(:,i2Idx));
                end
            end
            detects = sum(sum(t))/numel(t); % This is the place where we threshold
            if ( sum(sum(detects))>i2Idx*i2Idx/2 )
                detected = 1;
            else
                detected = 0;
            end
            
            detectMap(iIdx,jIdx) = detects;
            1;
            
            %             t = .5*window'*C0inv*Cest*window;
            %             testVals(iIdx,jIdx)=t;
            
            % FIXME - Need to figure out the RHS of the t equation to compare t against - look at notes
            
            %             if (plotThings)
            %                 subplot(1,3,1)
            %                 %         imshow(D*double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),[]);
            %                 imshow(double(currImage(i:windowSize+(i-1),j:windowSize+(j-1)))-reshape(mu,45,45),[]);
            %                 xlabel(['(' num2str(j+windowSize/2) ',' num2str(i+windowSize/2) ')'])
            %                 subplot(1,3,2)
            %                 imshow(reshape(sigWhite(:,idx),45,45),[]);
            %                 xlabel(['signal: ' num2str(idx)])
            %                 subplot(1,3,3)
            %                 imshow(reshape(window,45,45),[]);
            %                 xlabel('window')
            %                 pause(.0001)
            %             end
            
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
    centers = [85,117;...
        172,132;...
        193,130;...
        28,47;...
        99,102;...
        90,134;...
        190,240;...
        32,27;...
        156,123;...
        26,163;]
    
    [m,n]= size(detectMap);
    foldername = sprintf('images_p2/results/image%02d',picIdx);
    mkdir(foldername)
    filename = sprintf('images_p2/results/image%02d/result%02d_orig.png',picIdx,picIdx);
    imwrite(currImage ,filename)
    
    
    pIdx = 0;
    for threshInc=0:.05:1
        pIdx = pIdx + 1;
        detections = zeros(m,n);
        thresh =threshInc*max(max(detectMap));
        for i=1:m
            for j=1:n
                if (detectMap(i,j)>thresh)
                    detections(i,j) = 1;
                end
            end
        end
        
        % Get center and radius
        radius = windowSize/2;
        center = centers(picIdx,:) - [radius, radius];
        % fix center for tumors by edges
        if (center(1)<radius)
            center(1) = radius;
        end
        if (center(2)<radius)
            center(2) = radius;
        end
        radius = ceil(windowSize/2/2);
        % Count # of 1's ball around radius
        sigWindow = detections(center(2)-radius:center(2)+radius-1, center(1)-radius:center(1)+radius-1);
        figure(298);imshow(sigWindow,[]);
        numInSigWindow = sum(sum(sigWindow));
        Pds(picNum,pIdx) = numInSigWindow / (2*radius*2*radius);
        % break image into 4 out of signal sections
        numNotInSigWindow = sum(sum( detections(1:center(2)-radius-1, :) ));
        numNotInSigWindow = numNotInSigWindow + sum(sum( detections(center(2)+radius:end, :) ));
        numNotInSigWindow = numNotInSigWindow + sum(sum( detections(center(2)-radius:center(2)+radius-1, 1:center(1)-radius-1) ));
        numNotInSigWindow = numNotInSigWindow + sum(sum( detections(center(2)-radius:center(2)+radius-1, center(1)+radius:end) ));
        Pfas(picNum,pIdx) = numNotInSigWindow / (m*n - 2*2*radius*radius);
        
        
        %     filename = sprintf('images_p2/results/result%02d_t%.2f_detMap.png',picIdx,threshInc)
        %     imwrite(detectMap ,filename)
        filename = sprintf('images_p2/results/image%02d/result%02d_t%.2f_detections.png',picIdx,picIdx,threshInc);
        imwrite(detections ,filename)
        
        figure(2);clf;
        subplot(1,3,1)
        % imagesc(detectMap)
        imshow(detectMap,[])
        subplot(1,3,2)
        imshow(currImage,[])
        subplot(1,3,3)
        imshow(detections,[])
        
    end
    
end

%% At this point we should have all of the Pd and Pfa info, now just take
% the mean of each then plot!
Pdmean = mean(Pds);
Pfamean = mean(Pfas);
figure(13);clf;
plot(Pfamean,Pdmean)
xlabel('P_F_A')
ylabel('P_D')
xlim([-.05 1.05])
ylim([-.05 1.05])
grid minor