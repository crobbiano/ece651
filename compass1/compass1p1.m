%% Save the boobies!
%  ECE651 computer assignment 1 part 1
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
load('images\noise\stats.mat')
%% Read in tumors
cd('images\tumors');
ts = dir('tumor*.mat');
for i=1:length(ts)
    tumors(i) = load(ts(i).name);
end
cd('..\..')
%% make signals 
sig=zeros(12,windowSize*windowSize);
sig(1,:) = noiseImages(1,:)-mu;
figure(2);clf;
subplot(4,3,1);
imshow(reshape(sig(1,:),windowSize,windowSize),[])
for i=1:length(tumors)
%     sig(i+1,:) = reshape(double(tumors(i).currTumor),[],1);
    % Find the center of the tumor and extract outwards from there
    currTumor = double(tumors(i).currTumor);
    currTumor = currTumor( floor((45-windowSize)/2) + 1 : end- floor((45-windowSize)/2), floor((45-windowSize)/2) + 1 : end - floor((45-windowSize)/2));
    sig(i+1,:) = reshape(currTumor,[],1)-mu';
%     sig(sig<0)=0;
    subplot(4,3,i+1);
    imshow(reshape(sig(i+1,:),windowSize,windowSize),[])
end
%% Whighten data by y=Ds, D=C^{1/2}
[U E] = eig(inv(C));
E(E<0)=0;
D=U*sqrt(E)*U';
sigWhite=zeros(windowSize*windowSize,12);
sigThresh = zeros(1,12);
for i=1:12
    sigWhite(:,i) = D*sig(i,:)';
    currSig = sigWhite(:,i);
    sigProd = currSig'*currSig;
%     sigThresh(i) = 3.7*sqrt(sigProd);
%     sigThresh(i) = 3.7*sqrt(sigProd);
    subplot(4,3,i);
%     imshow(reshape(sigWhite(:,i),windowSize,windowSize),[])
end


%% Look at just the first image for now and scan windows for tumors
%  for each window, check the test statistics and determine what hypothesis
%  should be accepted
% currImage = images(1).Image;
currImage = images{9};
% Cheat and test against a signal
% currImage = reshape(sig(6,:),windowSize,windowSize);
% currImage = imread('test2.png');

detects = zeros(1,12);
detectMap = [];
plotThings = 0;
if (plotThings)
    figure(234);clf;
end
coords=[];
iIdx=0;jIdx=0;
for i=1:1:size(currImage,1)-windowSize+1
    iIdx= iIdx +1;
    jIdx=0;
    for j=1:1:size(currImage,2)-windowSize+1
        jIdx= jIdx +1;
%         window = D*(reshape(double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),windowSize*windowSize,1));
        window = D*(reshape(double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),windowSize*windowSize,1)-mu');
        
%         window = reshape(double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),[],1)-mu';
%         hist(window,255);
        

        tests = -inf(1,12);
%         for l=[2 3 4 5 6 7 8 9 10 11 12]
        for l=[ 6 7 8 9]
%             tests(l)=window'*sig(l,:)' - .5*sig(l,:)*sig(l,:)';
            tests(l)=window'*sigWhite(:,l) - .5*sigWhite(:,l)'*sigWhite(:,l);
%             tests(l)=window'*sigWhite(:,l);
        end
        
        if (sum(tests(2:end)>sigThresh(2:end))>0)
            detected = 1;
        else
            detected = 0;
        end
        
        [maxim,idx] = max(tests);
        detectMap(iIdx,jIdx) = tests(idx);
        
        if (idx ~= 1)
            coords(end+1,:) = [ (j+windowSize/2) (i+windowSize/2) idx];
        end
%         display(['max: ' num2str(idx)])
        detects(idx) = detects(idx) + 1;
        if (plotThings)
            subplot(1,3,1)
            %         imshow(D*double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),[]);
            imshow(double(currImage(i:windowSize+(i-1),j:windowSize+(j-1)))-reshape(mu,windowSize,windowSize),[]);
            xlabel(['(' num2str(j+windowSize/2) ',' num2str(i+windowSize/2) ')'])
            subplot(1,3,2)
            imshow(reshape(sigWhite(:,idx),windowSize,windowSize),[]);
            xlabel(['signal: ' num2str(idx)])
            subplot(1,3,3)
            imshow(reshape(window,windowSize,windowSize),[]);
            xlabel('window')
            pause(.0001)
        end
    end
    1;
end
detects

%% processes the detectMap
% look at single pixel, then if majority of neighbors passes threshold then
% decide there is a detection
[m,n]= size(detectMap);
detections = zeros(m,n);
% thresh = -2200;
thresh = -.65*max(max(detectMap));
thresh = 1.45*max(max(detectMap));
for i=2:m-1
    for j=2:n-1
        
        avgT = detectMap(i-1,j-1) + detectMap(i,j-1) + detectMap(i+1,j-1) + detectMap(i-1,j) + detectMap(i,j) + ...
            detectMap(i+1,j) + detectMap(i-1,j+1) +detectMap(i,j+1)+ detectMap(i+1,j+1);
        avgT = avgT / 9;
%         detections(i,j) = avgT;
        
        if (avgT>thresh)
           detections(i,j) = 1;
        end
    end
end
figure(33);imagesc(detections)