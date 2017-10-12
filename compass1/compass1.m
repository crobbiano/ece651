%% Save the boobies!
%  ECE651 computer assignment 1
clear all
clc
%% Read in test images
cd('images\original');
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
sig=zeros(2025,12);
sig(:,1) = reshape(meanImage,2025,[]);
for i=1:length(tumors)
    sig(:,i+1) = reshape(double(tumors(i).currTumor),[],1);
%     sig(:,i+1) = reshape(double(tumors(i).currTumor)-meanImage,[],1);
end
%% Whighten data by y=Ds, D=C^{1/2}
[U E] = eig(sampleCov);
D=U*sqrt(E)*U';
sigWhite=zeros(2025,12);
for i=1:11
    sigWhite(:,i+1) = reshape( D*reshape(sig(:,i+1),45,45),[],1);
end


%% Look at just the first image for now and scan windows for tumors
%  for each window, check the test statistics and determine what hypothesis
%  should be accepted
currImage = images{5};

% outs = conv2(double(currImage), reshape(sig(:,2),45,45)/sum(sig(:,2)));
% figure(83);clf;imshow(reshape(sig(:,2),45,45)/sum(sig(:,2)),[])
% figure(89);clf;
% subplot(2,1,1)
% imshow(outs,[])
% subplot(2,1,2)
% matchingRegions = abs(outs - 1) < 0.01;
% % Use axes() or figure() to switch to a new axes if you want.
% imshow(matchingRegions, []);

detects = zeros(1,12);
figure(234);clf;
for i=1:size(currImage,1)-windowSize
    for j=1:size(currImage,1)-windowSize
%         window = reshape(D*double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),[],1);
        window = reshape(double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),[],1);
%         hist(window,255);
        
        tests = [];
        for l=1:12
%             window'*inv(sampleCov)*reshape(sig(:,l),45,45)
            tests(l)=reshape(window,1,[])*sig(:,l);
%             tests(l)=reshape(window,1,[])*sigWhite(:,l);
%             tests(l)=reshape(window,1,[])*sigWhite(:,l) - .5*sigWhite(:,l)'*sigWhite(:,l);
        end
        [maxim,idx] = max(tests);
%         display(['max: ' num2str(idx)])
        detects(idx) = detects(idx) + 1;
        subplot(3,1,1)
%         imshow(D*double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),[]);
        imshow(double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),[]);
        ylabel(['(' num2str(i) ',' num2str(j) ')'])
        subplot(3,1,2)
        imshow(reshape(sig(:,idx),45,45),[]);
%         imshow(reshape(sigWhite(:,idx),45,45),[]);
        ylabel(['signal: ' num2str(idx)])
    end
    1;
end
detects