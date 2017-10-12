%% Save the boobies!
%  ECE651 computer assignment 1
clear all
clc
%% Read in test images
cd('images\cropped');
ts = dir('cropped*.png');
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
sig(:,1) = meanImage;
figure(2);clf;
subplot(4,3,1);
imshow(reshape(sig(:,1),45,45),[])
for i=1:length(tumors)
    sig(:,i+1) = reshape(double(tumors(i).currTumor),[],1);
%     sig(:,i+1) = reshape(double(tumors(i).currTumor),[],1)-meanImage;
    sig(sig<0)=0;
    subplot(4,3,i+1);
    imshow(reshape(sig(:,i+1),45,45),[])
end
%% Whighten data by y=Ds, D=C^{1/2}
[U E] = eig(sampleCov);
E(E<.00001)=0;
D=U*sqrt(E)*U';
sigWhite=zeros(2025,12);
for i=1:12
    sigWhite(:,i) = D*sig(:,i);
end


%% Look at just the first image for now and scan windows for tumors
%  for each window, check the test statistics and determine what hypothesis
%  should be accepted
currImage = images{10};

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
for i=400:size(currImage,1)-windowSize
    for j=1:size(currImage,2)-windowSize
        window = D*(reshape(double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),2025,1));
%         window = D*(reshape(double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),2025,1)-meanImage);
%         window = reshape(double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),[],1)-meanImage;
%         hist(window,255);
        
        tests = [];
        for l=1:3
%             window'*inv(sampleCov)*reshape(sig(:,l),45,45)
%             tests(l)=window'*sig(:,l) - .5*sig(:,l)'*sig(:,l);
%             tests(l)=window'*sig(:,l);
            tests(l)=window'*sigWhite(:,l)- .5*sigWhite(:,l)'*sigWhite(:,l);
%             tests(l)=reshape(window,1,[])*sigWhite(:,l);
%             tests(l)=reshape(window,1,[])*sigWhite(:,l) - .5*sigWhite(:,l)'*sigWhite(:,l);
        end
        [maxim,idx] = max(tests);
%         display(['max: ' num2str(idx)])
        detects(idx) = detects(idx) + 1;
        subplot(1,3,1)
%         imshow(D*double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),[]);
        imshow(double(currImage(i:windowSize+(i-1),j:windowSize+(j-1)))-reshape(meanImage,45,45),[]);
        ylabel(['(' num2str(i) ',' num2str(j) ')'])
        subplot(1,3,2)
        imshow(reshape(sig(:,idx),45,45),[]);
        ylabel(['signal: ' num2str(idx)])
        subplot(1,3,3)
        imshow(reshape(sig(:,1),45,45),[]);
        pause(.01)
    end
    1;
end
detects