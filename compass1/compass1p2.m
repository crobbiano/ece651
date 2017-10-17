%% Save the boobies!
%  ECE651 computer assignment 1 part 2
clear all
clc
%% Read in test images
cd('images\cropped');
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
%% make signals statistics
sig=zeros(11,2025);
signalMean = [];
signalCov = [];
for i=1:length(tumors)
    sig(i,:) = reshape(double(tumors(i).currTumor),[],1);
    subplot(4,3,i);
    imshow(reshape(sig(i,:),45,45),[])
end
%% Whighten data by y=Ds, D=C^{1/2}
[U E] = eig(inv(C));
E(E<0)=0;
D=U*sqrt(E)*U';
sigWhite=zeros(2025,12);
for i=1:12
    sigWhite(:,i) = D*sig(i,:)';
    subplot(4,3,i);
%     imshow(reshape(sigWhite(:,i),45,45),[])
end


%% Look at just the first image for now and scan windows for tumors
%  for each window, check the test statistics and determine what hypothesis
%  should be accepted
% currImage = images(1).Image;
currImage = images{6};
% Cheat and test against a signal
% currImage = reshape(sig(6,:),45,45);

detects = zeros(1,12);
plotThings = 0;
if (plotThings)
    figure(234);clf;
end
coords=[];
for i=1:1:size(currImage,1)-windowSize+1
    for j=1:1:size(currImage,2)-windowSize+1
        window = D*(reshape(double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),2025,1));
%         window = D*(reshape(double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),2025,1)-mu');
        
%         window = reshape(double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),[],1)-mu';
%         hist(window,255);
        
        tests = [];
        for l=1:4
%             tests(l)=window'*sig(l,:)' - .5*sig(l,:)*sig(l,:)';
            tests(l)=window'*sigWhite(:,l) - .5*sigWhite(:,l)'*sigWhite(:,l);
        end
        [maxim,idx] = max(tests);
        if (idx ~= 1)
            coords(end+1,:) = [ (j+windowSize/2) (i+windowSize/2) idx];
        end
%         display(['max: ' num2str(idx)])
        detects(idx) = detects(idx) + 1;
        if (plotThings)
            subplot(1,3,1)
            %         imshow(D*double(currImage(i:windowSize+(i-1),j:windowSize+(j-1))),[]);
            imshow(double(currImage(i:windowSize+(i-1),j:windowSize+(j-1)))-reshape(mu,45,45),[]);
            xlabel(['(' num2str(j+windowSize/2) ',' num2str(i+windowSize/2) ')'])
            subplot(1,3,2)
            imshow(reshape(sigWhite(:,idx),45,45),[]);
            xlabel(['signal: ' num2str(idx)])
            subplot(1,3,3)
            imshow(reshape(window,45,45),[]);
            xlabel('window')
            pause(.0001)
        end
    end
    1;
end
detects