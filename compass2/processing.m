%% Processing for compass2
clear all
clc

%% Read in images
cd('images\cropped')
mats = dir('Image*.png');
for i=1:length(mats)
    im{i} = imread(mats(i).name);
%     filename = sprintf('Image%02d.png',i)
%     imwrite(im(i).Image ,filename)
end
cd('..\..')

cd('images\tumors')
mats = dir('tumor*.png');
for i=1:length(mats)
    tu{i} = imread(mats(i).name);
%     filename = sprintf('Image%02d.png',i)
%     imwrite(im(i).Image ,filename)
end
cd('..\..')
%% Build the statistics for noise and signal based on independent sub strips
windowSize = 40;
stripSize = fix(windowSize/2);
saveNoise = 0;
saveSignal = 0;

%% Build noise matrix first
imgCnt = 0;
noiseImages = [];
for k=1:length(im)
    display(['Processing noise image ' num2str(k)])
    currImage = im{k};
    % extract columns from ROI of windowSize x windowSize
    for i=1:windowSize:size(currImage,1)-windowSize
        for j=1:windowSize:size(currImage,2)-windowSize
            imgCnt= imgCnt+1;
            ROI = (currImage(i:i+windowSize-1,j: j+ windowSize -1));
            % loop over ROI and get column and store
            for x=1:windowSize
                for y=1:stripSize:windowSize-stripSize
                    noiseImages(end+1,:) = ROI(y:y+stripSize-1,x);
                end
            end
            if (saveNoise)
                filename = sprintf('images/noise/noise%06d.png',imgCnt);
                imwrite(ROI ,filename)
            end
        end
    end
end
%% Calculate stats for noise
muN = mean(noiseImages);
sumNorm2 = 0;
for i=1:size(noiseImages,1)
    sumNorm2 = sumNorm2 + norm(noiseImages(i,:))^2;
end
varN = sumNorm2/(stripSize*size(noiseImages,1))

Cn = cov(noiseImages);
[Un En] = eig(inv(Cn));
En(En<0)=0;
Dn=Un*sqrt(En)*Un';
whiteCn = cov(noiseImages*Dn);

%% Now build signal matrix
muS = 0;
imgCnt = 0;
signalImages = [];
tumorSize = 45;

for k=1:length(tu)-1
    display(['Processing signal image' num2str(k)])
    currImage = tu{k};
%     figure(4);clf;imshow(currImage,[])
    currImage = currImage( floor((45-tumorSize)/2) + 1 : floor((45-tumorSize)/2) + tumorSize , floor((45-tumorSize)/2) + 1 :  floor((45-tumorSize)/2) + tumorSize);
    % extract columns from ROI of windowSize x windowSize
    for i=1:windowSize:size(currImage,1)-windowSize
        for j=1:windowSize:size(currImage,2)-windowSize
            imgCnt= imgCnt+1;
            ROI = (currImage(i:i+windowSize-1,j: j+ windowSize -1));
            % loop over ROI and get column and store
            for x=1:windowSize
                for y=1:stripSize:windowSize-stripSize
                    signalImages(end+1,:) = ROI(y:y+stripSize-1,x);
                end
            end
            if (saveNoise)
                filename = sprintf('images/signal/signal%06d.png',imgCnt);
                imwrite(ROI ,filename)
            end
        end
    end
end
%% Get signal stats
muS = mean(signalImages);

Cs = cov(signalImages);
[Us Es] = eig(inv(Cs));
Es(Es<0)=0;
Ds=Us*sqrt(Es)*Us';
whiteCs = cov(signalImages*Ds);

[UU DD VV] = svd(signalImages','econ');

% find p<stripSize such that sum_{i=1}^p DD(i) / sum_{i=1}^p DD(i) = .95
lamSquares = diag(DD).^2;
maxEnergy = sum(lamSquares);
partialEnergy = 0;
numPrincipleVectors = 0;
threshEnergy = .9999;
while (numPrincipleVectors<stripSize && partialEnergy < threshEnergy*maxEnergy)
    numPrincipleVectors = numPrincipleVectors + 1;
    partialEnergy = partialEnergy + lamSquares(numPrincipleVectors);
end
H = UU(:,1:numPrincipleVectors);

%% Save
save('images\stats.mat','Cn','windowSize','stripSize','varN','muN','noiseImages','muS','Cs','signalImages','numPrincipleVectors')