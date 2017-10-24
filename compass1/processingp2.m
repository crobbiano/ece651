%% Processing for part 2
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
%% Build the statistics for noise and signal based on independent sub windows
windowSize = 40;
subWindowSize = 4;
saveNoise = 0;
saveSignal = 0;

%% Build noise stats first

imgCnt = 0;
noiseImages = [];
for k=1:10
    display(['Processing image ' num2str(k)])
    currImage = im{k};
    for i=1:subWindowSize:size(currImage,1)-subWindowSize
        for j=1:subWindowSize:size(currImage,2)-subWindowSize
            imgCnt= imgCnt+1;
            currNoise = (currImage(i:i+subWindowSize-1,j: j+ subWindowSize -1));
            noiseImages(end+1,:) = reshape(currNoise,subWindowSize*subWindowSize,1);
            if (saveNoise)
                filename = sprintf('images_p2/signal/signal%06d.png',imgCnt);
                imwrite(currNoise ,filename)
            end
        end
    end
end
muN = mean(noiseImages);

% for i=1:size(noiseImages,1)
%     noiseImages(i,:) = noiseImages(i,:) - mu;
% end

Cn = cov(noiseImages);
[Un En] = eig(inv(Cn));
En(En<0)=0;
Dn=Un*sqrt(En)*Un';
whiteCn = cov(noiseImages*Dn);

%% Now build signal stats
muS = 0;
imgCnt = 0;
signalImages = [];
tumorSize = 34;
for k=1:9
    display(['Processing signal ' num2str(k)])
    currImage = tu{k};
%     figure(4);clf;imshow(currImage,[])
    currImage = currImage( floor((45-tumorSize)/2) + 1 : floor((45-tumorSize)/2) + tumorSize , floor((45-tumorSize)/2) + 1 :  floor((45-tumorSize)/2) + tumorSize);
    for i=1:subWindowSize:size(currImage,1)-subWindowSize
        for j=1:subWindowSize:size(currImage,2)-subWindowSize
            imgCnt= imgCnt+1;
            currSig = (currImage(i:i+subWindowSize-1,j: j+ subWindowSize -1));
            signalImages(end+1,:) = reshape(currSig,subWindowSize*subWindowSize,1);
            if (saveSignal)
                filename = sprintf('images_p2/signal/signal%06d.png',imgCnt);
                imwrite(currSig ,filename)
            end
        end
    end
end
muS = mean(signalImages);

% for i=1:size(noiseImages,1)
%     noiseImages(i,:) = noiseImages(i,:) - mu;
% end

Cs = cov(signalImages);
[Us Es] = eig(inv(Cs));
Es(Es<0)=0;
Ds=Us*sqrt(Es)*Us';
whiteCs = cov(signalImages*Ds);


%% Save
save('images_p2\noise\stats2.mat','Cn','windowSize','subWindowSize','muN','noiseImages', 'muS','Cs','signalImages')