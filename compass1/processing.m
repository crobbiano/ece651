%% Process boob shots
clear all
clc

%% Read in images
cd('images\cropped')
mats = dir('*.png');0
for i=1:length(mats)
    im{i} = imread(mats(i).name);
%     filename = sprintf('Image%02d.png',i)
%     imwrite(im(i).Image ,filename)
end
cd('..\..')

%%
windowSize = 45
mu = 0;
imgCnt = 0;
noiseImages = [];
for k=1:10
    currImage = im{k};
    for i=1:20:size(currImage,1)-windowSize
        for j=1:20:size(currImage,2)-windowSize
            imgCnt= imgCnt+1;
            noiseImages(end+1,:) = reshape(double(currImage(i:i+windowSize-1,j: j+ windowSize -1)),windowSize*windowSize,1);
        end
    end
end
mu = mean(noiseImages);

% for i=1:size(noiseImages,1)
%     noiseImages(i,:) = noiseImages(i,:) - mu;
% end

C = cov(noiseImages);
[U E] = eig(inv(C));
E(E<0)=0;
D=U*sqrt(E)*U';
whiteC = cov(noiseImages*D);
%% Setup operations and operate
makeNoise = 0;
if (makeNoise ==1 )
    operations = [1 891 477 30;...
        1 856 500 26;...
        2 599 525 33;...
        3 566 471 40;...
        4 343 538 29;...
        5 618 462 44;...
        6 875 432 20;...
        7 530 680 20;...
        8 265 470 29;...
        9 388 347 26;...
        10 671 337 45;];
    
    
    buff = ceil(max(operations(:,4))/2);
    mu = zeros(windowSize*windowSize,1);
    numValidNoise = 0;
    for i=1:length(operations)
        currIm = im(operations(i,1)).Image;
        row = operations(i,2);
        column = operations(i,3);
        tumor{i}=currIm(row - buff+1:row + buff-1,column - buff+1: column + buff-1);
        %     filename = sprintf('images/tumors/tumor%d.mat',i);
        currTumor = tumor{i};
        %     save(filename, 'currTumor');
        %
        %     filename = sprintf('images/tumors/tumor%d.png',i);
        %     imwrite(currTumor, filename);
        
        noiseIms{(i-1)*8 + 1} = currIm(row - 3*buff+1:row - buff-1,column - buff+1: column + buff-1);
        noiseIms{(i-1)*8 + 2} = currIm(row + buff+1:row + 3*buff-1,column - buff+1: column + buff-1);
        noiseIms{(i-1)*8 + 3} = currIm(row - buff+1:row + buff-1,column - 3*buff+1: column - buff-1);
        noiseIms{(i-1)*8 + 4} = currIm(row - buff+1:row + buff-1,column + buff+1: column + 3*buff-1);
        noiseIms{(i-1)*8 + 5} = currIm(row - 3*buff+1:row - buff-1,column - 3*buff+1: column - buff-1);
        noiseIms{(i-1)*8 + 6} = currIm(row - 3*buff+1:row - buff-1,column + buff+1: column + 3*buff-1);
        noiseIms{(i-1)*8 + 7} = currIm(row + buff+1:row + 3*buff-1,column - 3*buff+1: column - buff-1);
        noiseIms{(i-1)*8 + 8} = currIm(row + buff+1:row + 3*buff-1,column + buff+1: column + 3*buff-1);
        %     figure(54);clf
        for k=1:8
            %         subplot(4,2,k);
            %         imshow(noiseIms{(i-1)*8 + k})
            someval = sum(sum(double(noiseIms{(i-1)*8 + k})))
            if (someval > 20000)
                numValidNoise = numValidNoise + 1;
                mu = (mu + reshape(double(noiseIms{(i-1)*8 + k}),windowSize*windowSize,1));
                %             meanImage = (meanImage + reshape(noiseIms{(i-1)*8 + k},windowSize*windowSize,1));
            end
            %         filename = sprintf('images/noise/noise%d.mat',(i-1)*8 + k);
            %         currNoise = noiseIms{(i-1)*8 + k};
            %         save(filename, 'currNoise');
            %
            %         filename = sprintf('images/noise/noise%d.png',(i-1)*8 + k);
            %         imwrite(currNoise, filename);
        end
        
    end
    mu = mu/numValidNoise;
end
%% Save
save('images\noise\stats.mat','C','windowSize','mu','noiseImages')