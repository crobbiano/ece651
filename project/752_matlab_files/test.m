%# random data vector of integers
M = randi([50 200], [100 1]);

%# compute bins
nbins = 10;
binEdges = linspace(min(M),max(M),nbins+1);
aj = binEdges(1:end-1);     %# bins lower edge
bj = binEdges(2:end);       %# bins upper edge
cj = ( aj + bj ) ./ 2;      %# bins center

%# assign values to bins
[~,binIdx] = histc(M, [binEdges(1:end-1) Inf]);

%# count number of values in each bin
nj = accumarray(binIdx, 1, [nbins 1], @sum);

%# plot histogram
bar(cj,nj,'hist')
set(gca, 'XTick',binEdges, 'XLim',[binEdges(1) binEdges(end)])
xlabel('Bins'), ylabel('Counts'), title('histogram of measured bacterial')