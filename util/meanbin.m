function [xmean, ymean, ystd, nPoints] = meanbin(x, y)

    bins = unique(x);
    nBins = length(bins);
    [ymean,ystd, xmean, nPoints] = deal(zeros(nBins,1));
    for i=1:nBins
        inds = x==bins(i);
        ymean(i) =nanmean(y(inds));
        ystd(i)= nanstd(y(inds));
        xmean(i) = nanmean(x(inds));
        nPoints(i)= sum(inds);
      %  close all;figure; hist(y(inds), 20);pause;
    end
end