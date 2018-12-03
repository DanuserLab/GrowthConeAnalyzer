function [ filoIdx,sign ] = getFiloIDFromCoords(filoInfo, x,y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
test  = arrayfun(@(x) x.Ext_coordsXY,filoInfo,'uniformoutput',0);
%  mask = zeros(size(imgSize));
%  mask(idx) = 1;
%  bwdist(mask)

d = cellfun(@(z) min(((z(:,1)-x).^2 + (z(:,2)-y).^2)),test);
%idx = sub2ind(imgSize,y,x);


filoIdx =  find(d == min(d));

if length(filoIdx) >1
    % test just the end point.
    filoC = filoInfo(filoIdx);
    test2 = arrayfun(@(x) x.Ext_endpointCoordFitXY,filoC,'uniformoutput',0);
    d2 = cellfun(@(z) min(((z(:,1)-x).^2 + (z(:,2)-y).^2)),test2);
    filoIdx = filoIdx(d2==min(d2));
    %      filoIdx = filoIdx(1);
end

% get sign 
dFinal = cellfun(@(z) (z(:,1)-x).^2 + (z(:,2)-y).^2,test(filoIdx),'uniformoutput',0); 
x = dFinal{:}; 
idx = find(x==min(x));
pixIndices = filoInfo(filoIdx).Ext_pixIndices;
idxEnd = find(pixIndices == filoInfo(filoIdx).Ext_endpointCoordFitPix);
% pixIndicesPlot = pixIndices(1:idxEnd);
testForSign = idx>=idxEnd;
if testForSign 
    sign = 1; 
else 
    sign = -1; 
end 
end

