function [ backMask,beforeDilate] = GCACreateBackgroundMask(veilStemMask,filoInfo,img,varargin )
%GCACreateBackgroundMask 
%%   Check Input 
ip = inputParser;
ip.CaseSensitive = false;

ip.addRequired('veilStemMask'); 
ip.addRequired('filoInfo'); 
ip.addRequired('img'); 

ip.addParameter('dilateLocalRegion',true); % flag to dilate 
ip.addParameter('LRDilRad',10); % dilation radius in pixels. 


ip.CaseSensitive = false;
ip.parse(veilStemMask,filoInfo,img,varargin{:});
%%
imgSize = size(img); 
% collect all the xyCoords of the filoBranch reconstruction before
% fitting and put into a mask.
xyCoords = vertcat(filoInfo(:).Ext_coordsXY);
linIndices = sub2ind(imgSize,xyCoords(:,2),xyCoords(:,1));
linIndices = linIndices(~isnan(linIndices)); 
filoMask = zeros(imgSize);
filoMask(linIndices) = 1;
gcaMask= (filoMask|veilStemMask);
% dilate gcaMask: walking forward from initial detections/the interpolation 
% via the reconstruction likewise helps
% make sure no filopodia are included in the background estimation. 
% dilate around this estimate of the GC after reconstruction
maskObjGC = imdilate(gcaMask,strel('disk',ip.Results.LRDilRad));

% Now estimate the background (in case any high fluorescence background)
% thresh = mean + 2*std of the first pdf, plus some local dilation
% operations if user desires. 
% The below are potentially the same masks used in the first filtering step of 
% GCAReconstructFiloBranch so can maybe have an option to just save those 
% .tifs at that time and read in. 
if ip.Results.dilateLocalRegion
    [maskBackEst,~,~] = gcaEstimateBackgroundArea(img,'PostProcess',true);
    % dilate around this estimate of the object based on the intensity 
    % thresholding
    maskObjAll = imdilate(~maskBackEst,strel('disk',ip.Results.LRDilRad));
else
    maskObjAll = gcaEstimateBackgroundArea(img);
end

% combine the two masks: little bit overkill but should ensure no filo in
% background AND nor any other high intensity objects 
backMask = ~(maskObjGC | maskObjAll);
backMask = logical(backMask); % make sure logical

beforeDilate = (gcaMask | ~maskBackEst); 

end

