function [erodMask,saveMask] = gcaMorphologicalOpeningWithGeometryConstraints(mask,varargin)
% 
% Small function that will mark which CCs from the erosion can be removed
% and NOT break the cell body these will be marked as potential higher
% confidence linkages.
% the CCs that should not be removed because they will create major
% breakages in the overal body
%%  INPUT: 
% mask (REQUIRED) RxC logical array
%       of the original local thresholded image before morphological opening
%       : where R is the height (ny) and C is the width (nx) of the input image
%
% %% PARAMS: Morphological Opening %%
%
%    'DiskSizeLarge' (PARAM) : Positive scalar 
%       
%       Default = 6 for LifeAct: 4 for membrane markers
%       This parameter specifies the radius of the disk (in pixels) to be used for
%       the removal of thin objects from the input mask (ie the filopodia)
%       Larger values will remove thicker structures. Note for the lifeAct
%       channels that have very strong filopodia signal very often these
%       gradients for crossing filopodia tend not to be well segmented so
%       disks use are typically slightly bigger than if using a membrane
%       marker where the filpodia often exhibit weak signal relative to
%       entire image and are those not segmented at all.
%
%
%    'DiskSizeSmall' (PARAM) : Positive scalar 
%     
%     Default = 3 : 
%     Morphological 
%OUTPUT:
% erodMask: 
% 
% saveMask: 
%% INPUTPARSER
%%Input check
ip = inputParser;

ip.CaseSensitive = false;
ip.KeepUnmatched = true; 
%REQUIRED
ip.addRequired('mask');

% PARAMETERS

ip.addParameter('DiskSizeLarge',6,@(x) isscalar(x));

ip.addParameter('considerGeometry',true); 
ip.addParameter('DiskSizeSmall',3,@(x) isscalar(x)); 

ip.parse(mask,varargin{:});
%% Initiate
[ny,nx] = size(mask);
saveMask = zeros(ny,nx);
%% Perform morphological opening
erodForBody = imopen(mask,strel('disk',ip.Results.DiskSizeLarge,0));

if ip.Results.considerGeometry
    %% Check local morphological opening will result in an increase in the number of total connected components
    % if so reduce the morphological opening radius for these regions.
    
    diffMask = mask-erodForBody;
    CCDiff = bwconncomp(diffMask);
    CCStart = bwconncomp(mask);
    numStart = CCStart.NumObjects;
    
    for iPiece = 1:CCDiff.NumObjects
        
        maskPiece = zeros(ny,nx);
        
        maskPiece(CCDiff.PixelIdxList{iPiece}) = 1;
        
        CCEnd = bwconncomp(mask-maskPiece);
        numEnd= CCEnd.NumObjects;
        delta =  numStart - numEnd;
        
        if delta < 0 ; % the number of end pieces  more than the number of start pieces
            CCDiff.diffMark{iPiece}= 1;
        else
            CCDiff.diffMark{iPiece}= 0;
        end
    end % for iPiece
    
    %% Add back the connected component regions marked above
    filter = vertcat(CCDiff.diffMark{:});
    filter = logical(filter);
    
    % erosion pieces to save
    saveMask(vertcat(CCDiff.PixelIdxList{filter})) = 1;
    
    erodMask = (erodForBody|saveMask);
    %% perfrom a smaller radius morphological opening to remove any filopodia like-structures along these pieces
    erodMask = imopen(erodMask,strel('disk',ip.Results.DiskSizeSmall,0));
else 
    erodMask = erodForBody; 
end % if ip.Results.considerGeometry
end % function 