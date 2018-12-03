function [densitiesCell,filoBranchNew] = GCAAnalysisExtract_filoDensityAlongVeil(filoBranch,filoFilterSet,veilStem,varargin)
% INPUT : 
%filoBranch :  REQUIRED:  large rx1 structure
%                       where r is the number of frames
%                       that includes all the information one would want about the
%                       segmentation of the movie (including the filopodia information)
%                       Output of Step VII of GCA : GCAfitFilopodiaMovie.m
%
%filoFilterSet: PARAM:  rx2 cell array where r (row) is the number of
%                       frames and each cell is a rBundlex2 logical vector.
%                       where rBundle is the number of Actin Bundles
%                       detected in the frame.
%                       Column 1 of rBundle places a 1 for each extra-Veil
%                       filopodia
%                       to be included in the analysis based on filtering
%                       criteria
%                       Column 2 of rBundle places a 1 for each embedded
%                       detection of the actin bundle that is considered
%                       significant based on filtering criteria.
%                       Output of GCACreateFilopodiaFilterSet.m
%OUTPUT: 
%densitiesCell: rx1 cell where r is the number frames of the movie and 
%               each cell holds a single value corresponding to the number 
%               of filopodia per 10 um veil for a given frame. 
%               Note cell array format is employed to maintain consistency with
%               other feature extractions maintain more than one value per
%               frame. 
%% Check Parameters
ip = inputParser;
ip.KeepUnmatched = true;
ip.CaseSensitive = false;
ip.addParameter('umPerPixel',.216);
ip.parse(varargin{:});
%% Initialize
nFrames = numel(filoFilterSet); 
densitiesCell = cell(nFrames,1);
distBoundMicronC = cell(length(filoBranch),1);
filoBranchNew = []; 
%% START
for iFrame = 1:nFrames
    % get filopodia
    
    filterFrameC= filoFilterSet{iFrame};
    
    % get masks
    neuriteMask = veilStem(iFrame).finalMask;
    
    % sort pixels
    [ny,nx] = size(neuriteMask);
    roiYX = bwboundaries(neuriteMask);
    edgeMask = zeros([ny,nx]);
    
    pixEdge =  sub2ind([ny,nx],roiYX{1}(:,1),roiYX{1}(:,2));
    edgeMask(pixEdge) = 1;
    
    % add a thinning step
    edgeMask = bwmorph(edgeMask,'thin','inf');
    
    % take out border pixels
    boundaryMask = zeros([ny,nx]);
    boundaryMask(1:ny,1) =1;
    boundaryMask(1:ny,nx)=1;
    boundaryMask(1,1:nx)= 1;
    boundaryMask(ny,1:nx) =1;
    
    edgeMask(boundaryMask==1) = 0;
    pixEdgeMask = find(edgeMask==1);
    EPs = getEndpoints(pixEdgeMask,[ny,nx],0);
    
    
    pixIdxBack = nan(length(pixEdgeMask),1);
    edgeMask = logical(edgeMask);
    
    if isempty(EPs)
        EPs = getEndpoints(pixEdgeMask(2:end),[ny,nx],1);
    end
    
    transform = bwdistgeodesic(edgeMask,EPs(1,1),EPs(1,2));
    
    iPix = 0;
    while length(find(transform==iPix)) == 1
        pixIdxBack(iPix+1) = find(transform==iPix); % start at the endpoint
        iPix = iPix +1;
    end
    
    distBoundMicron =  calculateDistance(pixIdxBack,[ny,nx],0,'pixelSizeMic',ip.Results.umPerPixel);
    
    if isnan(distBoundMicron) && iFrame ~=1
        distBoundMicron  = distBoundMicronC{iFrame-1};
    end
    %% Calculate Density
    if ~isempty(filterFrameC)
        numFilo = sum(filterFrameC(:,1));  
        densitiesCell{iFrame,1} = numFilo/distBoundMicron*10;
        distBoundMicronC{iFrame,1} = distBoundMicron;
    else
        densitiesCell{iFrame,1} = []; 
        distBoundMicronC{iFrame,1} = [];
    end   
end  % for iFrame
end % function