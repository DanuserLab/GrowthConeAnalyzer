function [filoOrient,filoBranchNew] = GCAAnalysisExtract_filoOrient(filoBranch,filoFilterSet)
%% GCAAnalysisExtract_filoOrient
% Collects Filopodia Orientation Distributions using a
% filoFilterSet specified using GCACreateFilopodiaFilterSet.m
%
%INPUT:
%
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
%                       significant.
%                       Output of GCACreateFilopodiaFilterSet.m
%
% OUTPUT: filoOrient :  rx1 cell array where r (row) is the number of frames
%                    and each cell holds
%                    an rx1 double of the filopodia orientation for each
%                    filopodia
%                    where r is the number of filopodia in a given frame
%                    after the filoFilterSet filter is applied
%% Initiate
nFrames = numel(filoFilterSet); 
filoOrient = cell(nFrames,1);
filoBranchNew = []; 
%% START 
for iFrame = 1:nFrames
    filoInfo = filoBranch(iFrame).filoInfo;
    if ~isempty(filoInfo)
        filterFrameC= filoFilterSet{iFrame};
        filoInfoFilt = filoInfo(filterFrameC(:,1)) ;
        orient =  vertcat(filoInfoFilt.orientation);
        filoOrient{iFrame} = orient;
    else
        filoOrient{iFrame} = [];
    end
    clear orient
end
end % Function
