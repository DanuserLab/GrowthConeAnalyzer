function [filoAvgInt,filoBranchNew] = GCAAnalysisExtract_filoAvgIntensity(filoBranch,filoFilterSet,varargin)
%% GCAAnalysisExtract_filoAvgIntensity
% Collects Filopodia Intensity Distributions for a
% Filtered Set of Filopodia 
%%
%INPUT:
%
%filoBranch :  REQUIRED:  large rx1 structure
%                       where r is the number of frames
%                       that includes all the information one would want about the
%                       segmentation of the movie (including the filopodia information)
%                       Output of Step VII of GCA : GCAfitFilopodiaMovie.m
%
%filoFilterSet: REQUIRED:  rx2 cell array where r (row) is the number of
%                       frames and each cell is a rBundlex2 logical vector.
%                       where rBundle is the number of Actin Bundles
%                       detected in the frame.
%                       Column 1 of rBundle places a 1 for each extra-Veil
%                       filopodia to be included in the analysis based on filtering
%                       criteria
%                       Column 2 of rBundle places a 1 for each embedded
%                       detection of the actin bundle that is considered
%                       significant based on filtering criteria.
%                       Output of GCACreateFilopodiaFilterSet.m
%
%% Check Input
ip = inputParser;
ip.KeepUnmatched = true;

ip.CaseSensitive = false;

ip.addParameter('filoPart', 'Ext_'); 
ip.addParameter('normToVeil',true); 

ip.parse(varargin{:});
%% Initiate
if ip.Results.normToVeil
    field = 'IntensityNormToVeil'; 
else 
    field = 'Intensity'; 
end 

nFrames = numel(filoFilterSet); 

filoAvgInt = cell(nFrames,1);
filoBranchNew = []; 
%% START
for iFrame = 1:nFrames
    filoInfo = filoBranch(iFrame).filoInfo;
    if ~isempty(filoInfo)
        filterFrameC= filoFilterSet{iFrame};
        
        if strcmpi(ip.Results.filoPart,'Int_')
            filterFrameC = (filterFrameC(:,1)==1 & filterFrameC(:,2) ==1);
        end
        
        filoInfoFilt = filoInfo(filterFrameC(:,1)) ;
        
        intensity =  vertcat(filoInfoFilt(:).([(ip.Results.filoPart) field]));
        
        filoAvgInt{iFrame} = intensity;
    else
        filoAvgInt{iFrame} = [];
    end
    clear intensity
end % for 
end % function 