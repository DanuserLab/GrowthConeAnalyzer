function [ branchDensitiesCell,filoBranchNew] = GCAAnalysisExtract_filoDensityAlongBranch(filoBranch,filoFilterSet,varargin)
%%% GCAAnalysisExtract_filoDensityAlongBranch
% Collects filopodia density along a branch for an entire movie for a
% filopodia filter set
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
% OUTPUT:
% branchDensitiesCell:  rx1 cell array where r (row) is the number of frames
%                    and each cell holds
%                    an rx1 double of the feature: 
%                    number of braches per 10 um of 
%                    filopodia segment, where r is the number of branches in a given frame
%                    after the filoFilterSet filter is applied
%% Check Parameters
ip = inputParser;
ip.KeepUnmatched = true;

ip.CaseSensitive = false;
ip.addParameter('umPerPixel',.216);

ip.parse(varargin{:});
%% Initiate
nFrames = numel(filoFilterSet); 
branchDensitiesCell = cell(nFrames,1);
filoBranchNew = []; 
%% Start
for iFrame = 1:nFrames
    
    filoInfo = filoBranch(iFrame).filoInfo;
    
    filterFrameC= filoFilterSet{iFrame};
    if ~isempty(filterFrameC)
        filoInfoFilt  = filoInfo(filterFrameC(:,1)); % after this filter should have only N order branch and its corresponding N-1 branchstem
        %% get the branch stems
        types = vertcat(filoInfoFilt(:).type);
        % change type 0 (no filo) to type 1 as want to calculate density of
        % branching using ALL filo attached to veil.
        types(types==0) = 1;
        NTypes = unique(types);
        % NTypes(NTypes==0)= 1; % for now just switch the 0 order (no filo to 1)
        typeStem = min(NTypes);
        
        if length(NTypes)==2 % filter ok
            % get the stem lengths
            idxStem= vertcat(types)==typeStem;
            lengthStem = vertcat(filoInfoFilt(idxStem).Ext_length).*ip.Results.umPerPixel;
            
            filoInfoStem = filoInfoFilt(idxStem);
            
            % get the number of branches per stem - tricky part here is this need to
            % likewise be filtered by fit etc which it will not be in the length of the .conIdx.
            % get the number of filo
            IDsCurrentSet = find(filterFrameC);
            % for each stem get the conIdx and filter by filtInfo
            numFiloBeforeFilt = arrayfun(@(x) length(filoInfoStem(x).conIdx),1:sum(idxStem));
            numFilo = arrayfun(@(x) length(intersect(filoInfoStem(x).conIdx,IDsCurrentSet)),1:sum(idxStem)) ;
            
            density = numFilo'./lengthStem*10;
            %% Calculate Density
            branchDensitiesCell{iFrame,1} = density; % output
        else
            if ~(isempty(filoInfoFilt) || length(filoInfoFilt)==1);
                display('Check Branch Filter: N~=2');
            end
        end
    else
        branchDensitiesCell{iFrame,1} = [];
    end % if ~isempty
end % for iFrame
end % Function