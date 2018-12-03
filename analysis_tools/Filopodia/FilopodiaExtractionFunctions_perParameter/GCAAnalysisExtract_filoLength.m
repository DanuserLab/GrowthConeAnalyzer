function [out,filoBranchNew] = GCAAnalysisExtract_filoLength(filoBranch,filoFilterSet,varargin)
%% GCAAnalysisExtract_filoLength
% Filters and Collects Filopodia/Actin Bundle Length Features from a filoBranch.mat
% for further analysis 
%%
%INPUT:
%
%filoBranch:  (REQUIRED): large rx1 structure
%                       where r is the number of frames
%                       that includes all the information one would want about the
%                       segmentation of the movie (including the filopodia information)
%                       Output of Step VII of GCA : GCAfitFilopodiaMovie.m
%
%filoFilterSet: (REQUIRED):  rx2 cell array where r (row) is the number of
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
%umPerPixel:  (PARAM): scalar : as name suggests provides the pixelSize in
%                      um (Note this function always converts the lengths from 
%                      pixels to microns)
%                      (DEFAULT: 0.216)
%
%filoPart:    (PARAM): character 'Int_', 'Ext_','Tot'
%                      specifying the part of the actin bundle "filopodia"
%                      to use for calculation 
%                      (DEFAULT: 'Ext_')
%                     
%outPercent:  (PARAM): logical 
%                      (DEFAULT: false) 
%                      if true and filoPart = 'Tot' will
%                      output the percent each filopodia is embedded in
%                      veil, if false and filoPart = 'Tot' will output 
%                      the total length of the actin bundle. 
%                      if filoPart does not equal 'Tot' this parameter will
%                      be ignored. 
%
%filterZeroPercent : PARAM: logical
%                           (DEFAULT: true)
%                            if true will filter out any actin bundles that
%                            with 0 % veil embedment. 
%
% OUTPUT:
% out:  rx1 cell array where r (row) is the number of frames
%                    and each cell holds the distribution of filo lengths
%                    measurements (in microns) for a given frame. 
%                    which measurement output is dictated by the params
%                    'filoPart' and 'outPercent'
%
% filoBranchNew: some of the GCAAnalysisExtract_xFeature.m functions add new information to the
%                filoBranch structure if necessary. For this specific function
%                output is empty, but is included to maintain consistency
%                with these other functions. 
%% Check Parameters
ip = inputParser;
ip.KeepUnmatched = true;
ip.CaseSensitive = false;
ip.addParameter('filoPart', 'Ext_');
ip.addParameter('outPercent',false); % in the case of filoPart
% 'Tot' will output the percentage embedded for each filopodia
% instead of the total length of the actin bundle
ip.addParameter('filterZeroPercent',true);
ip.addParameter('umPerPixel',.216);

ip.parse(varargin{:});
%% Initiate
nFrames = numel(filoFilterSet); 
out = cell(nFrames,1);
filoBranchNew = []; % keep to maintain consistency with some of hte other functions
%% Extract
for iFrame = 1:nFrames
    
    filoInfo = filoBranch(iFrame).filoInfo;
    
    if ~isempty(filoInfo);
        % currently 2 columns of the filter - one for the external and one for the
        % internal based on fitting criteria.
        filterFrameAll= filoFilterSet{iFrame};
        
        if strcmpi(ip.Results.filoPart,'Int_');
            filterFrameC = (filterFrameAll(:,1) == 1 & filterFrameAll(:,2) == 1);
        else
            filterFrameC = filterFrameAll(:,1);
        end
        
        filoInfoFilt  = filoInfo(filterFrameC);
        
        % collect lengths: if just int or ext just use the respective filopodia
        % filter...
        if ~strcmpi(ip.Results.filoPart,'Tot');
            lengths =  vertcat(filoInfoFilt(:).([(ip.Results.filoPart) 'length'])).*ip.Results.umPerPixel;
        else
      
            filterInt = (filterFrameAll(:,1) == 1 & filterFrameAll(:,2) ==0 ); % get the ID of all non-fits internally this filter is the length of
            %             % the original filoInfo detection
            filterInt = filterInt(filterFrameC); % filter the above logical filter to make it the same
            %filterInt = filterFrameC(:,2) ==0;  % get all the embedded that do not pass the filter criteria
            
            lengthsInt =  vertcat(filoInfoFilt(:).Int_length).*ip.Results.umPerPixel;
            lengthsExt = vertcat(filoInfoFilt(:).Ext_length).*ip.Results.umPerPixel;
            % convert NaN lengths of internal to zero
            lengthsInt(isnan(lengthsInt))=0;
            lengthsInt(filterInt) = 0;
            
            lengths = lengthsInt + lengthsExt;
            if ip.Results.outPercent;
                percent = lengthsInt./lengths;
                if ip.Results.filterZeroPercent
                    percent = percent(percent~=0);
                end
            end
        end
        if (ip.Results.outPercent && strcmpi(ip.Results.filoPart,'Tot'));
            out{iFrame} = percent;
        else
            out{iFrame} = lengths;
        end
    else
        out{iFrame} = [];
    end % isempty
    clear lengths
end % iFrame
end % Function