function [ filoFilterSet,filterParams] = GCACreateFilopodiaFilterSet(filoBranch,varargin)
%GCACreateFilopodiaFilterSet:
% Create a logical filterset for filopodia analysis based on a number of
% user adaptable criteria.
%
% INPUT:
%      filoBranch:     REQUIRED: rx1 structure
%                                where r is the number of frames in the movie 
%                                output by Step VII of GCAfitFilopodiaMovie
%                                Includes a field called .filoInfo with the
%                                necessary information for filo filtering. 
%
%
%      filterType:   PARAM:    char   Currently Either
%      'ConnectVeil_LengthInt',
%      'ConnectVeil_DensityOrient',
%      'Branch2ndOrder_LengthInt'
%
%
% OUTPUT:
%       filoFilterSet:   rx2 cell of logical filters for the filopodia for each frame
%                        where r1 is the number of frames and column 1 is
%                        the traditional non-embedded filter and column2 is
%                        the embedded actin bundle filter.
%% Check Input
ip = inputParser;
ip.CaseSensitive = false;
ip.addParameter('pixelSizeNm',216,@(x) isscalar(x));
ip.addParameter('filterType',[]); % if empty input the filopodia filter parameters 
% directory: other wise flag to load pre-defined filters for analysis 
% Pre-defined filter choices: 
% 'ConnectToVeil_LengthInt' 
% 'ConnectToVeil_Density', 
% 'Branch2ndOrder_LengthInt'
% 'Branch2ndOrder_Density_WithZeroOrder'
% 'Branch2ndOrder_Density_NoZero'
% 'Validation'
% 'Validation_NoEmbed'

% If 'filterType' empty will use below input 
ip.addParameter('filoTypes',[0 Inf]); % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
% default [0-Inf] Plot all segments that pass fit criteria 

ip.addParameter('filterByBundleLength',[0.3,inf]); % in microns [minValue, maxValue]
ip.addParameter('saveFiloByLengthAndSig',[]);

ip.addParameter('filterByFit',true);
ip.addParameter('embedFitCriteria', 95);
ip.addParameter('filoFitCriteria', 95);

ip.addParameter('filterBasedOnGroupUpstream',0);
ip.addParameter('filterBasedOnBranchConnection',0);
ip.addParameter('filterIntNoConnect',true);


ip.addParameter('trunc',false); % quickfix 201801
ip.parse(varargin{:});
%% START FILTERING

endFrame = length(filoBranch);
if ip.Results.trunc 
    if length(filoBranch)==1 
        endFrame =1; 
    else 
    endFrame = endFrame -1 ; % quick fix for my botch up 201801
    end 
end 

%% Load the predefined filter parameters if applicable 
% if ip.Restuls.filterType is empty, use the user defined input 
if isempty(ip.Results.filterType) 
    filterParams.filoTypes = ip.Results.filoTypes; 
    filterParams.filterByBundleLength = ip.Results.filterByBundleLength; 
    filterParams.saveFiloByLengthAndSig = ip.Results.saveFiloByLengthAndSig; 
    
    filterParams.filterByFit = ip.Results.filterByFit;
    filterParams.filoFitCriteria = ip.Results.filoFitCriteria; 
    filterParams.embedFitCriteria = ip.Results.embedFitCriteria; 
    
    filterParams.filterBasedOnGroupUpstream = ip.Results.filterBasedOnGroupUpstream; 
    filterParams.filterBasedOnBranchConnection = ip.Results.filterBasedOnBranchConnection;
    filterParams.filterIntNoConnect = ip.Results.filterIntNoConnect; 

else % load the pre-defined filter settings using filter name
  filterParams=  GCAloadDefaultFilopodiaFilters(ip.Results.filterType);   
end 
%% 

for iFrame = 1:endFrame
    filoInfo = filoBranch(iFrame).filoInfo;
    
    if ~isempty(filoInfo);
        %% FILTER BY EMBEDDED
        
        if isfield(filterParams,'embedFitCriteria') && ~isempty(ip.Results.embedFitCriteria)
            p1Int  = arrayfun(@(x) filoInfo(x).Int_params(1),1:length(filoInfo),'uniformoutput',0)'; % get the amplitude of the fit
            resids = arrayfun(@(x) filoInfo(x).Int_resid,1:length(filoInfo),'uniformoutput',0)'; % get the residuals of the fit
            
            % Filter out the amplitude of signal that is less than the
            % 95th percentile of the residuals.
            filtInt  = cellfun(@(x,y) (y < prctile(x,filterParams.embedFitCriteria) | isnan(y)),resids,p1Int);
        else
            filtInt = false(length(filoInfo),1); % no filtering based on embedded
        end
        %% FILTER BY SIGMOIDAL FITTING
        if filterParams.filterByFit;
            
            % find those data where the fit produces a NaN
            endpointCoordTest = vertcat(filoInfo(:).('Ext_endpointCoordFitPix')); % want this to filter out internal
            % put NaN for each filoInfo if the length metric is NaN.
            
            nanToRemove = isnan(endpointCoordTest);
            
            % find those data which did not pass the exit flag
            test = vertcat(filoInfo(:).Ext_exitFlag) ; % check exit flag - less than 1 means a majure failure
            % very infrequent though - can happen for a crossing check to see
            % what type of filopodia typical look like this.
            
            %logicals for filtering of filopodia
            toKeepBasedOnExit = (test>0 ); % good fits
            
            %filterParams.filoFitCriteria
            p1Ext  = arrayfun(@(x) filoInfo(x).Ext_params(1),1:length(filoInfo),'uniformoutput',0)'; % get the amplitude of the fit
            residsExt = arrayfun(@(x) filoInfo(x).Ext_resid,1:length(filoInfo),'uniformoutput',0)';
            
            %filtInt = (p1Int<filterParams.embedFitCriteria(1) | isnan(p1Int));
            % for a quick test filter out the amplitude of signal that is less than the
            % 95th percentile of the residuals.
            filtExt  = cellfun(@(x,y) (y < prctile(x,filterParams.filoFitCriteria) | isnan(y)),residsExt,p1Ext);
            toKeepBasedOnExit =  toKeepBasedOnExit & ~filtExt;
            % eventually here want to create based on
        end % if filterByFit
        
        %% FILTER BY FILO TYPE
        
        type = vertcat(filoInfo(:).('type'));
        if length(filterParams.filoTypes)>1 && filterParams.filoTypes(2) == inf
            % find the max branch type for that frame
            filoTypesC = filterParams.filoTypes(1):max(type);
        else
            filoTypesC = filterParams.filoTypes;
        end
        % make cell array of logicals
        toKeepBasedOnTypeCell = arrayfun(@(x) type == filoTypesC(x) ,1:length(filoTypesC),'uniformoutput',0);
        
        % combine them
        toKeepBasedOnType = horzcat(toKeepBasedOnTypeCell{:});
        toKeepBasedOnType = sum(toKeepBasedOnType,2)~=0;
        %% FILTER BY TOTAL LENGTH
        if ~isempty(filterParams.filterByBundleLength)
            minLength = filterParams.filterByBundleLength(1);
            maxLength = filterParams.filterByBundleLength(2);
            % get all lengths
            lengthExt = vertcat(filoInfo(:).Ext_length);
            if isfield(filterParams,'embedFitCriteria') && ~isempty(filterParams.embedFitCriteria)
                lengthInt = vertcat(filoInfo(:).Int_length);
            else
                lengthInt = zeros(length(filoInfo),1);
            end
            % change NaNs into zero for addition
            lengthExt(isnan(lengthExt))= 0 ;
            lengthInt(isnan(lengthInt)) = 0 ;
            lengthInt(filtInt) = 0;
            lengthBundle = lengthExt + lengthInt;
            % convert to um
            lengthBundle = lengthBundle.* ip.Results.pixelSizeNm./1000; % most intuitive to set the length of the bundle in um
            
            toKeepBasedOnLength = lengthBundle>minLength & lengthBundle<maxLength; % logical marking all the filo that
            % make the bundle length criteria
        end  % isempty
        %% Save Long Filopodia with Strong Signal If Desired
        % note sometimes filopodia that cross will have a poor fit, these need to be
        % excluded from certain measurements (for instance length) and need
        % to be included in others (for instance in a density metric).
        
        if ~ isempty(filterParams.saveFiloByLengthAndSig);
            s = filterParams.saveFiloByLengthAndSig;
            % get filopodia that meet length cut-off..
            lengthExt = vertcat(filoInfo(:).Ext_length)*ip.Results.pixelSizeNm./1000; % convert
            
            savePop1 = lengthExt>s(1,1) & lengthExt < s(1,2) & toKeepBasedOnType ; %
            
            % get the full population
            intensities = vertcat(filoInfo(:).Ext_IntensityNormToVeil);
            intensitiesForPer = intensities(~isnan(intensities));
            cutoffMin = prctile(intensitiesForPer,s(2,1));
            cutoffMax = prctile(intensitiesForPer,s(2,2));
            
            savePop2 = intensities>cutoffMin & intensities<cutoffMax & ~isnan(intensities) & toKeepBasedOnType;
            
            savePop = savePop1 & savePop2;
        end
        
        %% Particle Filter
        if isfield(filterParams,'minPercentFiloAboveBack') % Quick fix for now to make sure does not crash
            if ~isempty(filterParams.minPercentFiloAboveBack)
                existMaxPosSlope = arrayfun(@(x) ~isempty(find(filoInfo(x).Ext_slopeMaxPos(:,1)<filoInfo(x).Ext_length)),1:length(filoInfo));
                % get the mean of hte local background
                meanLocalBack = arrayfun(@(x) mean(x.Ext_localBackValues),filoInfo);
                stdLocalBack = arrayfun(@(x) std(x.Ext_localBackValues),filoInfo);
                aboveMeanBack = arrayfun(@(x) filoInfo(x).Ext_weightedAvg>(meanLocalBack(x)+filterParams.nSTDParticleDetector*stdLocalBack(x)),1:length(filoInfo),'uniformoutput',0);
                % 
%                 aboveMeanBack2 = arrayfun(@(x) filoInfo(x).Ext_weightedAvg>meanLocalBack(x)+filterParams.nSTDParticleDetector*3,1:length(filoInfo),'uniformoutput',0); 
                disFiloIdx = arrayfun(@(x) x.Ext_distFilo<x.Ext_length,filoInfo,'uniformoutput',0);
                
%                 v2 = cellfun(@(x,y) x(y), aboveMeanBack2,distFiloIdx,'uniformoutput',0); 
%                 percentFiloAboveBack2 = cellfun(@(x) (sum(x)/length(x))*100,v2); 
                
                v= cellfun(@(x,y) x(y),aboveMeanBack,disFiloIdx,'uniformoutput',0);
                percentFiloAboveBack = cellfun(@(x) (sum(x)/length(x))*100,v);
             
                
                % if the first point around the connection is way above the
                % local background it is a good indication that it's a long
                % filo with decay into the background NOT a particle. 
                % save these from the particle detector. % could
                % technically make it such that if the first point is zero
                % for this up the std filter 
                firstPtZero = arrayfun(@(x) filoInfo(x).Ext_weightedAvg(2) >  (meanLocalBack(x)+9*stdLocalBack(x)),1:length(filoInfo)); 
                
               %  toDeleteParticleTest = (percentFiloAboveBack <
               %  filterParams.minPercentFiloAboveBack)' &
               %  existMaxPosSlope'; % original 
                toDeleteParticleTest = (percentFiloAboveBack < filterParams.minPercentFiloAboveBack)' & existMaxPosSlope' & ~firstPtZero'; % 
                %toDeleteParticleTest2 = ~firstPtZero' & (percentFiloAboveBack2 < filterParams.minPercentFiloAboveBack); 
               
                toKeepParticleTest = ~toDeleteParticleTest;
            else
                toKeepParticleTest = ones(length(filoInfo),1);
            end % filterParams.particleFilter
        else
            toKeepParticleTest = ones(length(filoInfo),1);
        end % particle filter
       %% Filter Based on Group Upstream
        toKeepBasedOnGroupUpstream = ones(length(filoInfo),1);
        if filterParams.filterBasedOnGroupUpstream
            % get all the first order filopodia filtered (type 1)
            % if filtered propogate the filter to all downstream branches by getting
            % the idx of the group ID.
            types = vertcat(filoInfo(:).type);
            typeID = unique(types);
            typeID(typeID==0) = [];  % remove all nonbranch filo
            %         typeID(typeID~=0) =[];
            nTypes= length(typeID)-1; % don't need to look at the last type as there is
            % no branches attached to this.
            %      for i = 2:nTypes
            if ~isempty(typeID)
                for i = 1:nTypes
                    idxBranchStruct= vertcat(filoInfo(:).type) == typeID(i);
                    if isfield(filterParams,'minPercentFiloAboveBack')
                        
                        idxBadFitWithAttach = ((~toKeepBasedOnExit | nanToRemove| ~toKeepParticleTest) & idxBranchStruct);
                    else
                        
                        idxBadFitWithAttach = ((~toKeepBasedOnExit | nanToRemove) & idxBranchStruct);
                    end
                    
                    % get the groupIDs of all the body attached filopodia
                    groupIDsToFilter = vertcat(filoInfo(idxBadFitWithAttach).groupCount) ;
                    groupIDsToFilter = unique(groupIDsToFilter);
                    
                    if ~isempty(groupIDsToFilter)
                        % want to filter out any filo branches with a nType greater
                        % than the low confidence branch previously filtered
                        % so for each  detection of the current frame
                        % ask if it is member of the groupIDstoFilter and if
                        % it has a ID greater than the branch stem with attachment
                        % under question (don't want to filter high confidence
                        % members of the group that are closer to the cell
                        % veil/stem)
                        toKeepBasedOnUpstreamNType(:,i) = arrayfun(@(x) ~(ismember(filoInfo(x).groupCount, groupIDsToFilter)...
                            & filoInfo(x).type > typeID(i)),1:length(filoInfo));
                    else
                        toKeepBasedOnUpstreamNType(:,i) = ones(length(filoInfo),1); % keep them all
                    end
                end % i = 1:nTypes
                
                % if there is no flag to remove keep
                toKeepBasedOnGroupUpstream = (sum(toKeepBasedOnUpstreamNType,2) == nTypes);
            else % no branches to check
                toKeepBasedOnGroupUpstream = ones(length(filoInfo),1);
            end
            clear toKeepBasedOnUpstreamNType
        end % filterBasedOnGroupUpstream
        %% Test for Connection after filter
        idxKeepBasedOnBranchConnection = ones(length(filoInfo),1);
        if filterParams.filterBasedOnBranchConnection
            % get the branch IDs
            
            % for each connectivity index test if conXYCoords (3) is significantly > than
            % length if it is get the corresponding conIdx
            
            % initiate logicalfilter vector: assume initially keep all filo
            idxKeepBasedOnBranchConnection = ones(length(filoInfo),1);
            
            if isfield(filoInfo,'conXYCoords') % if no branches in the frame will not have this field.
                
                for i = 1:length(filoInfo)
                    if ~isempty(filoInfo(i).conXYCoords)
                        % ask if the attachment distance is < the length metric
                        % if it is keep if not it will be set to zero
                        d = filoInfo(i).conXYCoords(:,3);
                        all = filoInfo(i).Ext_params;
                        conID = filoInfo(i).conIdx;
                        typeC = filoInfo(i).type;
                        
                        
                        if isnan(filoInfo(i).Ext_distFilo);
                            d = inf.*ones(length(d),1); 
                        else 
                        d= filoInfo(i).Ext_distFilo(d); 
                     
                        end 
                        if isnan(all)
                            l = 0;
                        else
                            l = filoInfo(i).Ext_params(2);
                        end
                        % get branches to remove
                        idxRemovePerFilo = d>(l+1);
                     
                        % this is poorly coded 
                        % 
                        branchIDs = []; 
                        if sum(idxRemovePerFilo)~=0
                            branchIDs = conID(idxRemovePerFilo); 
                           % for iBranch = 1:length(branchIDs)
                           downStreamC= vertcat(filoInfo(branchIDs).conIdx);
                            count = 1;
                            while count > 0
                                % get the IDs of all those branches connected downstream
                                % to the branch youi will be removing.
                               
                                % if empty nothing attached so break loop
                                if isempty(downStreamC)
                                  if count ==1 
                                      downStream = []; 
                                  end 
                                    count= 0;
                                else % check next
                                    downStream{count,1} = downStreamC; 
                                    count = count +1;
                                    downStreamC= vertcat(filoInfo(downStreamC).conIdx);
                                end % isempty
 
                            end % while
                            if ~isempty(downStream)
                                downStreamAll = vertcat(downStream{:});
                            else
                                downStreamAll = [];
                            end
                            
                            % get decendents
%                             groupC = filoInfo(conID(idxRemovePerFilo)).groupCount;
%                             groupIDs = vertcat(filoInfo(:).groupCount);
%                             typeAll = vertcat(filoInfo(:).type);
%                            
%                             idxRemoveDownstream = find(arrayfun(@(x) ismember(filoInfo(x).groupCount,groupC)...
%                                 & filoInfo(x).type>typeC+1,1:length(filoInfo)));
%                             
%                             toRemove{i} = idxRemoveDownstream';
                        else 
                            downStreamAll = []; 
                        end
                        toRemove{i,1} = [downStreamAll; branchIDs]; 
                    else toRemove{i,1} =[];
                    end
                end % i = 1:length(filoInfo)
                idxNumToRemoveAll = vertcat(toRemove{:});
                if ~isempty(idxNumToRemoveAll)
                    idxKeepBasedOnBranchConnection(idxNumToRemoveAll) = 0;
                end
                
                clear toRemove
            end % isfield(filoInfo,'conXYCoords')
        end % filterBasedOnBranchConnection
        %% Quick fix for embedded filo that hit a cell edge before being reconnected with partner.
        idxKeepBasedOnIntNoConnect =  ones(length(filoInfo),1);
        if isfield(filterParams,'filterIntNoConnect')
            if filterParams.filterIntNoConnect
                % get the end points of the two coords
                % if not within 3 pixels remove
                %         idxKeepBasedOnIntNoConnect = ones(length(filoInfo),1);
                y = zeros(length(filoInfo(:)),1);
                for i = 1:length(filoInfo);
                    if ~isnan(filoInfo(i).Int_coordsXY);
                        test1 = filoInfo(i).Int_coordsXY(1,:);
                        test2 = filoInfo(i).Ext_coordsXY(1,:);
                        if ~isnan(test1)
                            y(i) = sqrt((test1(1)-test2(1))^2 +(test1(2)-test2(2))^2);
                            
                        end
                    end
                end
                idxKeepBasedOnIntNoConnect = y<6;
                %          idxKeepBasedOnIntNoConnect  = idxKeepBasedOnIntNoConnect';
            end % if filterParams.filterIntNoConnect
        end % if isfield
        %% Add Biosensor Filter
        % Initiate
        toKeepBasedOnNumSigAtActinBundEnd = ones(length(filoInfo),1);
        toKeepBasedOnNumBackSampleSize = ones(length(filoInfo),1);
        if isfield(filterParams,'sigmaNumerator');
            % ok note there is a bit of a design problem as I had these as variable names
            % when collecting - for now they should be the same though for all
            % so keep here inflexible.
            
            % collect the background numerator values and estimate a local
            % threshold based on a simple 3*std of the background gaussian.
            backValuesNum = arrayfun(@(x) filoInfo(x).FRET_Detect.BackgroundValues,1:length(filoInfo),'uniformoutput',0) ;
            
            thresholds = cellfun(@(x) nanmean(x) + filterParams.sigmaNumerator*nanstd(x),backValuesNum);
            
            % for now just filter by the estimation at the end of the actin
            % bundle - make sure to flip the dimensions as saved the indices
            % of the full actin bundle fit from outside veil (traditional filo) to inside
            % veil. ie The first point should be the lowest.
            filoValuesNum =  arrayfun(@(x) flip(filoInfo(x).FRET_Detect.Ext_values),1:length(filoInfo),'uniformoutput',0);
            
            % unfortunately have to do this via a formal loop as NaN as a
            % index will error
            for i = 1:length(filoInfo)
                idxC = filoInfo(i).indicesFitFullBundle(1); % get the discrete index for the
                % end of the full actin bundle (ordered from external to
                % internal)
                if ~isnan(idxC)
                    % get the intensity at the end of the actin bundle and
                    % test if it is higher than the threshold
                    filoValuesNumAtPtC = filoValuesNum{i}(idxC);
                    toKeepBasedOnNumSigAtActinBundEnd(i) = filoValuesNumAtPtC > thresholds(i);
                else
                    toKeepBasedOnNumSigAtActinBundEnd(i) = 0;
                end
            end % for i = 1:length(filoInfo)
            toKeepBasedOnNumBackSampleSize = cellfun(@(x) length(x) > filterParams.backSampleSize,backValuesNum);
            toKeepBasedOnNumBackSampleSize =   toKeepBasedOnNumBackSampleSize';
        end
     
        %% Make Final Filo Filter Set Based on All the Above Criteria
        filoFilter = (toKeepBasedOnExit & toKeepBasedOnType & ~nanToRemove & toKeepBasedOnLength & toKeepBasedOnGroupUpstream & idxKeepBasedOnBranchConnection ...
            &  toKeepBasedOnNumSigAtActinBundEnd  & toKeepBasedOnNumBackSampleSize & toKeepParticleTest);
        if ~isempty(filterParams.saveFiloByLengthAndSig);
            if isfield(filoInfo,'minPercentFiloAboveBack');
                if ~isempty(filterParams.minPercentFiloAboveBack)
                    savePop  = savePop & toKeepParticleTest; % don't save if didn't make particle filter criteria
                end
            end
            filoFilter = (filoFilter | savePop);
        end
        filtInt = ~filtInt & idxKeepBasedOnIntNoConnect;
        %
        filoFilterSet{iFrame} = [filoFilter filtInt];
    else
        filoFilterSet{iFrame} = [];
        display(['No FiloInfo for frame ' num2str(iFrame)]);
    end % ~isempty
end % iFrame
end 