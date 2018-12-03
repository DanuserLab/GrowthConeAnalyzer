function [gcaGroupData] = GCAGroupAnalysisCollectFeatureDataForStats(varargin)
%GCAGroupAnalysisCollectFeatureDataForStats :
% INPUT: gcaGroupData structure (also internally referred to as toPlot
% for historic reasons... with field .info
%                                         .names : name of group
%                                         .projList : projList associated
%                                          with group with individual movie IDs
%                                         .grouping : vector of group
%                                          numbers
% OUTPUT: gcaGroupData structure with all features that have been run added
%% Input check
ip = inputParser;

ip.CaseSensitive = false;

ip.StructExpand = false;

ip.addOptional('gcaGroupData',[]); % optional : can enter the loaded structure 
% directly 
ip.addParameter('InputPath',[]); % or load it here with the full path name 

ip.addParameter('Interactive',true); 
ip.addParameter('OutputDirectory',[]);

ip.addParameter('clearOldFields',false);
ip.addParameter('splitMovie',false);
ip.addParameter('splitFrame', 62);  % last frame you want to include

ip.addParameter('FeatureFolder','GCAFeatureExtraction');
ip.addParameter('perFrame',false); % will collect the median value per frame

ip.parse(varargin{:});
%% Initiate

if  isempty(ip.Results.gcaGroupData)
    if isempty(ip.Results.InputPath)
        [filename,pathname]  =  uigetfile(pwd,'Please select a gcaGroupData.mat File');
        load([pathname filesep filename]);
    else
        load(ip.Results.InputPath);
    end    
    toPlot = gcaGroupData; % keep the current downstream variable name 'toPlot' for now (historical reasons)
else
    toPlot= ip.Results.gcaGroupData;
end

nGroups = numel(toPlot.info.names);
if ip.Results.clearOldFields
    params = fieldnames(toPlot);
    params = params(~strcmpi(params,'info')); % keep the info
    toPlot = rmfield(toPlot,params);
end

if isempty(ip.Results.OutputDirectory)
    if ip.Results.Interactive
        outDir= uigetdir(pwd,'Please select a directory to store your gcaGroupDataWithFeatureValues.mat file');
    else
        outDir = pwd ;
    end
else
    outDir = ip.Results.OutputDirectory;
end

if ~isdir(outDir)
    mkdir(outDir);
end

%% START
for iGroup = 1:nGroups
    
    projListC = toPlot.info.projList{iGroup}(:,1);
    
    nMovies = size(projListC,1);
    
    if ip.Results.splitMovie
        
        % Grouping Var1 : grouping per condition
        % create the grouping variable for pooling full group data
        % [1,(Repeated 2*nCellsProj1 times), 2(Repeated
        % 2*nCellsProj2)....[n,(Repeated 2*ncellsProjN times)]
        grpVar{iGroup} = repmat(iGroup,nMovies*2,1); %
        
        % Grouping Var2  :   grouping per cell
        % [1,1,2,2,3,3,...numCellsTotal,numCellsTotal]
        g = arrayfun(@(x) repmat(x,2,1),1:nMovies,'uniformoutput',0);
        grpVar2{iGroup} = size(vertcat(toPlot.info.projList{1:iGroup-1}),1) + vertcat(g{:});
        
        % Grouping Var3 : grouping per treatment
        g3 = repmat([1,2],1,nMovies)';
        grpVar3{iGroup} = g3 + 2*(iGroup-1);
        
    end % ip.Results.splitMovie
    
    % check for consistency among the parameters that were run.
    for iProj = 1:nMovies
        
        %load([projListC{iProj}  filesep 'analysis' filesep 'movieData.mat']);
        load([projListC{iProj}  filesep 'GrowthConeAnalyzer' filesep 'movieData.mat']);
        
        parameterDir = [MD.outputDirectory_ filesep ip.Results.FeatureFolder];
        
        % might also include ylabel name and ylim for each parameter and
        % read in each Descriptor directory to keep constant.
        
        % search all descriptor parameters.
        %localParamFiles = searchFiles('meas_',[],[parameterDir filesep 'Descriptor'],1);
        localParamFiles = searchFiles('meas_','csv',parameterDir,1);
        paramNamesC = cellfun(@(x) strrep(x,'meas_',''),localParamFiles(:,1),'uniformoutput',0);
        paramNamesC = cellfun(@(x) strrep(x,'.mat',''),paramNamesC,'uniformoutput',0);
        
        %
        for iParam = 1:numel(paramNamesC)
            
            % collect the data for each parameter.
            load([localParamFiles{iParam,2} filesep localParamFiles{iParam,1}]);
            
            if ip.Results.splitMovie
                
                % currently assumes only splitting movie in two pool these
                % values
                dataSetGroup.(paramNamesC{iParam}).valuesWholeMovie{2*(iProj-1)+1} = vertcat(measC{1:ip.Results.splitFrame});
                dataSetGroup.(paramNamesC{iParam}).valuesWholeMovie{2*(iProj-1)+2} = vertcat(measC{ip.Results.splitFrame+1:end});
                
            else
                dataSetGroup.(paramNamesC{iParam}).valuesWholeMovie{iProj} = vertcat(measC{:});
            end
            %% add the per frame values
            if ip.Results.perFrame
                if ip.Results.splitMovie
                    dataSetGroup.(paramNamesC{iParam}).valuesPerFrame{2*(iProj-1)+1} = vertcat(measC{1:ip.Results.splitFrame});
                    dataSetGroup.(paramNamesC{iParam}).valuesPerFrame{2*(iProj-1)+2} = vertcat(measC{ip.Results.splitFrame+1:end});
                else
                    
                    idx = find(cellfun(@(x) isempty(x),measC));
                    
                    for i = 1:length(idx)
                        measC{idx(i)} = NaN;
                    end
                    
                    frameEnd = numel(measC);
                    valuesCell = cellfun(@(x) nanmedian(x), measC,'uniformoutput',0); % take the median per frame
                    valuesCell = vertcat(valuesCell{:});
                    % for now always truncate from 1:119
                    valuesCell = valuesCell(1:frameEnd);
                    dataSetGroup.(paramNamesC{iParam}).valuesPerFrame{iProj} = valuesCell;
                end   % if splitMovie
            end     % if ip.Results.perFrame
        end % for iParam
    end % for iProj
    %%  reformat all features (params)
    paramsAll  = fieldnames(dataSetGroup);
    %
    reformat = arrayfun(@(i) reformatDataCell(dataSetGroup.(paramsAll{i}).valuesWholeMovie),1:numel(paramsAll),...
        'uniformoutput',0);
    if ip.Results.perFrame
        reformat2 = arrayfun(@(i) reformatDataCell(dataSetGroup.(paramsAll{i}).valuesPerFrame),1:numel(paramsAll),...
            'uniformoutput',0);
    end % perFrame
    
    for iParam = 1:numel(paramsAll)
        toPlot.(paramsAll{iParam}).dataMat{iGroup} = reformat{iParam};
        %         [~,hashTag] = gcaArchiveGetGitHashTag;
        %         toPlot.(paramsAll{iParam}).hashTag = hashTag;
        toPlot.(paramsAll{iParam}).timestamp = clock;
        toPlot.(paramsAll{iParam}).featFolder = ip.Results.FeatureFolder;
        if  ip.Results.perFrame
            %toPlot.(paramsAll{iParam}).dataMatPerFrame{iGroup} = reformat2{iParam};
            % QandD 20180422
            toPlot.(paramsAll{iParam}).dataMatPerFrame{iGroup} = reformat2{iParam}(1:119,:);
            toPlot.(paramsAll{iParam}).timestampPerFrame = clock;
            toPlot.(paramsAll{iParam}).featFolderPerFrame = ip.Results.FeatureFolder;
        end % if ip.Results.perFrame
        
    end % for iParam (param here  means feature)
    clear reformat paramsAll dataSetGroup
end  % for iGroup
%% Add the grouping variable information
if ~ip.Results.splitMovie
    toPlot.info.groupingPerCell = 1:size(vertcat(toPlot.info.projList{:}),1);
    toPlot.info.groupingPoolWholeMovie = toPlot.info.grouping;
else
    toPlot.info.groupingPoolWholeMovie = vertcat(grpVar{:});
    toPlot.info.groupingPerCell= vertcat(grpVar2{:});
    toPlot.info.groupingPoolBeginEndMovie = vertcat(grpVar3{:});
end
gcaGroupData = toPlot;
if ip.Results.splitMovie
    save([outDir filesep 'gcaGroupDataWithFeatureValues_SplitMoviesAtFrame' num2str(ip.Results.splitFrame) '.mat'],'gcaGroupData');
else
    
    save([outDir filesep 'gcaGroupDataWithFeatureValues.mat'],'gcaGroupData');
end
end % function