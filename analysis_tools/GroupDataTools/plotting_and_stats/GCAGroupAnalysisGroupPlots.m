function [ gcaGroupData ] = GCAGroupAnalysisGroupPlots(varargin)
%% GCAGroupAnalysisGroupPlots
% Function Designed for Data Exploration: 
% Plot one or more
%
% INPUT: gcaGroupData 
% (also sometimes referred to as 'toPlot' variable internally in these functions)
% a structure with fields of features to plot such that
% gcaGroupData.dataMat.featName{iGroup} = dataMat where each row is an observation and each
% column is data corresponding to an individual movie. 
% The key to this format is these
% can be easily horizontally cat to make a large array of all groups.
% this can be easily fed to boxplot to make individual boxplots and group
% boxplots
% gcaGroupData.info gives information about the project Lists, the group
% names, the colors to use for the plots etc.
%% Input check

ip = inputParser;

ip.CaseSensitive = false;
ip.StructExpand = false; 
ip.addOptional('gcaGroupData',[]); 
% Input/Output
ip.addParameter('OutputDirectory',[],@(x) ischar(x) || isempty(x));

% Features to examine
ip.addParameter('Interactive',true,@(x) islogical(x)); % if true will ask the user to select features to 
% examine if false will load the cell array entered by the feature
% parameter belwo 
ip.addParameter('Features',[]); % cell of features to be plotted: default is empty which flags to plot all the 
% feature fields included in the gcaGroupData structure. 

% Analysis Type
ip.addParameter('plotType',[]); % Options, 'perCell',  'pooled' 'perCellDistrb'
ip.addParameter('perNeuriteStat','nanmean' ); % if plotType is perCell will calculate this statistic per movie. 


% Plot Characteristics
ip.addParameter('makePlots',true); % if false will just add the stats. 
ip.addParameter('yLimOff',false,@(x) islogical(x));
ip.addParameter('splitMovie',false);
ip.addParameter('order',[]); % order to plot data from the cluster object (input cell array of names) 
ip.addParameter('FontSize',20);
ip.addParameter('groupIDs',[]);
ip.addParameter('yLims',[]); 

% Statistical options: by default will run a non-parametric pairwise permutation t-test 
% comparing each group to the control. 
ip.addParameter('multCompare',false); % flag to run multiple comparisons statistics 

% Data Export
ip.addParameter('writeCSV',true); 
ip.parse(varargin{:});
%% Initiate

if  isempty(ip.Results.gcaGroupData)
    [filename,pathname]  =  uigetfile(pwd,'Please select a gcaGroupData.mat File');
    load([pathname filesep filename]);
    toPlot = gcaGroupData; % keep the current downstream variable name 'toPlot' for now (historical reasons)
else 
   toPlot= ip.Results.gcaGroupData; 
end

if isempty(ip.Results.OutputDirectory)
   saveDir = [pwd filesep 'GroupDataPlots' filesep ip.Results.plotType];  
else
   saveDir = [ip.Results.OutputDirectory filesep 'GroupDataPlots' filesep ip.Results.plotType]; 
end 

if ~isdir(saveDir)
    mkdir(saveDir)
end 

if ~isempty(ip.Results.Features);
    params = ip.Results.Features;
else
    params = fieldnames(toPlot);
    params = params(~strcmpi('info',params)) ;
    
    if ip.Results.Interactive
        paramSelect  = listSelectGUI(params,[],'move');
        params  = params(paramSelect);
    end
end

if isempty(ip.Results.plotType) 
    plotType{1} = 'perCell'; 
    plotType{2} = 'pooled'; 
    plotType{3} = 'perCellDistrb'; 
    
    plotSelect =  listSelectGUI(plotType,1,'move'); 
    plotType = plotType{plotSelect}; 
else 
    plotType = ip.Results.plotType; 
end 

perNeuriteStat = str2func(ip.Results.perNeuriteStat);
toPlot = helperAddYLabels(toPlot);

if isempty(ip.Results.groupIDs)
    groupIDs = 1:numel(toPlot.info.names);
else
    groupIDs = ip.Results.groupIDs;
    grouping = arrayfun(@(i) repmat(i,size(toPlot.info.projList{groupIDs(i)},1),1),1:length(groupIDs),'uniformoutput',0);
    grouping = vertcat(grouping{:});
end
%% START
for iParam = 1:length(params)
    
    % get the dataMat (all groups horizontally catenated)
    toPad= max(cellfun(@(x) length(x),toPlot.(params{iParam}).dataMat(groupIDs)));
    forDataMat = cellfun(@(x) [x ;nan(toPad-length(x(:,1)),length(x(1,:)))],toPlot.(params{iParam}).dataMat(groupIDs),'uniformoutput',0);
    
    dataMatLargeC = horzcat(forDataMat{:});

    switch ip.Results.splitMovie
        case true
            
            %% Make the boxplot of the individual cells
            % Start boxplots : if user desires plot begin and end of movie separately for each neurite
            %setAxis('on');
            %setupFigure;
            projListAll = vertcat(toPlot.info.projList{:}) ;
            names = projListAll(:,2);
            names = strrep(names,'_',' ');
            
            pre=  arrayfun(@(x) [names{x} ' PreTreat'],1:numel(names),'uniformoutput',0);
            post = arrayfun(@(x) [names{x} ' PostTreat'],1:numel(names),'uniformoutput',0);
            
            labels =   arrayfun(@(x,y) [x; y],pre,post,'uniformoutput',0);
            labels = vertcat(labels{:});
            labels = cell(labels) ;
            
            % create color shades
            for iGroup= 1: numel(toPlot.info.names)
                colorShades{iGroup}=  [toPlot.info.colorShades{iGroup}(end,:) ; toPlot.info.colorShades{iGroup}(7,:)];
            end
            colorShadesFinal = vertcat(colorShades{:});
            
            h1 =  boxplot(dataMatLargeC,'colorGroup',toPlot.info.groupingPoolBeginEndMovie,...
                'colors',colorShadesFinal,'notch','on','outlierSize',1,'symbol','+','Labels',labels,'labelorientation','inline');
            if (~isfield(toPlot.(params{iParam}),'ylim') || ip.Results.yLimOff == true) ;
            else
                axis([0.5 size(dataMatLargeC,2)+ 0.5 0 toPlot.(params{iParam}).ylim]);
            end
            if isfield(toPlot.(params{iParam}),'yLabel');
                ylabel(toPlot.(params{iParam}).yLabel);
            else
                ylabel(strrep(params{iParam},'_',' '));
            end
            
            valuesM = nanmedian(dataMatLargeC,1);
            
            % perform the individual stats
            nums = [1,3,5];
            [idxPerm(:,1),pValuesPerm(:,1)] = arrayfun(@(i) permTest(dataMatLargeC(:,i) , dataMatLargeC(:,i+1),'CmpFunc',@nanmedian,'nRep',10000),nums);
            % lot the line median
            arrayfun(@(i) line([i,i+1], [valuesM(i),valuesM(i)],'Color',colorShades{1}(1,:)),nums,'uniformoutput',0);
            
            nums2 = [7,9,11];
            [idxPerm(:,2),pValuesPerm(:,2)] = arrayfun(@(i) permTest(dataMatLargeC(:,i) , dataMatLargeC(:,i+1),'CmpFun',@nanmedian,'nRep',10000),nums2);
            % lot the line median
            arrayfun(@(i) line([i,i+1], [valuesM(i),valuesM(i)],'Color',colorShades{2}(1,:)),nums2,'uniformoutput',0);
            
            % h1 = boxplot(dataMatLargeC,'colorGroup',toPlot.info.grouping,'notch','on','outlierSize',1,'colors',colors,'Labels',names,'labelorientation','inline');
            set(h1(:),'Linewidth',2);
            Results.(params{iParam}).pValuesPerm = pValuesPerm;
            Results.(params{iParam}).idxPerm = idxPerm;
            save([saveDir filesep 'perCell_TimeSplitStats_Results.mat'],'Results');
            
            saveas(gcf,[saveDir filesep 'perCellBox_TimeSplit' params{iParam} '.fig']);
            saveas(gcf,[saveDir filesep 'perCellBox_TimeSplit' params{iParam} '.png']);
            saveas(gcf,[saveDir filesep 'perCellBox_TimeSplit' params{iParam} '.eps'],'psc2');
            close gcf
            %% Percent Change
%             setAxis('on')
            % get the median value of each group before and after
            values = nanmedian(dataMatLargeC,1);
            values = values';
            ID =  toPlot.info.groupingPoolBeginEndMovie;
            percentChange(:,1) = (values(ID==2) - values(ID==1))./values(ID == 1);
            percentChange(:,2) = (values(ID ==4) -values(ID==3))./values(ID ==3);
            
            h = notBoxPlot(percentChange);
            
            %set the colors
            arrayfun(@(i) set(h(i).data,'markerFaceColor', colorShades{i}(1,:)),1:2);
            arrayfun(@(i) set(h(i).data,'markerEdgeColor','w'),1:2);
            arrayfun(@(i) set(h(i).mu,'color',colorShades{i}(1,:)),1:2);
            arrayfun(@(i) set(h(i).semPtch,'faceColor',colorShades{i}(2,:)),1:2);
            arrayfun(@(i) set(h(i).sdPtch,'faceColor',colorShades{i}(2,:)),1:2);
            
            forLabel = strrep(params{iParam},'_',' ');
            
            groupNames = toPlot.info.names;
            
            ylabel({forLabel ; 'Percent Change In Median'});
            %ylabel(toPlot.(params{iParam}).yLabel);
            set(gca,'XTick',1:numel(groupNames));
            set(gca,'XTickLabel',groupNames,'FontSize',20);
            
            h1 = get(gcf,'CurrentAxes');
            yLim = h1.YLim; % use whatever they used
            axis([0.5 2.5 yLim(1) yLim(2)]);
            
            saveas(gcf,[saveDir filesep 'percentChange_TimeSplit' params{iParam} '.fig']);
            saveas(gcf,[saveDir filesep 'percentChange_TimeSplit' params{iParam} '.png']);
            saveas(gcf,[saveDir filesep 'percentChange_TimeSplit' params{iParam} '.eps'],'psc2');
            
            close gcf
            %% Pooled boxplots split movie            
%             %% pooled boxplots
%             setAxis('off')
%            
%             
%             prePool=  cellfun(@(x) [x ' PreTreat'],toPlot.info.names,'uniformoutput',0);
%             postPool = cellfun(@(x) [x ' PostTreat'],toPlot.info.names,'uniformoutput',0);
%             labelsPool =   arrayfun(@(x,y) [x; y],prePool,postPool,'uniformoutput',0);
%             labelsPool = vertcat(labelsPool{:});
%             labelsPool = cell(labelsPool) ;
%             
%             
%             h = boxplot(dataMatLargeC,toPlot.info.groupingPoolBeginEndMovie,'notch','on',...
%                 'symbol', '+','outliersize',1,'colors',colorShadesFinal, ...
%                 'Labels',labelsPool,'labelorientation','inline');
%                        
%             set(h(:),'Linewidth',2);
%             nPooledGrps = length(unique(toPlot.info.groupingPoolBeginEndMovie));
%                         
%             if (~isfield(toPlot.(params{iParam}),'ylim') || ip.Results.yLimOff == true) ;
%             else
%                 axis([0.5 nPooledGrps + 0.5 0 toPlot.(params{iParam}).ylim]);
%             end
%             
%             if isfield(toPlot.(params{iParam}),'yLabel');
%                 ylabel(toPlot.(params{iParam}).yLabel);
%             else
%                 ylabel(strrep(params{iParam},'_',' '));
%             end
%             
%             set(gca,'XTick',1:nPooledGrps);
%             set(gca,'XTickLabel',labelsPool,'FontSize',10);
%             %     ylabel(params{iParam},'FontSize',30)
%             % get the median of group 1
%             %     values = toPlot.(params{iParam}){1}; % first group all projects
%             %     forLine = nanmedian(values(:));
%             %   if ~isempty(regexp(params{iParam},'Vel','Once'));
%             %   forLine =   forLine.*216./5;
%             % end
%             
%             %     line([0.5 numel(toPlot.info.names)+0.5] ,[forLine, forLine],'Linewidth',2,'color','k');
%             
%             % get the values for each group before and after treat
%             % forN will be a cell 2*the number of conditions long (for before and after the condition)
%             % each forN will contain a double matrix : max observation number x n neurite matrix
%             forN  = arrayfun(@(x) dataMatLargeC(:,toPlot.info.groupingPoolBeginEndMovie==x),1:nPooledGrps,'uniformoutput',0);
%             % combine the data for all neurites in the condition.
%             forN = cellfun(@(x) x(:),forN,'uniformoutput',0);
%             %forN = cellfun(@(x) x(:),toPlot.(params{iParam}).dataMat,'uniformoutput',0);
%             
%             medianValuesAll = cellfun(@(x) nanmedian(x),forN) ;
%             
%             N = cellfun(@(x) length(x(~isnan(x))),forN);
%             Nstring = num2cell(N);
%             Nstring = cellfun(@(x) ['N = '  num2str(x) ] ,Nstring,'uniformoutput',0);
%             title(Nstring);
%             
%             pairs = unique(toPlot.info.groupingPoolBeginEndMovie);
%             % row is before after treat, column is the condition
%             pairs = reshape(pairs,2,length(pairs)/2);
%             
%             for ipair = 1:size(pairs,2)
%                 [hit(ipair),pValues(ipair)] = permTest(forN{pairs(1,ipair)},forN{pairs(2,ipair)},'CmpFunction',@nanmedian);
%                 
%                 % plot line
%                 forLine = nanmedian(forN{pairs(1,ipair)});
%                 
%                 line([pairs(1,ipair)-0.25, pairs(2,ipair)+0.25],[forLine,forLine],...
%                     'color',colorShadesFinal(pairs(1,ipair),:),'linewidth',2);
%                 
%                 if hit(ipair)==0
%                     text(pairs(1,ipair),double(nanmedian(forN{pairs(1,ipair)}+2)),'NS','FontSize',10);
%                     % text(pairs(2,ipair),nanmedian(forN{pairs(2,ipair)}),'NS','FontSize',10);
%                 else
%                     text(pairs(1,ipair),double(nanmedian(forN{pairs(1,ipair)}+2)),num2str(pValues(ipair),4),'FontSize',10);
%                     
%                 end
%                 % test the difference between groups
%                 percentChange(ipair) = ( nanmedian(forN{pairs(2,ipair)}) -nanmedian(forN{pairs(1,ipair)})) / nanmedian(forN{pairs(1,ipair)});
%                 
%             end
%             
%             toPlot.(params{iParam}).percentChange= percentChange;
%             toPlot.(params{iParam}).pValues = pValues; % add a percent change field
%             clear percentChange
%             
%             saveas(gcf,[saveDir filesep 'pooled_TimeSplit' params{iParam} '.fig']);
%             saveas(gcf,[saveDir filesep 'pooled_TimeSplit' params{iParam} '.png']);
%             saveas(gcf,[saveDir filesep 'pooled_TimeSplit' params{iParam} '.eps'],'psc2');
            
            %% START if splitMovie false
        case false
            if isfield(toPlot.(params{iParam}),'yGroup');
                cDir = [saveDir filesep toPlot.(params{iParam}).yGroup];
                if ~isdir(cDir)
                    mkdir(cDir);
                end
                
            else
                cDir = saveDir;
            end % if isfield
            switch plotType
                case 'pooled'
                    grpNames = vertcat(toPlot.info.names(groupIDs));
                    
                    idxC = arrayfun(@(x) toPlot.info.groupingPoolWholeMovie(toPlot.info.groupingPoolWholeMovie==groupIDs(x)),1:length(groupIDs),...
                        'uniformoutput',false);
                    groupingPoolWholeMovie = vertcat(idxC{:});
                    forN  = arrayfun(@(x) dataMatLargeC(:,toPlot.info.groupingPoolWholeMovie==x),1:numel(grpNames),'uniformoutput',0);
                    forN = cellfun(@(x) x(:),forN,'uniformoutput',0);
                       if ip.Results.makePlots
                           % regular boxplot
                            close all 
                            
                            colors  = vertcat(toPlot.info.color{groupIDs});
                            grpNames = vertcat(toPlot.info.names(groupIDs));
                           
                            
                            h1 =  boxplot(dataMatLargeC,groupingPoolWholeMovie,...
                                'colors',colors,'notch','on','outlierSize',1,'symbol','+','Labels',grpNames,'labelorientation','inline');
                            
                            set(h1(:),'Linewidth',2);
 
                            cmpFunc = @nanmedian;
                            % for N median
                            percentChange = arrayfun(@(x) (cmpFunc(forN{x})-cmpFunc(forN{1}))./cmpFunc(forN{1}),2:numel(forN));
                            values = arrayfun(@(x) cmpFunc(forN{x}),1:numel(forN));
                            toPlot.(params{iParam}).percentChange= percentChange;
                            
                            
                            N = cellfun(@(x) length(x(~isnan(x))),forN);
                            Nstring = num2cell(N);
                            Nstring = cellfun(@(x) ['N = '  num2str(x) ] ,Nstring,'uniformoutput',0);
                            title(Nstring);
                            line([0.5,numel(grpNames)+0.5],[nanmedian(forN{1}),nanmedian(forN{1})],'color','k','linewidth',2);
                            
                            % Perform some simple stats against control
                            for i = 2:numel(grpNames)
                                [hit(i),pValue(i)] =   permTest(forN{1},forN{i},'CmpFunction',@median);
                                if hit(i) == 0
                                    text(i,nanmedian(forN{i}),'NS','color','k','FontSize',10);
                                else
                                    if pValue(i) == 0
                                        text( i,double(nanmedian(forN{i})),num2str(pValue(i),'%01d'),'color','k','FontSize',10);
                                    else
                                        text(i,double(nanmedian(forN{i})),num2str(pValue(i),'%04d'),'color','k','FontSize',10);
                                    end
                                end % if hit
                            end % for iGroup
                            
                            
                            ylabel(strrep(params{iParam},'_',' '));
 
                            saveas(gcf,[cDir filesep 'pooled_' params{iParam} '.fig']);
                            saveas(gcf,[cDir filesep 'pooled_' params{iParam} '.png']);

                        end % if makePlots 
                        
                        if ip.Results.writeCSV
                            % first take out the NaNs 
                            forN_noNans = cellfun(@(x) x(~isnan(x)),forN,'uniformoutput',0); 
                            dataMatPooled = reformatDataCell(forN_noNans);
                            
                            t = mat2dataset(dataMatPooled,'VarNames',vertcat(toPlot.info.names(:)));
                            
                            export(t,'file',[cDir filesep 'pooled_' params{iParam} '.csv']);
                            
                        end % write CSV
                      
                    %% plot per movie (cell)
                case 'perCell'
                    % create the group folder (Veil/Branch/Filo)
                    
                    if ip.Results.makePlots
                       % setAxis('on',0.95,20);
                       %setupFigure
                    end
                    if length(groupIDs) == numel(toPlot.info.names);
                        nCells = length(unique(toPlot.info.groupingPerCell));
                        % collect the data forN (note 2016030 not sure why I
                        % formated this way-  it it seems odd
                        % I compiled the groups only to again separate them...
                        % might just be due to history... check the split movie
                        % format)
                        
                        % put into a 1xc cell where c is the total number of
                        % neurite movies sampled over all groups
                        % each cell contains N features compiled per movie
                        
                        forN  = arrayfun(@(x) dataMatLargeC(:,toPlot.info.groupingPerCell==x),1:nCells,'uniformoutput',0);
                        forN = cellfun(@(x) x(:),forN,'uniformoutput',0);
                        
                    else
                        nCells = size(vertcat(toPlot.info.projList{groupIDs}),1);
                        forN = arrayfun(@(x) dataMatLargeC(:,x),1:nCells,'uniformoutput',0);
                        
                    end
                    % get the perNeurite stat of each neurite for the entire dataMat (condition info
                    % lost unfortunately) 1xc double where c the total number
                    % of neurite movies : the single row is distribution
                    % stat per neurite movie
                    
                    neuriteStatEachMovieCompile = cellfun(@(x) perNeuriteStat(x), forN);
                    
                    % add back the condition grouping info ( again I know a
                    % bit circuitous) 1xc cell where c is the number of
                    % perturbation groups, each cell is an rx1 double
                    % where r is the number of movies sampled per
                    % perturbation- and each row is a distribution stat as
                    % calculated for the entire neurite movie
                    if length(groupIDs) == numel(toPlot.info.names)
                        neuriteStatEachMovieGrouped =  arrayfun(@(x) neuriteStatEachMovieCompile(toPlot.info.grouping==x)',1:numel(toPlot.info.names),'uniformoutput',0);
                    else
                        neuriteStatEachMovieGrouped =  arrayfun(@(x) neuriteStatEachMovieCompile(grouping ==x)',1:length(groupIDs),'uniformoutput',0);
                    end
                    % reformat cell to a padded double for notBoxplot
                    % rxc matrix where r is the single measurement per movie (defined by the perNeuriteStat)
                    % and c is the perturbation group
                    perCellDataMat = reformatDataCell(neuriteStatEachMovieGrouped);
 
                    % sort the order of the perturbation groups if required
                    % (note this is helpful so that the order of groups reflects the
                    % clustergram output order for easy comparison)
                    if ~isempty(ip.Results.order)
                        namesC = ip.Results.order;
                        IDSort = cellfun(@(x) find(strcmpi(x,toPlot.info.names)),namesC);
                        IDSort = [1, IDSort];
                        perCellDataMat = perCellDataMat(:,IDSort);
                        names = ['KD Control' ; namesC'];
                        %names = ['Control' ; namesC'];
                        colors = toPlot.info.color(IDSort);
                        colorShades = toPlot.info.colorShades(IDSort);
                        neuriteStatEachMovieGrouped = neuriteStatEachMovieGrouped(IDSort);
                    else
                        colors = toPlot.info.color(groupIDs);
                        names = toPlot.info.names(groupIDs);
                        colorShades = toPlot.info.colorShades(groupIDs);
                    end
        %% Mult compare stats 
                    % perform some quick stat tests
                    if ip.Results.multCompare
                        [p,~,stats] = anova1(perCellDataMat,names,'off');
%                         setupFigure;
                        % make a folder
                        [c,~] =  multcompare(stats,'display','off');
                        
                        [c2,~] = multcompare(stats,'ctype','bonferroni','display','off'); 
                        
                        multSaveDir = [ cDir  filesep 'MultipleComparisons'];

                        if ~isdir(multSaveDir)
                            mkdir(multSaveDir);
                        end
                        title([]);
                        
                        if isfield(toPlot.(params{iParam}),'yLabel');
                            ylabel(toPlot.(params{iParam}).yLabel);
                        end
                        
                        save([multSaveDir filesep 'ANOVA1Stats_' params{iParam}],'stats','p');
                        
                        for i = 2:length(groupIDs)
                            [hit(i-1),pValues(i-1,1)] = permTest(neuriteStatEachMovieGrouped{1},neuriteStatEachMovieGrouped{i});
                        end
                        
                        fdr = mafdr(pValues,'BHFDR',true); 
                        fullMat = [pValues,fdr,c(1:(length(groupIDs)-1),end),c2(1:(length(groupIDs)-1),end)]; 
                        testMult = mat2dataset(fullMat); 
                        testMult.Properties.VarNames = {'ttestPerm','FDR','TukeyKramer','Bonferroni'} ; 
                        
                        testMult.Properties.ObsNames = names(2:end)' ; 
                        
                        export(testMult,'file',[multSaveDir filesep 'fullMatMultCompare_' params{iParam} '.csv']); 
                        save([multSaveDir filesep 'fullMatMultCompare_' params{iParam}],'fullMat'); 
                        clear pValues
%                         save([multSaveDir filesep 'ANOVA1Stats_Tukey-Kramer_' params{iParam}],'c'); 
%                         save([multSaveDir filesep 'ANOVA1Stats_Bon_' params{iParam}],'c2'); 
                      
%                         saveas(gcf,[multSaveDir filesep 'MultCompare_' params{iParam} '.fig']);
%                         saveas(gcf,[multSaveDir filesep 'MultCompare_' params{iParam} '.eps'],'psc2');
%                         saveas(gcf,[multSaveDir filesep 'MultCompare_' params{iParam} '.png']); 
%                         close gcf
                        pValues = c2([1:7],end); 
                        hit = pValues<0.1; 
                        pValues = [1;pValues]; 
                        hit = [0;hit]; 
                    else % run the pair-wise permutation t-test of the means 
                        for i = 2:length(groupIDs)
                            [hit(i),pValues(i)] = permTest(neuriteStatEachMovieGrouped{1},neuriteStatEachMovieGrouped{i});
                        end
                    end % mult compare
                    if ip.Results.makePlots
%                         setupFigure
                        h = notBoxPlot(perCellDataMat);
                        
                        %set the colors
                        arrayfun(@(i) set(h(i).data,'markerFaceColor', colors{i}),1:numel(names));
                        arrayfun(@(i) set(h(i).data,'markerEdgeColor','w'),1:numel(names));
                        arrayfun(@(i) set(h(i).mu,'color',colors{i}),1:numel(names));
                        arrayfun(@(i) set(h(i).semPtch,'faceColor',colorShades{i}(4,:)),1:numel(names));
                        arrayfun(@(i) set(h(i).sdPtch,'faceColor',colorShades{i}(1,:)),1:numel(names));
                        arrayfun(@(i) set(h(i).data,'markerSize',5),1:numel(names));
 
                    for i = 2:length(groupIDs)
                       
                            if hit(i) == 0
                                text(i,double(mean(neuriteStatEachMovieGrouped{i})),'NS','FontSize',10);
                            else
                                if pValues(i) == 0
                                    text(i,double(mean(neuriteStatEachMovieGrouped{i})),num2str(pValues(i),1),'FontSize',10);
                                else
                                    text(i,double(mean(neuriteStatEachMovieGrouped{i})),num2str(pValues(i),4),'FontSize',10);
                                end
                            end
                        
                    end % length groups
       
                        forLabel = strrep(params{iParam},'_',' ');

                        ylabel(forLabel);
                        %ylabel(toPlot.(params{iParam}).yLabel);
                        set(gca,'XTick',1:numel(names));
                        %set(gca,'XTickLabel',names,'FontSize',ip.Results.FontSize);
                        set(gca,'XTickLabel',names); 
                        
%                         axis([0.5 length(toPlot.info.names)+ 0.5 0 toPlot.(params{iParam}).ylim]);
                        
                        if isfield(toPlot.(params{iParam}),'yLabel');
                            ylabel(toPlot.(params{iParam}).yLabel);
                        end
                        set(gca,'XTickLabelRotation',45);
                        if isempty(ip.Results.yLims)
                        h1 = get(gcf,'CurrentAxes');
                        yLim = h1.YLim; % use whatever they used
                        else 
                            yLim = ip.Results.yLims; 
                        end 
                        axis([0.5 length(groupIDs)+ 0.5  yLim(1) yLim(2)]);
                        
                        saveas(gcf,[cDir filesep 'perCell' params{iParam} '.fig']);
                        saveas(gcf,[cDir filesep 'perCell' params{iParam} '.png']);
                        %saveas(gcf,[cDir filesep 'perCell' params{iParam} '.eps'],'psc2');
                        saveas(gcf,[cDir filesep 'perCell' params{iParam} '.psc2']); 
                        close gcf
                    end % makePlots
                    
                    % if saveCSV 
                    if ip.Results.writeCSV
                      t =  mat2dataset(perCellDataMat,'VarNames',vertcat(toPlot.info.names(:))); 
                      export(t,'file',[cDir filesep 'perCell' params{iParam} '.csv']); 
                    end 
                   
                    if exist('pValues')                        
                        toPlot.(params{iParam}).pValues.perCell = pValues;
                    end
                    %%
                case 'perCellDistrb'
                    
                    close gcf
                    setAxis('on',0.95,20);
                    movieIDs = vertcat(toPlot.info.projList{groupIDs});
                    movieIDs = movieIDs(:,2);
                    movieIDs = cellfun(@(x) strrep(x,'_',' '),movieIDs,'uniformoutput',0);
                    forLabel = strrep(params{iParam},'_',' ');
                    
                    ylabel(forLabel);
                    idxC = arrayfun(@(x) toPlot.info.grouping(toPlot.info.grouping==groupIDs(x)),1:length(groupIDs),'uniformoutput',0); 
                    grouping = vertcat(idxC{:}); 
                    colors = vertcat(toPlot.info.color{groupIDs});
                    boxplot(dataMatLargeC,'ColorGroup',grouping, ...
                        'colors',colors,'notch','on','outlierSize',1,'symbol','+');
                    
                    set(gca,'XTick',1:numel(movieIDs));
                    set(gca,'XTickLabel',movieIDs,'FontSize',12);
                    set(gca,'XTickLabelRotation',45);
                    
                    % get the axis limits
                    data = dataMatLargeC(:);
                    data = data(~isnan(data));
                    maxY = prctile(data,99.5);
                    minY = prctile(data,0.05);
                    axis([0.5 size(dataMatLargeC,2) + 0.5, minY, maxY]);
                    
                  
                    
                    saveas(gcf,[cDir filesep 'perCellDistr' params{iParam} '.fig']);
                    saveas(gcf,[cDir filesep 'perCellDistr' params{iParam} '.png']);
                    saveas(gcf,[cDir filesep 'perCellDistr' params{iParam} '.eps'],'psc2');
                    close gcf
            end  % switch plotType
    end % switch splitMovie
    clear pValues
end % iParam (ie Feat)
%%
gcaGroupData = toPlot; 
save([saveDir filesep 'gcaGroupDataWithStats.mat'],'gcaGroupData');
end % function