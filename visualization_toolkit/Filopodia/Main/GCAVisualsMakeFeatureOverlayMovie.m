function GCAVisualsMakeFeatureOverlayMovie(varargin)
%% GCAVisualsMakeFeatureOverlayMovie
%
% This function is a wrapper to visualize features.
% Default is to search for all features already run for the movie of interest 
% and then ask the user interactively which type of visualization they
% would like to view.
%
%% 
% INPUT:
%
%  MD movieData object (Optional): 
%  If no input will ask the user to interactively load 
%  else please load the MD file and input it into the function directly. 
%  GCAVisualsMakeFeatureOverlayMovie(MD). 
%
% 
% %%INPUT/OUTPUT%%
% 'FeatureDirectory': (PARAM) : character array
%                              Full path to the current feature extraction directory
%                              (DEFAULT: [MD.outputDirectory_ filesep 'GCAFeatureExtraction'];
%                                Note the function by default will search for all features run in this folder
%                                by looking for meas_*.mat files. Where * indicates the name of the
%                                feature.  (Note these files do not have to be directly in the feature
%                                folder, and can be nested)). 
% 
% 'OutputDirectory': (PARAM) : character array 
%                     Full path where you would like your overlay images stored 
%   (DEFAULT: save the movie in the given feature folder of the feature 
%   directory) 
%
% %% FEATURE SELECTION %%
% 'interactive' (PARAM) :  logical                  
%                          Flag to ask the user to input features for
%                          overlays interactively. If true a listSelectGUI 
%                          will pop up with all possible overlay options. 
%  (DEFAULT: true) 
% 
% %% MOVIE DETAILS %% 
% 'frames':  (PARAM) : numeric vector
%                     vector of the frame numbers for which the user would
%                     like to make overlays.
% (DEFAULT: make full movie of overlays). 

%% Check input
ip = inputParser;

ip.CaseSensitive = false;

ip.StructExpand = false;

ip.addOptional('MD',[]); % if empty will ask to load interactively

%%% INPUT/OUTPUT %%%
ip.addParameter('ChannelIndexOverlay',1); % Channel over which to overlay
ip.addParameter('ChannelIndexDetection',1); % Channel of the filopodia
ip.addParameter('ChannelIndexVeil',1);

ip.addParameter('FeatureDirectory',[],@(x) ischar(x) || isempty(x));
ip.addParameter('InputDirectory', [],@(x) ischar(x) || isempty(x)); % directory where 
% the original filoBranch.mat used for the feature extraction is stored. 
ip.addParameter('VeilDirectory',[],@(x) ischar(x) || isempty(x)); % directory where 
% the veil/stem information is stored (necessary information for the
% density overlays)
ip.addParameter('OutputDirectory',[],@(x) ischar(x) || isempty(x));

ip.addParameter('NonGCImage',false); 

%%% FEATURE SELECTION %%%
ip.addParameter('interactive',true,@(x) islogical(x)); % allows user to
% click on a feature that has been extracted and make the movie

ip.addParameter('features',[],@(x) iscell(x)); % cell array of feature names
% if interactive is turned to false these features will be selected.
% If empty will make feature overlays for all features extracted.

%%% MOVIE DETAILS %%%
ip.addParameter('visible','off'); % if on show figure while plotting

ip.addParameter('frames', []);

ip.addParameter('ScaleBar',false,@(x) islogical(x) );
ip.addParameter('ScaleBarSize',10); % in um
ip.addParameter('Timer',false,@(x) islogical(x));
ip.addParameter('FontSizeTimer',18);

ip.addParameter('TreatmentFrame',[]);
ip.addParameter('TreatmentTitle','CK666 Added');

ip.addParameter('plotText',false); % flag to plot individual feature values
ip.addParameter('plotIDs',false) % flag to ignore feature values and color/write
% text according to the filopodia segment ID number in the data structure.

ip.addParameter('cMapLimits','defaults'); % color map limits : if defaults will use below

% defaults for color map limits
defaults{1,1} = 'filoDensityAlongVeil'; defaults{1,2} = [0,10];
defaults{2,1} = 'filoOrientation'; defaults{2,2} = [0,180];
defaults{3,1} = 'filoIntensityEmbedded';defaults{3,2} = [0.5,2];
defaults{4,1} = 'filoIntensityToVeil'; defaults{4,2} = [0.5,2];
defaults{5,1} = 'filoLengthEmbedded';defaults{5,2} = [0,10];
defaults{6,1} = 'filoLengthFullActinBundle';defaults{6,2} = [0,10];
defaults{7,1} = 'filoLengthToVeil'; defaults{7,2} = [0,10];
defaults{8,1} = 'filoMaxCurvature'; defaults{8,2} = [0,1.5]; 
defaults{9,1} = 'branchLength_2ndOrder'; defaults{9,2} =[0,10];
defaults{10,1} = 'branchOrientation_2ndOrder' ; defaults{10,2} = [0,180];
defaults{11,1} = 'branchIntensity_2ndOrder' ; defaults{11,2} = [0 2];
defaults{12,1} = 'validation' ; defaults{12,2} = [0,10];
defaults{13,1} = 'percentEachActinBundleEmbed' ; defaults{13,2} = [0,1];
defaults{14,1} = 'branchMaxCurvature_2ndOrder'; defaults{14,2} = [0,0.5];
defaults{15,1} = 'ForMainMovie'; defaults{15,2} = [0,10];
defaults{16,1} = 'branchDensity_2ndOrder' ; defaults{16,2} = [0,10];
defaults{17,1} = 'filoIntensityEmbedded_Norm'; defaults{17,2} = [0.5,2];
defaults{18,1} = 'filoIntensityToVeil_Norm'; defaults{18,2} = [0.5,2];

ip.addParameter('minMaxDefaults',defaults);

ip.addParameter('colorbarOverlay',true);
ip.addParameter('colorByValue',true); % color individual filopodia by value
% if applicable.
ip.addParameter('extraColor',[]); % adds an extra color to the color bar so that NaN values are white
ip.addParameter('colorFiloBranch',[0,0,0]); % default color is black
ip.addParameter('colorVeilStem',[0,0,0]); % default color is black

ip.addParameter('otherFiles',true);
% ip.addParameter('SubRegionFlag',false,@(x) islogical(x))

ip.addParameter('UseSmoothedCoords',false);

ip.addParameter('firstFrameLimits',false);

ip.addParameter('plotSigma',[]);

ip.addParameter('filoIDs',[]); % if empty will just use the filter

ip.addParameter('screen2png',true);

ip.addParameter('plotFilteredExt',false);

%%% FLAG FOR MASK %%%
ip.addParameter('createMask',false);
ip.addParameter('OutputDirectoryMask',[]);

ip.addParameter('OutputDirectoryMaskBack',[]);
ip.addParameter('createBackMask',false);
ip.addParameter('BMDilRad',10);

ip.parse(varargin{:});
%% Initiate
if isempty(ip.Results.MD);
    [MDname, path] =  uigetfile(pwd,'Please Select a movieData.mat File' );
    load([path filesep MDname]);
else
    MD = ip.Results.MD;
end

if isempty(ip.Results.FeatureDirectory)
    measDir = [MD.outputDirectory_ filesep 'GCAFeatureExtraction' ];
else
    measDir = ip.Results.FeatureDirectory;
end

if isempty(ip.Results.InputDirectory)
    
    inDir = [MD.outputDirectory_ filesep 'SegmentationPackage' ...
        filesep 'StepsToReconstruct' filesep 'VII_filopodiaBranch_fits' filesep 'Channel_' num2str(ip.Results.ChannelIndexDetection)];
else
    inDir = ip.Results.InputDirectory;
end % isempty

if isempty(ip.Results.VeilDirectory)
    if ~ip.Results.NonGCImage
        veilDir = [MD.outputDirectory_ filesep 'SegmentationPackage' ...
            filesep 'StepsToReconstruct' filesep 'IV_veilStem_length' filesep 'Channel_' num2str(ip.Results.ChannelIndexVeil)];
    else
        veilDir = [MD.outputDirectory_ filesep 'SegmentationPackage' ...
            filesep 'StepsToReconstruct' filesep 'III_veilStem_reconstruction' filesep 'Channel_' num2str(ip.Results.ChannelIndexVeil)];
    end
else
    veilDir = ip.Results.VeilDirectory;
end % isempty

if ip.Results.createMask
    if isempty(ip.Results.OutputDirectoryMask)
        outDirMask = [MD.outputDirectory_ filesep 'SegmentationPackage' filesep 'GCAMasks' filesep 'Channel_' num2str(ip.Results.ChannelIndexDetection)];
    else
        outDirMask = ip.Results.OutputDirectoryMask;
    end
    
    if ~isdir(outDirMask)
        mkdir(outDirMask);
    end
    
    if ip.Results.createBackMask
        if isempty(ip.Results.OutputDirectoryMaskBack)
            outDirMaskBack = [MD.outputDirectory_ filesep 'GCABackgroundMasks' filesep 'Channel_' num2str(ip.Results.ChannelIndexDetection)];
        else
            outDirMaskBack = ip.Results.OutputDirectoryMaskBack;
        end
        if ~isdir(outDirMaskBack)
            mkdir(outDirMaskBack);
        end
    end % if createBackMask
end % isempty
%% go into each folder and look for measC feature .mat

% might also include ylabel name and ylim for each parameter and
% read in each Descriptor directory to keep constant.
if exist([measDir filesep 'Descriptor' filesep 'Filopodia'],'dir')~=0
    searchDir = [measDir filesep 'Descriptor' filesep 'Filopodia'];
else
    searchDir = measDir;
end

% search all features in folder
localParamFiles = searchFiles('meas_','.csv',searchDir,1);

paramNamesC = cellfun(@(x) strrep(x,'meas_',''),localParamFiles(:,1),'uniformoutput',0);
paramNamesC = cellfun(@(x) strrep(x,'.mat',''),paramNamesC,'uniformoutput',0);

% if interactive
% ask the user for which feature they would like to make a sanity movie
% most measurments are going to be simple it will just be plotting
% the filopodia of the filter set with the value saved
if ip.Results.interactive
    paramSelect  = listSelectGUI(paramNamesC,[],'move');
    selected  = paramNamesC(paramSelect);
else % make all movies
    if ~isempty(ip.Results.features);
        selected = ip.Results.features;
        % find those features from the parameter names
        test = arrayfun(@(x) strcmpi(selected{x},paramNamesC),1:numel(selected),'uniformoutput',0);
        paramSelect = cellfun(@(x) find(x),test);
    else % loop through all found
        selected = paramNamesC ;
        paramSelect = 1:numel(selected);
    end
end % interactive
%% load the filoInfo and veil/stem
load([inDir filesep 'filoBranch.mat']);
imgSize = MD.imSize_;
load([veilDir filesep 'veilStem.mat']);

if isempty(ip.Results.frames)
    frames = 1:MD.nFrames_;
else
    frames = ip.Results.frames;
end

%% will check to make sure the filter and the number of frames is equivalent

nFrames = length(frames);

%% Loop over selected features
for iSelect = 1:numel(selected)
    
    if strcmpi(selected{iSelect},'filoDensityAlongVeil') || ...
            strcmpi(selected{iSelect}, 'filoBranchComplexity') || ...
            strcmpi(selected{iSelect},'branchDistanceFromVeil') || ...
            strcmpi(selected{iSelect},'percentTotalActinBundlesVeilEmbedded');
        colorByValue = false;
    else
        colorByValue = ip.Results.colorByValue;
    end

    if isempty(ip.Results.OutputDirectory)
        % make the directory to save the measurement movie
        outDirPre = [localParamFiles{paramSelect(iSelect),2} filesep selected{iSelect} '_Feature_Movie' filesep ...
            'Channel' num2str(ip.Results.ChannelIndexDetection) 'Detect_OverlaidOnChannel' num2str(ip.Results.ChannelIndexOverlay)];
    else
        outDirPre = ip.Results.OutputDirectory;
    end
    
    if ip.Results.plotIDs
        add1 = 'PlotIDs';
    else
        add1 = [];
    end
    
    if ip.Results.plotText
        add2 = 'PlotText';
    else
        add2=  [];
    end
    
    outDir  = [outDirPre '_' add1 '_' add2 ];
    
    if ~isdir(outDir)
        mkdir(outDir) ;
    end
    
    load([localParamFiles{paramSelect(iSelect),2} filesep localParamFiles{paramSelect(iSelect),1}]);
    
    x = upDirectory(localParamFiles{paramSelect(iSelect),2},1);
    
    load([x filesep 'filoFilterSet.mat']);
    
    % check that measC is the same size as nFrames
    if nFrames > numel(measC);
        nFrames = numel(measC);
        %         display(['The number of frames selected is'...
        %             'greater than the number of measurement frames: Truncating']);
        frames = 1:numel(measC);
    end
    
    % Loop overs frames in the movie
    for iFrame = 1:nFrames
        frameC = frames(iFrame);
        
        filoInfo = filoBranch(frameC).filoInfo;
        
        if ip.Results.plotIDs
            plotValues = 'IDs';
        else
            plotValues = measC{frameC};
        end
        
        if  iFrame ==1
            if ip.Results.firstFrameLimits
                lims = [min(-img(:)) max(-img(:))];
            else
                lims = [];
            end
        end
        if ~isempty(plotValues);
            if colorByValue && ~strcmpi(plotValues,'IDs')
                if strcmpi(ip.Results.cMapLimits,'defaults')
                    [cMapLimits] =  defaults{strcmpi(selected{iSelect},defaults(:,1)),2} ;
                elseif isempty(ip.Results.cMapLimits)
                    cMapLimits(1) = min(plotValues);
                    cMapLimits(2) = max(plotValues);
                else
                    cMapLimits = ip.Results.cMapLimits;
                end
            elseif strcmpi(plotValues,'IDs')
                cMapLimits = [1,length(filoInfo)];
            else
                cMapLimits(1) = min(plotValues);
                cMapLimits(2) = max(plotValues);
            end
            
            if isempty(ip.Results.filoIDs)
                filterSetC= filoFilterSet{frameC};
            else
                x =  zeros(length(filoInfo));
                x(ip.Results.filoIDs)=1;
                filterSetC = logical(x);
            end
            
            img = double(imread([MD.getChannelPaths{ip.Results.ChannelIndexOverlay} filesep MD.getImageFileNames{ip.Results.ChannelIndexOverlay}{frameC}]));
            
            
            setFigure(imgSize(2),imgSize(1),ip.Results.visible);
            
            imshow(-img,lims) ;
            hold on
            
            %% Check for branch feature input
            branchMode = (~isempty(regexpi(selected{iSelect},'branchDistanceFrom'))  || ~isempty(regexpi(selected{iSelect},'branchDensity')));
            
            plotTextAtBranches = (~isempty(regexpi(selected{iSelect},'branchDistanceFrom')));
            %% Call Curvature Plotting Function if Necessary
            
            if strcmpi(selected{iSelect},'filoMaxCurvature') || strcmpi(selected{iSelect}, 'branchMaxCurvature_2ndOrder')
                
                GCAVisualsColorCodeByCurvature(filoInfo,'filoFilterSet',filterSetC,'cMapLimits',cMapLimits,'pix2Micron',MD.pixelSize_/1000);
                
                %% plot percentActinBundlesVeilEmbedded
            elseif strcmpi(selected{iSelect},'percentTotalActinBundlesVeilEmbedded');
                
                plotText{1} = true;
                plotText{2} = false;
                filterSetForPlot{1}=  filterSetC(:,1);
                filterSetForPlot{2} = filterSetC(:,1) & filterSetC(:,2);
                
                for iFilter = 1:2
                    
                    GCAVisualsFilopodiaMeasurementOverlays(filoInfo,imgSize,...
                        'filoFilterSet',filterSetForPlot{iFilter},'plotValues',plotValues,...
                        'branchMode',branchMode,'colorByValue',false,'plotText',plotText{iFilter},'justExt',...
                        iFilter,'extraColor',ip.Results.extraColor,'cMapLimits',cMapLimits,'UseSmoothedCoords',ip.Results.UseSmoothedCoords);
                end               
                
            elseif (strcmpi(selected{iSelect},'filoLengthFullActinBundle') || strcmpi(selected{iSelect},'forMainMovie'))
                % plot each filopodia by the color of the full actin
                % bundle length : do not plot embedded bundles that do
                % not pass the criter
                
                plotText{1} = ip.Results.plotText;
                plotText{2} = false;
                %filterSetForPlot{1}=  filterSetC(:,1)   ;  % both have to be significant
                %filterSetForPlot{2} = filterSetC(:,1) & filterSetC(:,2);
                filterFrameC = filterSetC(:,1);
                
                filterInt = (filterSetC(:,1) == 1 & filterSetC(:,2) ==0 ); % get the ID of all non-fits internally this filter is the length of
                %             % the original filoInfo detection
                filterInt = filterInt(filterFrameC); % keep only the
                filoInfoExtBund = filoInfo(filterFrameC);
                filoInfoIntBund = filoInfoExtBund(~filterInt);

                filoInfoFilt{1} = filoInfoExtBund; % should be the same size as the numbers
                filoInfoFilt{2} = filoInfoIntBund;
                if strcmpi(plotValues,'IDs')
                    plotValues = 1:length(filoInfo);
                    plotValues = plotValues';
                    plotValuesSub{1} = plotValues(filterFrameC);
                    plotValuesSub{2} =[];
                    cMapLimits = [min(plotValuesSub{1}),max(plotValuesSub{1})];
                else
                    
                    plotValuesSub{1} = plotValues;
                    plotValuesSub{2} = plotValues(~filterInt);
                end
                %
                filoMaskCell = cell(2,1);
                for i = 1:2
                    
                    filoMaskCell{i} =    GCAVisualsFilopodiaMeasurementOverlays(filoInfoFilt{i},imgSize, ...
                        'plotValues',plotValuesSub{i},'colorByValue',true,'plotText',plotText{i},'justExt',i,...
                        'extraColor',ip.Results.extraColor,'cMapLimits',cMapLimits,'UseSmoothedCoords',ip.Results.UseSmoothedCoords,...
                        'plotSigma',ip.Results.plotSigma,'createMask',ip.Results.createMask);
                    
                end
                filoMask = filoMaskCell{1}; % for now just get the external filo for the mask
                if ip.Results.plotFilteredExt
                    GCAVisualsFilopodiaMeasurementOverlays(filoInfo(~filterFrameC),imgSize,...
                        'justExt',1,'plotFiltered',true);
                end
            elseif  strcmpi(selected{iSelect},'percentEachActinBundleEmbed')
                
                
                filoInfoFilt= filoInfo(filterSetC(:,1) & filterSetC(:,2));
                
                for i = 1:2
                    
                    GCAVisualsFilopodiaMeasurementOverlays(filoInfoFilt,imgSize, ...
                        'plotValues',plotValues,'colorByValue',true,'plotText',ip.Results.plotText,'justExt',i,...
                        'extraColor',ip.Results.extraColor,'cMapLimits',cMapLimits,'UseSmoothedCoords',ip.Results.UseSmoothedCoords,...
                        'plotSigma',ip.Results.plotSigma);
                end
                %% filoBranch complexity metric visualization
            elseif strcmpi(selected{iSelect},'filoBranchComplexity');
                
                % filter based only on external filoBranch network
                filoInfoFilt= filoInfo(filterSetC(:,1));
                
                filoMask1 = GCAVisualsFilopodiaMeasurementOverlays(filoInfoFilt,imgSize,...
                    'plotValues',plotValues,  'UseSmoothedCoords',ip.Results.UseSmoothedCoords,'createMask',ip.Results.createMask,...
                    'plotSigma',ip.Results.plotSigma);
                
                % plot the branches in red
                types = vertcat(filoInfoFilt(:).type);
                filoInfoBranch = filoInfoFilt(types>1);
                filoMask2 = GCAVisualsFilopodiaMeasurementOverlays(filoInfoBranch,imgSize,...
                    'UseSmoothedCoords',ip.Results.UseSmoothedCoords,'colorFiloBranch',[1,0,0],'createMask',ip.Results.createMask,...
                    'plotSigma',ip.Results.plotSigma) ;
                
                filoMask = (filoMask1 | filoMask2);
%                 
%                 NBranches = sum(types>1); % get the total number of branches
%                 
%                 lengths = vertcat(filoInfoFilt(:).Ext_length);
%                 totalLength = sum(lengths(~isnan(lengths)));
                
%                 totalLength = totalLength.*.216*10;
                %% All other cases
            else
                
                if ~isempty(regexpi(selected{iSelect},'Embedded'));
                    filoPlotType  = 2;
                    filterSetC = (filterSetC(:,1)==1 & filterSetC(:,2) ==1);
                    
                else
                    filoPlotType = 1;
                    filterSetC = filterSetC(:,1) ==1;
                end
                
                [filoMask] = GCAVisualsFilopodiaMeasurementOverlays(filoInfo,imgSize,...
                    'filoFilterSet',filterSetC,'plotValues',plotValues,...
                    'branchMode',branchMode,'colorByValue',colorByValue,'plotText',ip.Results.plotText,'justExt',...
                    filoPlotType,'extraColor',ip.Results.extraColor,'cMapLimits',cMapLimits,'plotTextAtBranches',...
                    plotTextAtBranches,'UseSmoothedCoords',ip.Results.UseSmoothedCoords,'createMask',ip.Results.createMask,...
                    'colorFiloBranch',ip.Results.colorFiloBranch,'plotSigma',ip.Results.plotSigma);
                
                if strcmpi(selected{iSelect},'filoOrientation') ;
                    overlay = zeros(imgSize);
                    backbone = veilStem(frameC).neuriteLongPathIndices;
                    overlay(backbone) = 1;
                    spy(overlay,'k');
                end
                
            end % if strcmpi
            %% Plot Veil/Stem
            hold on
            
            veilStemMask = veilStem(frameC).finalMask;
            
            roiYX = bwboundaries(veilStemMask);
            cellfun(@(x) plot(x(:,2),x(:,1),'color',ip.Results.colorVeilStem),roiYX);
            
            %% Plot Movie Extras
            
            if ip.Results.ScaleBar
                pixSizeMic = MD.pixelSize_./1000;
                plotScaleBar(ip.Results.ScaleBarSize/pixSizeMic,2,'Color',[0,0,0],'Location','SouthEast');
            end
            
            if ~isempty(ip.Results.TreatmentFrame)
                if frameC >=ip.Results.TreatmentFrame;
                    text(10,10,[ip.Results.TreatmentTitle '( + ' num2str((frameC-ip.Results.TreatmentFrame)*5) ' (s)) '],'color','k');
                end
            end
            
            if ip.Results.Timer
                text(10,30,[num2str(frameC*5-5) ' s'],'color','k',...
                    'FontSize',ip.Results.FontSizeTimer,'FontWeight','Bold');
            end
            %% Create the mask
            if ip.Results.createMask
                gcaMask = (veilStemMask | filoMask);
                gcaMask = logical(getLargestCC(gcaMask));
                
                %imwrite(gcaMask,[ outDirMask filesep 'GCAMask_' num2str(iFrame,'%03d') '.tif'])
                imwrite(gcaMask,[outDirMask filesep 'GCAMask_' num2str(frameC,'%03d') '.tif'],'compression','none');
                %                 [x] = upDirectory(outDirMask,1);
                %                 filoDir = [x filesep 'FilopodiaMasks'];
                %                 if ~isdir(filoDir)
                %                     mkdir(filoDir)
                %                 end
                %                 imwrite(filoMask,[filoDir filesep 'filoMask_' num2str(frameC,'%03d') '.tif'],'compression','none');
                %
                if ip.Results.createBackMask
                    
                    % dilate gcaMask
                    maskObjGC = imdilate(gcaMask,strel('disk',ip.Results.BMDilRad));
                    
                    
                    [maskBackEst,~,~] = gcaEstimateBackgroundArea(img,'PostProcess',true);
                    maskObjAll = imdilate(~maskBackEst,strel('disk',ip.Results.BMDilRad));
                    
                    maskBackFinal = ~(maskObjGC | maskObjAll);
                    imwrite(maskBackFinal ,[outDirMaskBack filesep 'backMask_' num2str(frameC,'%03d') '.tif'],'compression','none');
                    % for now also write sanity check overlays...
                    overlayDir = [outDirMaskBack filesep 'Overlays'];
                    if ~isdir(overlayDir)
                        mkdir(overlayDir);
                    end
                    [ny,nx] = size(img);
                    setFigure(nx,ny);
                    imshow(-img,[]);
                    hold on
                    spy(maskBackFinal,'g');
                    saveas(gcf,[overlayDir filesep num2str(frameC,'%03d') '.png']);
                    close gcf
                end
                
            end % if createMaskFlag
            %% Save Image
            if isempty(ip.Results.filoIDs)
                name = [outDir filesep num2str(frameC,'%03d')];
            else
                name = [outDir filesep num2str(frameC,'%03d') '_' num2str(ip.Results.filoIDs) ];
            end
            
            if ip.Results.screen2png
                helperScreen2png([name '.png']);
            else
                saveas(gcf,[name  '.png']);
            end
            
            if ip.Results.otherFiles
                saveas(gcf,[name  '.eps'],'psc2');
                saveas(gcf,[name '.fig']);
            end
            
            close gcf
            %% Plot Colorbar
            % color bar on figure will correspond to the image (black and
            % white) this will make a separate color bar figure
            % corresponding to the feature values that can be overlaid
            % using movie making editors such as premiere.
            if ip.Results.colorbarOverlay && colorByValue
                figure('visible',ip.Results.visible)
                test = zeros(imgSize);
                test(1,1) = cMapLimits(1);
                test(1,2) = cMapLimits(2);
                imagesc(test);
                colormap(jet(128));
                colorbar
                
                saveas(gcf,[outDir filesep 'ColorBarOverlay_' selected{iSelect} '.eps'],'psc2');
                saveas(gcf,[outDir filesep 'ColorBarOverlay_' selected{iSelect} '.fig']);
                close gcf
                
            end % colorbarOverlayFlag
        end % isempty(plotValues)
    end % for iFrames
end % for iSelect
end % function