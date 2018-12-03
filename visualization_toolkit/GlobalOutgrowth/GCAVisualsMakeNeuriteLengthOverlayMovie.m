function [ output_args ] = GCAVisualsMakeNeuriteLengthOverlayMovie(MD,varargin)
%GCAVisualsMakeNeuriteLengthOverlayMovie
% Make neurite length overlay movie as in Video9 and Video10 of GCA manuscript
% 




%% Input check
ip = inputParser;

defaultParDir = [MD.outputDirectory_ filesep 'GCAFeatureExtraction' filesep 'Partition_Outgrowth_Trajectory'];
ip.addParameter('ChannelIndexVeil',1);
ip.addParameter('OutputDirectory',[]);
ip.addParameter('VeilDirectory',[],@(x) ischar(x) || isempty(x));

ip.addParameter('colorByOutgrowthState',false); % flag to color the neurite length 
% by outgrowth state- otherwise will be color specified by 'colorLength'
% and 'colorVeil' parameters
% red: accelerate 
% pink: decelerate
% cyan: paused
% blue: retracting

ip.addParameter('PartitionDirectory',defaultParDir, @(x) ischar(x));

ip.addParameter('firstFrameLimits',false);
ip.addParameter('Timer',false);
ip.addParameter('FontTimerSize',18);
ip.addParameter('ScaleBar',false);

ip.addParameter('colorVeil','k'); 
ip.addParameter('colorLength','k');

ip.addParameter('writeStateName',true);

ip.addParameter('screen2png',true);
ip.parse(varargin{:});
%% Initiate

if isempty(ip.Results.OutputDirectory)
    outDir = [MD.outputDirectory_ filesep 'GCAFeatureExtraction' filesep 'GlobalFunctional' filesep 'neurite_outgrowth_measurements' ...
        filesep 'neuriteOutgrowth_FeatureMovie'];
else
    outDir = ip.Results.OutputDirectory;
end

partitionDir = ip.Results.PartitionDirectory;
if ~isdir(partitionDir) && ip.Results.colorByOutgrowthState
    colorByOutgrowthState = false; % force to be false
    display('No neurite trajectory directory found: Cannot color by Outgrowth State');
else
    colorByOutgrowthState = ip.Results.colorByOutgrowthState;
end

if isempty(ip.Results.VeilDirectory)
    
    neuriteOutgrowthDir  = [MD.outputDirectory_ filesep 'SegmentationPackage' filesep 'StepsToReconstruct' filesep ...
        'IV_veilStem_length' filesep 'Channel_' num2str(ip.Results.ChannelIndexVeil) ];
else
    neuriteOutgrowthDir = ip.Results.VeilDirectory;
end

if colorByOutgrowthState
    load([partitionDir filesep 'globalMeas.mat']);
    vels = horzcat(globalMeas.outgrowth.groupedVelSmoothed{:});
    states = globalMeas.outgrowth.stateAllFrames;
end

load([neuriteOutgrowthDir filesep 'veilStem.mat']);

if ~isdir(outDir)
    mkdir(outDir)
end

ny = MD.imSize_(1);
nx = MD.imSize_(2);
%% Start
for iFrame = 1:MD.nFrames_
    setFigure(nx,ny,'off');
    img =double(imread([MD.getChannelPaths{1} filesep MD.getImageFileNames{1}{iFrame}]));
    
    if  iFrame ==1
        if ip.Results.firstFrameLimits
            lims = [min(-img(:)) max(-img(:))];
        else
            lims = [];
        end
    end
    
    % load longPath mask
    mask = zeros(MD.imSize_);
    mask(veilStem(iFrame).neuriteLongPathIndices)= 1;
    
    imshow(-img,lims);
    hold on
    if colorByOutgrowthState
        stateC = states(iFrame);
        switch stateC
            case 1
                c = 'c';
                colorVeil = 'c';
                stateName = 'Pause';
            case 2
                c = 'b';
                colorVeil = 'b';
                stateName = 'Retract';
            case 3
                c = 'r';
                colorVeil = 'r';
                stateName = 'Growth : Accelerate';
            case 4
                c = 'm';
                colorVeil = 'm';
                stateName = 'Growth : Decelerate';
        end
        
        if ip.Results.writeStateName
            text(20,20,{'Outgrowth State: '; stateName}, 'FontSize',16);
        end
    else
        colorVeil = ip.Results.colorVeil;
        c = ip.Results.colorLength ;
    end
    veilStemMask = veilStem(iFrame).finalMask;
    roiYX = bwboundaries(veilStemMask);
    cellfun(@(x) plot(x(:,2),x(:,1),'color',colorVeil),roiYX);
    
    spy(mask,c);
    
    if ip.Results.Timer
        text(nx-70,ny-20,[num2str(iFrame*MD.timeInterval_ - MD.timeInterval_),' s'] ,'Color', 'k',...
            'FontSize',ip.Results.FontTimerSize,'FontWeight','Bold');
    end
    
    if ip.Results.ScaleBar
        pixSizeMic = MD.pixelSize_/1000;
        width  = 10/pixSizeMic;
        plotScaleBar(width,2,'Location','SouthWest','Color',[0 0 0]);
    end
    
    if ip.Results.screen2png
        helperScreen2png([outDir filesep num2str(iFrame,'%03d') '.png']);
    else
        saveas(gcf,[outDir filesep num2str(iFrame,'%03d') '.png']);
    end
    saveas(gcf,[outDir filesep num2str(iFrame,'%03d') '.eps'],'psc2');
    saveas(gcf,[outDir filesep num2str(iFrame,'%03d') '.fig']);
    close gcf
end % for iFrame
end % function