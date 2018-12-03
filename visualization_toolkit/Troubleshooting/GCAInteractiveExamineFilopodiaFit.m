function [ output_args ] = GCAInteractiveExamineFilopodiaFit(MD,varargin)
%GCAInteractiveExamineFilopodiaFit

if nargin<1 
    [name,path] = uigetfile(pwd,'Please Select a MovieData.mat File'); 
    load([path filesep name]); 
end 

%% Check input
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('OutputDirectory',[]);
ip.addParameter('InputDirectoryFilo', []);
ip.addParameter('InputDirectoryVeil',[]);
ip.addParameter('ChannelIndexFilo',1);
ip.addParameter('ChannelIndexVeil',1);
ip.addParameter('cMapLimits',[0,10]) % in microns
ip.addParameter('embedded',false);

ip.addParameter('frames',[]);

ip.parse(varargin{:});
%%
if isempty(ip.Results.OutputDirectory)
    outDir = [MD.outputDirectory_  filesep 'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'VII_filopodiaBranch_fits' filesep ...
        'Channel_' num2str(ip.Results.ChannelIndexFilo) filesep 'Linescans'];
else
    outDir = ip.Results.OutputDirectory;
end

if ip.Results.embedded
    outDir = [outDir filesep 'Embedded_Actin_Bundle'];
else
    outDir = [outDir filesep 'Filopodia'];
end


if ~isdir(outDir)
    mkdir(outDir)
end

if isempty(ip.Results.InputDirectoryFilo)
    inDirFilo =  [MD.outputDirectory_ filesep 'SegmentationPackage' filesep 'StepsToReconstruct'...
        filesep 'VII_filopodiaBranch_fits' filesep 'Channel_' num2str(ip.Results.ChannelIndexFilo)];
else
    inDirFilo =  [ip.Results.InputDirectoryFilo filesep 'Channel_' num2str(ip.Results.ChannelIndexFilo)];
end


if isempty(ip.Results.InputDirectoryVeil)
    inDirVeil =  [MD.outputDirectory_ filesep  'SegmentationPackage' filesep 'StepsToReconstruct'...
        filesep 'III_veilStem_reconstruction' filesep 'Channel_' num2str(ip.Results.ChannelIndexVeil)];
else
    inDirVeil =  [ip.Results.InputDirectoryVeil filesep 'Channel_' num2str(ip.Results.ChannelIndexVeil)];
end

if isempty(ip.Results.frames)
    frames = [1:(MD.nFrames_ -1) ];
else
    frames = ip.Results.frames;
end

cMapLimits = ip.Results.cMapLimits;

load([inDirVeil filesep 'veilStem.mat']);
load([inDirFilo filesep 'filoBranch.mat']);

if ip.Results.embedded
    filoPart = 'Int_';
else
    filoPart = 'Ext_';
end

%%
for iFrame = 1:length(frames)
    fsFigure(0.75)
    frame = frames(iFrame);
    img = double(imread([MD.getChannelPaths{ip.Results.ChannelIndexFilo} filesep MD.getImageFileNames{ip.Results.ChannelIndexFilo}{frame}]));
    if ip.Results.embedded
        
        filoInfo = filoBranch(iFrame).filoInfo;
        x = vertcat(filoInfo(:).Int_pixIndicesBack);
        x = x(~isnan(x));
        embedded = zeros(size(img));
        embedded(x) =true;  
    end
    
    hclick = subplot(1,2,1);
    
    imshow(-img,[]);
    hold on
    filo = filoBranch(frame).reconstructInfo.output{end}.finalReconstruct;
    
    if ip.Results.embedded
        spy(filo | embedded,'r');
    else
        spy(filo,'r'); % to plot the reconstruction before fitting.
    end
    
    filoInfo = filoBranch(frame).filoInfo;
    % subplot(2,2,2);
    % Make the visual and create a corresponding mask
    %     GCAVisualsMakeFeatureOverlayMovie(MD,'cMapLimits',[0,5],'screen2png',true,'interactive',false,'FeatureDirectory',featDirVis,'ScaleBar',true);
    subplot(1,2,2);
    cMap{1} = jet(128);
  
    bodyFinal = veilStem(frame).finalMask;
    edgeYX = bwboundaries(bodyFinal);
    filoBranchC = filoBranch(frame);
    
    if ip.Results.embedded
        [filterSet,filoParams] = GCACreateFilopodiaFilterSet(filoBranch,'filterType','Validation2');
        filterSetC = filterSet(frame);
        pixSizeMic = MD.pixelSize_/1000;
        filoLengths = GCAAnalysisExtract_filoLength(filoBranchC,filterSetC,'umPerPixel',pixSizeMic,'FiloPart','Tot','outPercent',false);
        plotValues = filoLengths{1};
        plotText{1} = false;
        plotText{2} = false;
        filterSetC= filterSetC{1};
        filterFrameC = filterSetC(:,1);
        
        filterInt = (filterSetC(:,1) == 1 & filterSetC(:,2) ==0 ); % get the ID of all non-fits internally this filter is the length of
        %             % the original filoInfo detection
        filterInt = filterInt(filterFrameC); % keep only the
        filoInfoExtBund = filoInfo(filterFrameC);
        filoInfoIntBund = filoInfoExtBund(~filterInt);
        imgSize = size(img);
        
        filoInfoFilt{1} = filoInfoExtBund; % should be the same size as the numbers
        filoInfoFilt{2} = filoInfoIntBund;
        
        plotValuesSub{1} = plotValues;
        plotValuesSub{2} = plotValues(~filterInt);
        
        imshow(-img,[]);
        hold on
        cellfun(@(x) plot(x(:,2),x(:,1),'k'),edgeYX);
        hold on
        
        for i = 1:2
            
            filoMaskCell{i} =    GCAVisualsFilopodiaMeasurementOverlays(filoInfoFilt{i},imgSize, ...
                'plotValues',plotValuesSub{i},'colorByValue',true,'plotText',plotText{i},'justExt',i,'cMapLimits',cMapLimits);
        end
    else
        [filoFilterSet,filoParams] = GCACreateFilopodiaFilterSet(filoBranch,'filterType','Validation_NoEmbed_NoSave');
        
        filoFilterSetC = filoFilterSet(frame);
        pixSizeMic = MD.pixelSize_/1000;
        filoLengths = GCAAnalysisExtract_filoLength(filoBranchC,filoFilterSetC,'umPerPixel',pixSizeMic);
        
        plotValues = filoLengths{1};
        
        
        for iMap = 1:1
            
            
            imshow(-img,[]);
            hold on
            cellfun(@(x) plot(x(:,2),x(:,1),'k'),edgeYX);
            
            
            imgSize = size(img);
            GCAVisualsFilopodiaMeasurementOverlays(filoInfo,imgSize,...
                'filoFilterSet',filoFilterSetC{1},'plotValues',plotValues,...
                'justExt',1,'cMapLimits',cMapLimits,'colorByValue',true, ...
                'plotText',false,'colorMap',cMap{iMap});
            
            %     if ip.Results.plotScaleBar
            %         plotScaleBar(pixels,pixels/10,'Color',textColor);
            %     end
            
            
        end % iMap
    end
    clickHappy = true;
    while clickHappy == true
        reply1 = questdlg(['View Filopodium Intensity Decay for Frame ' num2str(frame,'%03d') '?']);
        
        switch reply1
            case 'Yes'
                
                h=impoint(hclick);
                position = wait(h);
                idx = sub2ind(size(img),round(position(2)),round(position(1)));
                %     coords(ptCount,2) = position(2);
                %     coords(ptCount,1) = position(1);
                x = position(1);
                y= position(2);
                
                [ filoIdx,~ ] = getFiloIDFromCoords(filoInfo, x,y );

                run = 1;
                if ~isempty(filoIdx)
                    if ip.Results.embedded
                        run = ~(sum(isnan(filoInfo(filoIdx).Int_pixIndicesBack))>0 && length(filoInfo(filoIdx).Int_pixIndicesBack)==1);
                    end
                    if run ==true;
                        %%
                        
                        GCAVisualsPlotFilopodiaFit( img,filoInfo,'idxToPlot',filoIdx,'FiloPart',filoPart);
                        reply3 = questdlg('Save FilopodiaFit Information?');
                        
                        if strcmpi(reply3,'yes')
                            saveas(gcf,[outDir filesep 'FilopodiaFit_Frame_' num2str(frame,'%03d') '_' 'Filopodia_' num2str(filoIdx,'%03d') '.png']);
                            saveas(gcf,[outDir filesep 'FilopodiaFit_Frame_' num2str(frame,'%03d') '_' 'Filopodia_' num2str(filoIdx,'%03d') '.fig']);
                            saveas(gcf,[outDir filesep 'FilopodiaFit_Frame_' num2str(frame,'%03d') '_' 'Filopodia_' num2str(filoIdx,'%03d') '.eps'],'psc2');
                            close(gcf)
                        else
                            close(gcf)
                        end
                    else
                        uiwait(msgbox('Could not ID an embedded actin bundle: Please try again'));
                    end
                else
                    uiwait(msgbox('Could not ID a filopodia: Please try again'));
                end
            case 'No'
                clickHappy = false;
            case 'Cancel'
                
                break
        end % if switch reply
    end % while
    if strcmpi(reply1,'Cancel')
        close(gcf)
        display('User Terminated: Relaunch to Continue');
        break
    end
end % iFrame
close all
end % Function