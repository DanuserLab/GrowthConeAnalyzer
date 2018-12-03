function [ output_args ] = GCACreateBackgroundMaskMovie(movieData,varargin)
% GCACreateBackgroundMaskMovie


%% 
ip = inputParser;
ip.CaseSensitive = false;

ip.addParameter('OutputDirectory',[]); 
ip.addParameter('InputDirectoryFiloBranch',[]); 
ip.addParameter('InputDirectoryVeilStem',[]); 

ip.addParameter('ChannelIndex',1); % vector of channel indexes 
ip.addParameter('ChannelIndexVeil',1); 

ip.addParameter('useFiloDetectBackEst',true); % load and use the same dilation parameters 
% that were originally used for the filopodia reconstruction. 

ip.addParameter('TSOverlays',true); 
ip.addParameter('maskToWrite','backMask'); % Options: backMask, or beforeDilate
ip.CaseSensitive = false;
ip.parse(varargin{:});

%% Set up
nChannels = length(ip.Results.ChannelIndex);

if isempty(ip.Results.InputDirectoryFiloBranch)
    inDirFiloBranch = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep 'StepsToReconstruct' filesep ...
        'VI_filopodiaBranch_reconstruction'];
else
    inDirFiloBranch = ip.Results.InputDirectoryFiloBranch;
end

if isempty(ip.Results.InputDirectoryVeilStem)
    inDirVeilStem = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep 'StepsToReconstruct' filesep ...
        'III_veilStem_reconstruction'];
else
    inDirVeilStem = ip.Results.InputDirectoryVeilStem;
end

if isempty(ip.Results.OutputDirectory)
    switch ip.Results.maskToWrite
        case 'backMask'
            outDir = [movieData.outputDirectory_ filesep 'GCABackSubtract' filesep 'BackgroundMasks'];
        case 'beforeDilate'
            outDir = [movieData.outputDirectory_ filesep 'GCABackSubtract' filesep 'ForBackgroundMasks'];
    end
else
    outDir = ip.Results.OutputDirectory;   
end

if ~isdir(outDir)
    mkdir(outDir);
end

%% Start

for iCh = 1:nChannels
    cChannel= ip.Results.ChannelIndex(iCh);
    
    load([inDirFiloBranch filesep 'Channel_' num2str(cChannel) filesep 'filoBranch.mat']);
 
    if ip.Results.useFiloDetectBackEst
        load([inDirFiloBranch filesep 'Channel_' num2str(cChannel) filesep 'params.mat']);
    end
    
    load([inDirVeilStem filesep 'Channel_' num2str(ip.Results.ChannelIndexVeil) filesep 'veilStem.mat']);
    
    imgSize = movieData.imSize_;
    nFrames = length(filoBranch); % truncate the filobranch sometimes... so don't use the original nFrames stored in movieData
    for iFrame = 1:nFrames
        img = double(imread([movieData.getChannelPaths{cChannel} filesep movieData.getImageFileNames{cChannel}{iFrame}])); 
        filoInfo = filoBranch(iFrame).filoInfo;
        veilStemMask = veilStem(iFrame).finalMask;
 
        [backMask,beforeDilate] =  GCACreateBackgroundMask(veilStemMask,filoInfo,img,'dilateLocalRegion',p(iFrame).dilateLocalRegion, 'LRDilRad',p(iFrame).LRDilRad);
       
        switch ip.Results.maskToWrite
            case  'backMask' 
        imwrite(backMask,[outDir filesep num2str(iFrame,'%03d') '.tif']);
        
        if ip.Results.TSOverlays
            overlayDir = [outDir  filesep 'Overlays'];
            if ~isdir(overlayDir);
                mkdir(overlayDir);
            end
            
            setFigure(imgSize(2),imgSize(1));
            imshow(-img,[]);
            hold on
            spy(backMask,'r');
            saveas(gcf,[overlayDir filesep num2str(iFrame,'%03d') '.fig']);
            saveas(gcf,[overlayDir filesep num2str(iFrame,'%03d') '.png']);
        end  % if ip.Results.TSOverlays 
  %%      
            case 'beforeDilate'
                imwrite(beforeDilate,[outDir filesep num2str(iFrame,'%03d') '.tif']);
                
                if ip.Results.TSOverlays
                    overlayDir = [outDir  filesep 'Overlays'];
                    if ~isdir(overlayDir);
                        mkdir(overlayDir);
                    end
                    
                    setFigure(imgSize(2),imgSize(1));
                    imshow(-img,[]);
                    hold on
                    spy(beforeDilate,'r');
                    saveas(gcf,[overlayDir filesep num2str(iFrame,'%03d') '.fig']);
                    saveas(gcf,[overlayDir filesep num2str(iFrame,'%03d') '.png']);
                    % save an extra file to fool the movieData biosensors
                    % package as the output of the GCAReconstruct for now is
                    % nFrames -1 (due to matching the protrusion normals)
                    if iFrame == nFrames
                         imwrite(beforeDilate,[outDir filesep num2str(iFrame+1,'%03d') '.tif']);
                    end 
                    
                end  % if ip.Results.TSOverlays
  
        end % switch 

    end % for iFrame
end % for iChannel
