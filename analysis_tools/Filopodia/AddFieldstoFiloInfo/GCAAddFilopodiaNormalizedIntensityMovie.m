function [ output_args ] = GCAAddFilopodiaNormalizedIntensityMovie(movieData,varargin)
%GCAFilopodiaNormalizedIntensityMovie(movieData)

% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS

ip.addParameter('InputDirectory',[]); % Where the filopodia information is saved
ip.addParameter('InputDirectoryVeilStem',[]); % Where the VeilStem information is saved
ip.addParameter('ChannelIndex',1); % channel index for the filopodia channel 
ip.addParameter('ChannelIndexVeil',1); % channel index for the veil channel 

ip.addParameter('frames',[]);

ip.addParameter('saveNormFactor',false); % flag to save a list of the norm
% factors per frame.
ip.addParameter('FeatureDirectory',[]);

ip.addParameter('UseBackSubtractImages',true); % flag to use the background subtracted images 

ip.parse(varargin{:});

%% Initiate
if isempty(ip.Results.InputDirectory)
    inDirFilo = [movieData.outputDirectory_ filesep 'SegmentationPackage' ...
        filesep 'StepsToReconstruct' filesep...
        'VII_filopodiaBranch_fits' ];
else
    inDirFilo = ip.Results.InputDirectory;
end

if isempty(ip.Results.InputDirectoryVeilStem)
    inputDirVeilStem = [movieData.outputDirectory_ filesep ...
        'SegmentationPackage' filesep 'StepsToReconstruct' filesep ...
        'III_veilStem_reconstruction'];
else
    inputDirVeilStem = ip.Results.InputDirectoryVeilStem;
end

load([inDirFilo filesep 'Channel_' num2str(ip.Results.ChannelIndex) filesep  'filoBranch.mat']);
load([inputDirVeilStem filesep 'Channel_' num2str(ip.Results.ChannelIndexVeil) filesep 'veilStem.mat']);
%
if isempty(ip.Results.frames)
    frames = 1:length(filoBranch);
else
    frames = ip.Results.frames;
end
%%
if ip.Results.UseBackSubtractImages
    imDir = [movieData.outputDirectory_ filesep 'GCABackSubtract' filesep 'BackgroundSubtracted_Images' filesep 'Channel_1'];
    imgFilenames = searchFiles('.tif',[],imDir,0,'all',1);
    toAdd = 'BS'; 
else % use the channel data
        imgDir = [movieData.channels_(ip.Results.ChannelIndex).channelPath_];
        imgFilenames = searchFiles('.tif','memo',imgDir,0,'all',1);
        toAdd = [];
end
%% START
for iFrame = 1:length(frames)
    
    % extract the veil
    veilMask = veilStem(frames(iFrame)).finalMask;
    
    img = double(imread(imgFilenames{frames(iFrame)}));
    
    % extract the filo info to read into the function
    filoInfo = filoBranch(frames(iFrame)).filoInfo;
    % add the metric to the filo info - NOTE in the future might want to just
    % calculate automaticaly at the time of fitting to be more efficient.
    [filoInfo,normFactPerFrame] = GCAAddFilopodiaNormalizedIntensity(img,veilMask,filoInfo);
    filoBranch(frames(iFrame)).filoInfo = filoInfo;
    measC{frames(iFrame)} = normFactPerFrame;
end
% resave the values
save([inDirFilo filesep 'Channel_' num2str(ip.Results.ChannelIndex) filesep  'filoBranch.mat'],'filoBranch','-v7.3');

%% save the normalization factor in the measurements folder
if ip.Results.saveNormFactor
    if isempty(ip.Results.FeatureDirectory)
        expFolder = [movieData.OutputDirectory_ filesep 'GCAFeatureExtraction' filesep 'Descriptor' filesep  toAdd ];
    else
        expFolder = ip.Results.FeatureDirectory;
    end
    
    if ~isdir(expFolder)
        mkdir(expFolder);
    end
    save([expFolder filesep 'meas_ExpressionNormalization.mat'],'measC');
end % ip.Results.saveNormFactor
end % function