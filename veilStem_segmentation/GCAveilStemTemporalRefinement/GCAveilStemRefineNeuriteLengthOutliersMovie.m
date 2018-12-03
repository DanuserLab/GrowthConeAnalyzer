function [ output_args ] = GCAveilStemRefineNeuriteLengthOutliersMovie( movieData,varargin)
%%% GCAveilStemRefineNeuriteLengthOutliersMovie :  movieDataWrapper function
% STEP V (Optional) in Veil/Stem estimation of GCA Segmentation
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%
%   movieData (REQUIRED)  - The MovieData object describing the movie, as created using
%   movieSelectorGUI.m
%
% Generic Fields:
%       ('Input/Output Fields Needed for Wrapper' -> Possible Values)
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the backboneInfo structure to.
%       If not input, the backboneInfo will be saved in the same directory
%       as the movieData, in a sub-directory called
%       "neurite_orientation"
%       % PERSONAL NOTE : You might need
%       to truncate these names for the windows paths. %
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) on which to perform the neurite orientation
%       estimation.
%       This index corresponds to the channel's location in the array
%       movieData.channels_. If not input, all channels will be analyzed
%
%       ('ProcessIndex' -> Positive integer scalar or vector)
%       This parameter specifies the output process to use for performing the
%       estimation
%       This allows a previously processed image (ie shade corrected or
%       background subtracted to potentially be input). If not input, the
%       backbone information will be calculated from the channels (ie raw
%       images)
%
%       ('Restart'  -> structure with fields = 'auto';   % 'auto'
%       paramsIn.restart.endFrame = 'auto'; % 'auto' or number, default auto.
%       paramsIn.plots = 1;
%
%  SPECIFIC INPUT
%% INPUTPARSER
% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('InputDirectory' ,[]);
ip.addParameter('OutputDirectory',[]);
ip.addParameter('ChannelIndex',1);
ip.addParameter('ProcessIndex',0);

% Specific
ip.addParameter('TSOverlays',true,@(x) islogical(x));

ip.addParameter('nFramesBack',5);
ip.addParameter('nFramesFor',5);
ip.addParameter('nFramesPixOutlier',2);
ip.addParameter('redilationRadius',4);

ip.parse(varargin{:});
p = ip.Results;
%% Initiate
nChan = numel(p.ChannelIndex);

if isempty(ip.Results.OutputDirectory)
    outDir = [movieData.outputDirectory_ filesep...
    'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'V_veilStemOutlierRefinement'];
else
    outDir = ip.Results.OutputDirectory; 
end

%% Start Channel Wrapper

for iCh = 1:nChan
    
    display(['Refining Segmentation of Outlier Frames for Channel ' num2str(p.ChannelIndex(iCh))]);
    %% Get Start and End Frames Based on Restart Choice
    
    % make final output dir where backboneInfo will be saved
    outDirC =  [outDir filesep 'Channel_' num2str(p.ChannelIndex(iCh))];
    
    if ~isdir(outDirC)
        mkdir(outDirC);
    end
    %%    Load Veil/Stem information with outlier detection
    
    if isempty(ip.Results.InputDirectory)
        inDir = [movieData.outputDirectory_ filesep  'SegmentationPackage'...
            filesep 'StepsToReconstruct' filesep 'IV_veilStem_length' filesep 'Channel_' num2str(p.ChannelIndex(iCh))];
    else
        inDir = ip.Results.InputDirectory;
    end
    
    load([inDir filesep 'veilStem.mat']);
    %% Main Function:
    
    veilStem = GCAveilStemRefineNeuriteLengthOutliers(veilStem,'nFramesFor',ip.Results.nFramesFor,...
        'nFramesBack',ip.Results.nFramesBack,'nFramesPixOutlier',...
        ip.Results.nFramesPixOutlier,'redilationRadius',ip.Results.redilationRadius);
    
    save([outDirC filesep 'veilStem.mat'],'veilStem');
    
end % for iCh
end % function