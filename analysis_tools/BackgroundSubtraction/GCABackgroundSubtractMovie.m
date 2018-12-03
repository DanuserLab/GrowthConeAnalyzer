function movieData = backgroundSubtractMovie(movieData,varargin)
%BACKGROUNDSUBTRACTMOVIE corrects for background by subtracting the average background value
%
% movieData = backgroundSubtractMovie(movieData)
%
% movieData = backgroundSubtractMovie(movieData,paramsIn)
%
% This function performs background subtraction on the movie described by
% the input movieData. This is accomplished by averaging the intensity in
% the areas covered by background masks, as created using
% createMovieBackgroundMasks.m. This average value is determined for each
% frame and then subtracted from each pixel in that frame. Negative values
% are converted to zero. 
% 
% 
% Input:
% 
%   movieData - The MovieData object describing the movie, as created using
%   setupMovieDataGUI.m
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the corrected images to.
%       Corrected images for different channels will be saved as
%       sub-directories of this directory. If not input, the corrected
%       images will be saved to the same directory as the movieData, in a
%       sub-directory called "background_subtracted_images"
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) to perform shade correction on. This
%       index corresponds to the channel's location in the array
%       movieData.channels_. If not input, all channels will be corrected.
%      
%       ('MaskChannelIndex' -> Positive integer scalar or vector)
%       This parameter specifies the channels to use masks from when
%       performing background subtraction on each channel. This allows
%       masks generated for one channel to be used in performing background
%       subtraction on other channels, which may or may not themselves have
%       masks. This vector or scalar must be the same size as ChannelIndex.
%       Optional. If not input, each channel will use it's own masks. 
%
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output and user
%       interaction is suppressed.
%
%
%
% Output:
%
%   movieData - the updated movieData object with the correction
%   parameters, paths etc. stored in it, in the field movieData.processes_.
%
%   The corrected images are written to the directory specified by the
%   parameter OuptuDirectory, with each channel in a separate
%   sub-directory. They will be stored as bit-packed .tif files. 
%
% 
% Hunter Elliott, 11/2009
% Revamped 5/2010
%

%%  --------- Parameters ------- %%

pString = 'bs_'; %The string to prepend before the background-subtracted image directory & channel name
saveName = 'background_subtraction_values_for_channel_'; %File name for saving subtracted values
dName = 'background_subtracted_images_for_channel_';%String for naming the directories for each corrected channel

%% ----------- Input ------------ %%

%Check that input object is a valid moviedata
if ~isa(movieData,'MovieData')
    error('The first input argument must be a valid MovieData object!')
end

%%Input check
ip = inputParser;

ip.CaseSensitive = false;
%Get the indices of any previous background subtraction processes from this
ip.addParameter('ChannelIndex',1);
ip.addParameter('OutputDirectory',[]);
ip.addParameter('BackMaskDirectory',[]);

ip.parse(varargin{:});


%% ---------- Init ---------- %%

if isempty(ip.Results.OutputDirectory)
    outDir = [movieData.outputDirectory_ filesep 'GCABackSubtract' filesep 'BackgroundSubtracted_Images'];
else
    outDir = ip.Results.OutputDirectory;
end

if isempty(ip.Results.BackMaskDirectory)
    backMaskDir = [movieData.outputDirectory_ filesep 'GCABackSubtract' filesep 'BackgroundMasks'];  
else
    backMaskDir = ip.Results.BackMaskDirectory;
end

listOfBackMaskPaths = searchFiles('.tif','',backMaskDir,0,'all',1);
nChan = length(ip.Results.ChannelIndex);
%% ---- Background Subtraction ---- %%
%Go through each image and subtract the average value behind the background
%mask from the image.

backgroundValues = cell(1,nChan);

%Go through each image and apply the appropriate shade correction
for iChan = 1:nChan
    %outDirC = [outDir filesep p.ChannelIndex(iChan) ];
    %     inDir  = movieData.processes_{iProc}.inFilePaths_{1,p.ChannelIndex(iChan)};
    %     outDir = movieData.processes_{iProc}.outFilePaths_{1,p.ChannelIndex(iChan)};
    %     corrDir = movieData.processes_{iProc}.inFilePaths_{2,p.ChannelIndex(iChan)};
    channelC = ip.Results.ChannelIndex(iChan);
    outDirC = [outDir filesep 'Channel_' num2str(channelC)];
    
    if ~isdir(outDirC)
        mkdir(outDirC)
    end
    %
    %     disp(['Background subtracting channel ' num2str(p.ChannelIndex(iChan)) '...'])
    %     disp(['Background subtracting images from "' inDir '", results will be stored in ' outDir]);
    %     disp(['Using background masks from ' corrDir])
    
    nImages = numel(listOfBackMaskPaths); % note for now use the backMask paths as truncated the original movie when performing the filo/branch
% reconstruct and masks 
    backgroundValues{iChan} = zeros(1,nImages);
    
    for iImage = 1:nImages
        
        % currIm = imread([inDir filesep inNames{iChan}{iImage}]);
        currIm =  imread([movieData.getChannelPaths{channelC} filesep movieData.getImageFileNames{channelC}{iImage}]);
        
        ogClass = class(currIm); %Determine class so it is not altered later
        
        %Load the background mask
        %currBackMask = imread([corrDir filesep bakNames{1}{iImage}]);
        currBackMask = imread([listOfBackMaskPaths{iImage}]);
        
        %Get average background intensity
        backgroundValues{iChan}(iImage) = ...
            mean(double(currIm(currBackMask(:))));
        currIm = double(currIm) - backgroundValues{iChan}(iImage);
        currIm(currIm < 0) = 0; %Clip negative values to zero
        currIm = cast(currIm,ogClass); %Return to original data type
        
        %Write it to disk
        %         imwrite(currIm,[outDir filesep pString ...
        %             inNames{iChan}{iImage}]);
        imwrite(currIm,[outDirC filesep num2str(iImage,'%03d') '.tif']);
    end
end

save([outDirC filesep 'backgroundValues.mat'],'backgroundValues');
setAxis
plot(backgroundValues{1}); 
ylabel('Mean Background Intensity Value (AU)'); 
xlabel('Frames'); 

saveas(gcf,[outDirC filesep 'meanBackgroundOverTime.fig']); 
saveas(gcf,[outDirC filesep 'meanBackgroundOverTime.png']); 
close gcf
end
