function [ output_args ] = GCARunFullPipeline(varargin)
% GCARunFullPipeline 
%% Input 
ip = inputParser;

ip.CaseSensitive = false;


ip.addParameter('ChannelDirectory',[]); % directory where your .tif files are located: 
% if empty will ask you to load 

ip.addParameter('OutputDirectory',[]); % directory where GCA will save your analysis: 
% if empty will ask you to load 

ip.addParameter('parameterMFile',[]); % 
% if empty will ask you to load : otherwise enter full path name 

ip.addParameter('multChannels',false); 

ip.addParameter('makeMovieData',true); 

ip.addParameter('runVeilStem',true); 

ip.addParameter('runFilopodia',true); 

ip.addParameter('runFeatureExtract',true); 

ip.addParameter('runWindAnalysis',false); 

ip.parse(varargin{:});
%% 
    x =pwd;
if ip.Results.makeMovieData

    if isempty(ip.Results.ChannelDirectory)
        channelDir = uigetdir(x,'Please Select a Channel Directory');
    else
        channelDir = ip.Results.ChannelDirectory;
    end
end

if ip.Results.makeMovieData
    if isempty(ip.Results.OutputDirectory)
        if ip.Results.multChannels
            %outDir = uigetdir(x,'Please Select an Output Directory');
            x = upDirectory(channelDir,1);
        else
            x = upDirectory(channelDir,2);
        end 
            outDir = [x filesep 'GrowthConeAnalyzer'];        
    else
        outDir = [ip.Results.OutputDirectory filesep 'GrowthConeAnalyzer'];
    end

if ~isdir(outDir)
    mkdir(outDir)
end 
else 
   [MDname, path] =  uigetfile(pwd,'Please Select a movieData.mat File' ); 
   load([path filesep MDname]); 
    
end

    if isempty(ip.Results.parameterMFile)
        
        [pMFile,pathName] = uigetfile(pwd,'Select a Parameter.m File');
       
    else
        fullPath = ip.Results.parameterMFile;
        [pathName,pMFile] = upDirectory(fullPath,1,1);         
    end
    addpath(genpath(pathName)); 

%% Create the movie data object 

if ip.Results.makeMovieData
    movieInfoMFile = strrep(pMFile,'.m',''); 
    toRun= str2func(movieInfoMFile);
    m = toRun('step2load','movie');
    
    if ip.Results.multChannels
        f = dir(channelDir);
        f = f(3:end);
        % make the multiple channel objects
        for iCh = 1:length(f)
            
            channel(iCh) = Channel([channelDir filesep f(iCh).name]);
        end
    else    % only a single channel object
        
        channel = Channel(channelDir);
    end
    % Constructor needs an array of channels and an output directory (for analysis)
    MD = MovieData(channel,outDir);
    % Set the path where to store the MovieData object.
    MD.setPath(outDir);
    MD.setFilename('movieData.mat');
    % Run sanityCheck on MovieData.
    % Check image size and number of frames are consistent.
    % Save the movie if successfull
    MD.sanityCheck;
    
    fields = fieldnames(m);
    for ifield = 1:numel(fields)
        % Set some additional movie properties
        %  MD.numAperture_=1.4;
        MD.(fields{ifield})= m.(fields{ifield});
    end
    MD.save;
end
%% Run the Veil/Stem Segmentation

if ip.Results.runVeilStem
    segParameterMFile = strrep(pMFile,'.m','');
    GCARunVeilStemSingleMovie(MD,'segParameterMFile',segParameterMFile);
end
%% Run the Filopodia Segmentation
if ip.Results.runFilopodia
    segParameterMFile = strrep(pMFile,'.m','');
    GCARunFilopodiaSingleMovie(MD,'segParameterMFile',segParameterMFile);
end

%% Run the Feature Extraction
if ip.Results.runFeatureExtract
    featureFile = strrep(pMFile,'.m','');
    GCARunFeatureExtractionSingleMovie(MD,'featureParameterMFile',featureFile);
end
%%  Run the Veil Velocity Analysis
if ip.Results.runWindAnalysis
    segParameterFile = strrep(pMFile,'.m','');
    GCARunVeilVelocityAnalysisSingleMovie(MD,'segParameterMFile',segParameterFile);
end
end % function  