function [ output_args ] = GCAAnalysisPartitionTimeSeriesMovie(MD,varargin)
%GCAAnalysisPartitionTimeSeriesMovie 
% Movie Wrapper for GCAAnalysisPartitionTimeSeries and 
% GCAfindPausingInNeuriteOutgrowthTrajectory 
% collect local params partitioned if user requires. 
% INPUT: 
% Pointer to a folder holding 
% neuriteLength: 
% rx1 double containing the neurite lengths at each frame (output of STEP
% IV of the GCA PACKAGE: in um) 
%% CHECK Input

ip = inputParser;
ip.addParameter('OutputDirectory',[]);
ip.addParameter('splineParam',0.01, @(x) isscalar(x) || isempty(x) ); 
ip.addParameter('threshPause',0.5, @isscalar); 

ip.addParameter('makePlots',true); % % flag to make neurite outgrowth 
% plots colored by outgrowth state as shown in Fig. of the manuscript 
ip.addParameter('title',[]); 

ip.parse(varargin{:});
%% Initiate 
if isempty(ip.Results.OutputDirectory)
    outDir  = [MD.outputDirectory_ filesep 'GCAFeatureExtraction' filesep ...
        'Partition_Outgrowth_Trajectory'];
else
    outDir = ip.Results.OutputDirectory;
end


splineParam = ip.Results.splineParam; 
threshPause = ip.Results.threshPause; 

% load the neurite Length array 
outgrowthFile = [MD.outputDirectory_ filesep 'GCAFeatureExtraction' ... 
        filesep 'GlobalFunctional' filesep 'neurite_outgrowth_measurements'... 
        filesep 'neuriteLengthOutput.mat'];
% Set up the output directory   
if exist(outgrowthFile)~= 0 ;
    load(outgrowthFile);
    
    if ~isdir(outDir);
        mkdir(outDir);
    end
    
    % Get up titles for plots if necessary 
    if isempty(ip.Results.title) % flag for automatic title based on filenames 
        title = gcaGetNeuriteID(MD.outputDirectory_);
        
        % just make sure no underscores
        title = strrep(title,'_',' ');
    else 
        title = ip.Results.title; 
    end
    %% MAIN 
    % Partition based on pausing in the neurite
    % trajectory (when the velocity of outgrowth reaches below a user
    % selected threshold)
    globalMeas =  GCAfindPausingInNeuriteOutgrowthTrajectory(neuriteLength,'outPath',outDir,... 
        'forTitle',title,'splineParam',splineParam,'threshPause',threshPause,... 
        'secPerFrame',MD.timeInterval_,'makePlot',ip.Results.makePlots) ;
else
    display('Please Run the Neurite Outgrowth Measurement to Continue');
end
end %function 