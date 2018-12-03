function [ output_args ] = GCARunGlobalSegmentationProcess( MD ,varargin)
% GCARunGlobalSegmentationProcess
%% Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
% Restart Options 
ip.addParameter('method','otsu');
ip.addParameter('GaussFilterSigma',1); 
ip.addParameter('OutputDirectory',[]); 
ip.addParameter('ChannelIndex',1); 
ip.parse(varargin{:});

%% 
if isempty(ip.Results.OutputDirectory)
    outDir = [MD.outputDirectory_ filesep 'SegmentationPackage' filesep ...
        'StepsToReconstruct' filesep 'III_veilStem_reconstruction' filesep 'InitialGlobalThresholding'];
else 
    outDir = ip.Results.OutputDirectory; 
end 

if ~isdir(outDir)
    mkdir(outDir)
end 
    
switch ip.Results.method
    case 'minMax'
        m = 1; 
    case 'otsu' 
        m = 2; 
    case 'Rosin' 
        m = 3; 
    case 'gradient'
        m = 4; 
end 

    thresProc = ThresholdProcess(MD);
    %
    % % % Save the process in the movie object
    MD.addProcess(thresProc);
 
    % Create a segmentation package
%     segPackage = SegmentationPackage(MD);
    
    %  Save the package in the movie object
%     MD.addPackage(segPackage);
    %
    % % Associate the threshold process to the package
%     MD.packages_{1}.setProcess(1,thresProc);
    
    params = MD.processes_{end}.funParams_;
    
    params.MethodIndx = m;
    params.GaussFilterSigma = ip.Results.GaussFilterSigma;
    params.OutputDirectory = outDir; 
    params.ChannelIndex = ip.Results.ChannelIndex; 
    
    parseProcessParams(MD.processes_{end},params);
    %
    % % % Run the process
    MD.processes_{end}.run; % run the gradient based
    
    MD.save

end

