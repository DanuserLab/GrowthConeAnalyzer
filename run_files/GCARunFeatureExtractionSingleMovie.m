function [ output_args ] = GCARunFeatureExtractionSingleMovie(MD,varargin)
% GCARunFeatureExtractionSingleMovie

%% Input
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('featureParameterMFile',[]);

ip.parse(varargin{:});
%% Perform the Background Subtraction Using the Segmentation Information
if isempty(ip.Results.featureParameterMFile)
    % if myParamsName empty just run the defaults internally
    GCACreateBackgroundMaskMovie(MD);
else % load any parameter changes from your personalized file
    cFunc = str2func([ip.Results.featureParameterMFile]);
    p = cFunc('step2load','backMask');
    GCACreateBackgroundMaskMovie(MD,p);
end
%% Perform the Background Subtraction
if isempty(ip.Results.featureParameterMFile)
    % if myParamsName empty just run the defaults internally
    GCABackgroundSubtractMovie(MD);
else % load any parameter changes from your personalized file
    cFunc = str2func([ip.Results.featureParameterMFile]);
    p = cFunc('step2load','backSubtract');
    GCABackgroundSubtractMovie(MD,p); 
end

%% Add filopodia intensity information
if isempty(ip.Results.featureParameterMFile)
    % if myParamsName empty just run the defaults internally
    GCAAddFilopodiaNormalizedIntensityMovie(MD);
else % load any parameter changes from your personalized file
    cFunc = str2func([ip.Results.featureParameterMFile]);
    p = cFunc('step2load','normFilo');
    GCAAddFilopodiaNormalizedIntensityMovie(MD,p);
end
%% Perform the Main Visualization

if isempty(ip.Results.featureParameterMFile)
    GCAVisualsMainOverlayMovie(MD);
else % load any parameter changes from your personalized file
    cFunc = str2func([ip.Results.featureParameterMFile]);
    p = cFunc('step2load','mainVis');
    GCAVisualsMainOverlayMovie(MD,p);
end

if isempty(ip.Results.featureParameterMFile)
    % Extract the default features for calc.
    GCAAnalysisExtractFilopodiaFeaturesMovie(MD);
else % load any parameter changes from your personalized file
    cFunc = str2func([ip.Results.featureParameterMFile]);
    p = cFunc('step2load','featExtract');
    if p.run; 
        GCAAnalysisExtractFilopodiaFeaturesMovie(MD,p);
    end
end

end

