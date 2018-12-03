function [ output_args ] = GCAVisualsMainOverlayMovie( MD, varargin)
% GCAVisualsMainOverlayMovie
% Make the primary segmentation overlays whereby each filopodia segment is colored by
% its length. 
% Just combines:  
% GCAAnalysisExtractFilopodiaFeaturesMovie
% GCAVisualsMakeFeatureOverlayMovie with the appropriate input to make
% easier for the user to incorporate into a pipeline. 

if nargin < 1 || ~isa(MD,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('InputDirectoryFilo',[]);
ip.addParameter('InputDirectoryVeil',[]); 

ip.addParameter('OutputDirectory',[]);

ip.addParameter('ChannelIndexFilo',1);
ip.addParameter('ChannelIndexVeil',1); 

ip.addParameter('Rewrite',false);

ip.addParameter('AnalysisType','MainMovie'); %  Flag to load a predifined analinput.mat designed for visualization
% MainMovie or MainMovieNoEmbed

ip.addParameter('cMapLimits',[0,10]); % in microns

ip.addParameter('createMask',true); % flag to create a mask of the full segmentation 

ip.addParameter('ScaleBar',true); 

ip.addParameter('NonGCImage',false); 
% ip.addParameter('MainMovie',false); % flag to make the output for the
% % primary visualizations (all ext filo color coded by length);
% ip.addParameter('MainMovieNoEmbed',false);

ip.addParameter('writeCSV',true); 

ip.addParameter('frames',[]);

% filter input : use if ip.Results.AnalysisType =[]
ip.addParameter('filoTypes',[0 Inf]); % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
% default [0-Inf] Plot all segments that pass fit criteria 

ip.addParameter('filterByBundleLength',[0.3,inf]); % in microns [minValue, maxValue]
ip.addParameter('saveFiloByLengthAndSig',[]);

ip.addParameter('filterByFit',true);
ip.addParameter('embedFitCriteria', 95);
ip.addParameter('filoFitCriteria', 95);

ip.addParameter('filterBasedOnGroupUpstream',0);
ip.addParameter('filterBasedOnBranchConnection',0);
ip.addParameter('filterIntNoConnect',true);



ip.parse(varargin{:});

%% 
if isempty(ip.Results.OutputDirectory)
    featDirVis = [MD.outputDirectory_ filesep 'GCAMainVisualization'];
else
    featDirVis = ip.Results.OutputDirectory;
end

pExt.pixelSizeNm= MD.pixelSize_;


pExt.filoTypes = ip.Results.filoTypes;
pExt.filterByBundleLength = ip.Results.filterByBundleLength; % in microns [minValue, maxValue]
pExt.saveFiloByLengthAndSig=ip.Results.saveFiloByLengthAndSig;

pExt.filterByFit= ip.Results.filterByFit;
pExt.embedFitCriteria= ip.Results.embedFitCriteria;
pExt.filoFitCriteria=ip.Results.filoFitCriteria;

pExt.filterBasedOnGroupUpstream= ip.Results.filterBasedOnGroupUpstream;
pExt.filterBasedOnBranchConnection =ip.Results.filterBasedOnBranchConnection;
pExt.filterIntNoConnect=ip.Results.filterIntNoConnect;  

pExt.InputDirectory= ip.Results.InputDirectoryFilo  ; 
pExt.OutputDirectory = featDirVis; 
pExt.AnalysisType = ip.Results.AnalysisType; 
pExt.Rewrite = ip.Results.Rewrite; 

pExt.VeilDirectory = ip.Results.InputDirectoryVeil; 
pExt.NonGCImage = ip.Results.NonGCImage; 
% Always two steps to the visualizations:
% First we always specify some analysisType
% which specifies one or more
% filopodia filters and the associated GCAAnalysisExtract functions.
% In the case of the Main Visualization: We want to extract all viable
% filopodia/branch pieces as well as the calculated lengths of these
% pieces (so we can color the segmentation by filopodia/branch length)

GCAAnalysisExtractFilopodiaFeaturesMovie(MD,pExt);
%GCAAnalysisExtractFilopodiaFeaturesMovie(MD,'OutputDirectory',featDirVis,'MainMovie',true);

% Once the filopodia information to make the visual is extracted: we call the function to make the visual and create a corresponding mask
GCAVisualsMakeFeatureOverlayMovie(MD,'cMapLimits',ip.Results.cMapLimits,'screen2png',true,'interactive',false,...
    'FeatureDirectory',featDirVis,'createMask',ip.Results.createMask,... 
    'ScaleBar',ip.Results.ScaleBar,'frames',ip.Results.frames,'NonGCImage',ip.Results.NonGCImage,...
    'ChannelIndexVeil',ip.Results.ChannelIndexVeil);

end

