function [p] = GCALoadParameters_Defaults(varargin)
% GCALoadParameters_Defaults
%
% This function serves as a template and simply allows the user to more 
% easily modify/organize GCA parameter segmentation changes if necessary.
% The user is encouraged simply to modify this code and save as 
% myParameterChanges.m where myParameterChanges can be any string the user chooses. 
% In the GCARunFullPipeline.m function the user has the option to load manually or call 
% directly (for batch processing) this modified template function by its given name. 
% (If you call the .m file make sure the file is in the matlab path and has a unique .m name,
%  if you load it manually the path containing the .m file will be set automatically) 
% 
% Note:All GCA functions also have the defaults set internally and you 
% can reset any single parameter by name any time you call that function.
% (All GCA functions use Matlab's inputparser: see inputparser help page for details). 
% These loadParameter.m files are for organizational purposes. 
%% Input check
ip = inputParser;

ip.CaseSensitive = false;
ip.addParameter('step2load',1);
ip.parse(varargin{:});

%% Start
switch ip.Results.step2load
    case 'movie'
        %% Movie parameters
        p.pixelSize_ = 215;
        p.timeInterval_ = 5;
        %% Segmentation Parameters
    case 0
        p.NonGCImage = false;
        %%  Step I GCAGetNeuriteOrientMovie
    case 1
        % Input/Output
        p.OutputDirectory = []; % (Empty) or full pathname to save output
        % If empty: uses default OutputDirectory:
        % [movieData.outputDirectory_ filesep...
        %'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'I_neurite_orientation'];
        
        p.ChannelIndex= 1; % Channel Index (Note [filesep Channel_X] will always added to the outputDirectory
        % name.
        
        p.StartFrame = 'auto';
        p.EndFrame = 'auto';
        
        % Troubleshooting overlays
        p.TSOverlays = true;
        
        %%% CURRENTLY IN PIXELS %%%: Needs to be modified based on pixel size
        p.BBScale = [5 6 7 8 9 10];
        p.MaxRadiusLargeScaleLink = 10;
        p.MaxDistNoGeoTerm = 3;
        
        p.GeoThresh=0; %in cos(thetaRL)
        
        p.ThreshNMSResponse = 25;
        p.MinCCRidgeBeforeConnect=3;
        p.MinCCRidgeAfterConnect =5;
        
        p.filterBackEst = true;
        p.dilateLocalRegion=false;
        p.LRDilRad = 10;
        
        p.MinCCEntranceRidgeFirstTry=10;
        p.MaxDistBorderFirstTry=10;
        %% Step II GCAneuriteOrientConsistencyCheckMovie
    case 2
        % Input/Output
        p.OutputDirectory = [];
        p.InputDirectory = []; % (empty) or pathname where backboneInfo.mat is stored
        % If Empty, flags default InputDirectory:
        % [movieData.outputDirectory_ filesep ...
        % 'SegmentationPackage' filesep 'StepsToReconstruct' filesep ...
        % 'I_neurite_orientation']
        
        p.ChannelIndex=1;
        
        % In pixels
        p.SizeOfConsistencyRestraint=5;
        
        % Troubleshooting
        p.TSOverlays = true;
        p.CheckOrient = false;
        %%   Step III GCAReconstructVeilStemMovie
    case 3
        % Input/Output
        p.OutputDirectory=[];
        p.InputDirectoryStep1 =[];
        p.InputDirectoryStep2 = [];
        p.useCheckOrient= true;
        p.NonGCImage = false;
        
        p.ChannelIndex=1;
        
        % Restart Options
        p.StartFrame = 'auto';
        p.EndFrame = 'auto';
        
        % Troubleshoot Overlays
        p.TSOverlays = true;
        p.TSMovie = false;
        
        p.threshType = 'local'; %(local), global, or external
        p.threshMethod = 'otsu';
        p.GaussFilterSigma = 1;
        % Detection of amorphous veil pieces : Initial
        p.maskDirectory=[]; % directory with another input mask
        % if empty will use local thresholding
        p.LocalThresholdPatchSize= 75;
        
        % Detection of amorphous veil pieces : Morphological Opening
        p.DiskSizeLarge= 6; %in pixels:
        p.DiskSizeSmall= 3; %in pixels:
        % see gcaMorphologicalOpeningWithGeometryConstraints.m
        
        % Stem Radius Estimate
        p.reDilationRadius =4 ; % in pixels : see gcaResolveVeilStemCycles
        
        % Extend Path Search Via Bridging
        p.MaxRadiusBridgeRidges =5;
        
        %%
    case 4
        p.OutputDirectory =[] ;
        p.ChannelIndex = 1;
        
        p.StartFrame = 'auto';
        p.EndFrame = 'auto';
        
        p.TSOverlays = true;
        
        p.neuriteElongTS_medFiltWindSize=10; % see gcaFindOutliersFromMedFilt.m
        p.neuriteElongTSOutlier_outlierDef_k=9; % see gcaFindOutliersFromMedFilt.m
        
        p.extractVeilStemThickness= true;
        
    case 5
        %%  Step VI
    case 6
        p.OutputDirectory = [];
        
        p.InputDirectory = [];
        p.ChannelIndex = 1; % channelIndex for Filo
        p.ChannelIndexVeil=1;
        
        p.StartFrame = 'auto';
        p.EndFrame = 'auto';
        
        
        p.TSOverlaysReconstruct = false;
        p.TSOverlaysRidgeCleaning = true;
        
        
        % Option to detect Embedded actin bundles (for LifeAct and actin staining
        % only)
        p.detectEmbedded = true;
        
        p.rotateVeilStemNormals = true; % set to false if non-growth cone image
        
        p.FiloScale = 1.5;
        
        % Estimate background to localize region of interest
        p.filterBackEst = true; % flag to estimate high confidence
        % background using the image intensity histogram : < mean + 2std
        % is considered background
        p.dilateLocalRegion = 'test'; % (true,false, or 'test) flag to dilate further the
        % local region of interest estimation (ie decrease the background
        % estimation), if 'test' will check the percent background estimated and
        % only dilate if the estimated background is larger than that percentage
        
        p.LRDilRad = 10; % dilation radius of the structuring element
        % applied to inital guess of the localized region of interest.
        
        p.percentBackCutoff = 50; % if background is estimated to be
        % greater than this percentage of the image and dilateLocalRegion is set to
        % 'test' the background estimation will be dilated by LRDilRad (ie the
        % estimation is assumed to be too high, and too much of the image is
        % likely estimated as background
        
        % Cleaning Response of Steerable Filter
        p.multSTDNMSResponse =3;
        p.minCCRidgeOutsideVeil = 3;
        %ip.addParameter('filterVeilFromNMSRidgeHist',false); Add as option to
        %final version
        p.filterBasedOnVeilStemAttachedDistr=true;
        p.filterBasedOnFloatingVeilPieces=false; %
        
        % Linking Parameters Embedded
        p.geoThreshEmbedded = 0.5;
        p.maxRadiusLinkEmbedded=10;
        p.curvBreakCandEmbed=0.05;
        
        % Linking Parameters Candidate Building
        p.maxRadiusLink = 5; %
        p.geoThresh = 0.9;
        
        % Linking Parameters Traditional Filopodia/Branch Reconstruct
        p.maxRadiusConnectFiloBranch= 15;
        p.geoThreshFiloBranch = 0.5;
        %%
    case 7
        
        p.OutputDirectory = [];
        p.InputDirectory=[];
        p.ChannelIndex=1;
        
        p.StartFrame = 'auto';
        p.EndFrame ='auto';
        
        % Specific
        p.TSOverlays = false;
        
        p.InternalFiloOn = 3;
        
        p.PSFSigma  = 0.43 ; %% NOTE CHANGE THIS TO BE READ IN FROM MD.
        p.fitLengthInPix = 10;
        %% Feature Parameters
    case 'backMask'
        %% Background Mask Creation
        p.OutputDirectory=[];
        p.InputDirectoryFiloBranch=[];
        p.InputDirectoryVeilStem=[];
        p.ChannelIndex= 1;  % vector of channel indexes
        
        p.useFiloDetectBackEst = true; % load and use the same dilation parameters
        % that were originally used for the filopodia reconstruction.
        
        p.TSOverlays= true;
        p.maskToWrite= 'backMask';
    case 'backSubtract'
        %% Background subtraction
        p.ChannelIndex = 1;
        p.OutputDirectory=[];
        p.BackMaskDirectory=[];
        %% Normalization of filpodia intensities
        % PARAMETERS
    case 'normFilo'
        p.InputDirectory =[];
        p.InputDirectoryVeilStem=[];
        p.ChannelIndex =1;
        p.ChannelIndexVeil=1;
        
        p.frames=[];
        
        p.saveNormFactor=false; % flag to save a list of the norm
        % factors per frame.
        p.FeatureDirectory=[];
        
        p.UseBackSubtractImages=true;
    case 'mainVis'
        %%Input check
        
        p.InputDirectoryFilo = [];
        p.InputDirectoryVeil=[];
        
        p.OutputDirectory=[];
        
        p.ChannelIndexFilo = 1;
        p.ChannelIndexVeil= 1;
        
        p.Rewrite=false;
        
        p.AnalysisType= 'MainMovie'; %  Flag to load a predifined analinput.mat designed for visualization
        % MainMovie or MainMovieNoEmbed
        
        p.cMapLimits=[0,10]; % in microns
        
        p.createMask=true; % flag to create a mask of the full segmentation
        
        p.ScaleBar=true;
        
        p.writeCSV=true;
        
        p.frames =[];
        
        p.NonGCImage=false;
        
        % Use below information if p.AnalysisType is set to empty
        % if not it will load the defaults for the specified analysis
        p.filoTypes=[0 Inf];
        
        p.filterByBundleLength=[0.3,inf]; % in microns [minValue, maxValue]
        p.saveFiloByLengthAndSig=[];
        
        p.filterByFit=true;
        p.embedFitCriteria= 95;
        p.filoFitCriteria= 95;
        
        p.filterBasedOnGroupUpstream=0;
        p.filterBasedOnBranchConnection=0;
        p.filterIntNoConnect=true;
        
    case 'featExtract'
        p.run = true;
        
        p.InputDirectory=[];
        p.OutputDirectory=[];
        
        
        p.ChannelIndex =1;
        
        p.Rewrite =false;
        
        p.AnalysisType='defaultAnalysis';
        
        p.writeCSV  = true;
        
        % if use AnalysisType = 'defaultAnalysis' below ignored (use the
        % predefined filters
        p.filoTypes=[0 Inf];
        
        p.filterByBundleLength=[0.3,inf]; % in microns [minValue, maxValue]
        p.saveFiloByLengthAndSig=[];
        
        p.filterByFit=true;
        p.embedFitCriteria= 95;
        p.filoFitCriteria= 95;
        
        p.filterBasedOnGroupUpstream=0;
        p.filterBasedOnBranchConnection=0;
        p.filterIntNoConnect=true;
        
      case 'partTraj'
        
        p.labelNeuriteOutgrowthState = false;
        
            p.threshPause = 0.5 ;
            p.splineParam= 0.01;
            p.OutputDirectoryP = [];
            p.makePlots = true;
            p.title = [];
 
        p.visualMovie = true;
            p.colorByOutgrowthState = true;
            p.colorVeil = 'k';
            p.colorLength = 'k';
            p.Timer = false;
            p.ScaleBar = false;
            p.OutputDirectoryV = [];
            
    case 'windParams'
        p.InputDirectory =[];
        p.ChannelIdx = 1;
        p.SegProcessIndex = 1;
        
        p.PerpSize = 1;
        p.MinSize = 500; % min size of the object to window
        p.MethodName = 'ConstantNumber';
        p.StartContour = 1;
        p.StartPointPropag = false;
        p.ParaSize = 5;
        p.ReInt= 61;
        
        p.Visual = true;
            p.OutputDirectoryV = []; 
            p.overlayLongPath = true; 
            p.neuriteElongDir = []; 
            p.StartFrame = 1; 
            p.EndFrame = []; 
            p.MakeProtMap = true; 
            p.SmoothActivityMap =true;    
        
        
        
end % switch
