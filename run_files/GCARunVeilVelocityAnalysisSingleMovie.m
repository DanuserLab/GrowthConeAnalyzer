function [ output_args ] = GCARunVeilVelocityAnalysisSingleMovie(MD,varargin )
%GCARunVeilVelocityAnalysisSingleMovie
%% Input check
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('segParameterMFile',[]);

ip.parse(varargin{:});
%% Run the windowing
wind = WindowingProcess(MD);
MD.addProcess(wind);
windParams = MD.processes_{end}.funParams_;
if ~isempty(ip.Results.segParameterMFile)
    
    cFunc = str2func( ip.Results.segParameterMFile);
    c = cFunc('step2load',0); % test if gc or non gc image
    wP = cFunc('step2load','windParams');
    if isempty(wP.InputDirectory)
        inputDirectory = [MD.outputDirectory_ filesep 'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'III_veilStem_reconstruction']; 
    else 
        inputDirectory = wP.InputDirectory; 
    end 
    if ~c.NonGCImage % if growth cone use the neurite entrance point as the start point location for the windows.
        load([inputDirectory filesep 'Channel_' num2str(wP.ChannelIdx) filesep 'veilStem.mat']);
        SPIdx = veilStem(1).idxEnterNeurite;
        [SPy,SPx] = ind2sub(MD.imSize_,SPIdx);
        windParams.StartPoint = [SPx SPy];
    end
    windParams.OutputDirectory = [MD.outputDirectory_ filesep 'GCALocalVeilVelAnalysis' filesep 'windowing'];
    windParams.ChannelIdx = wP.ChannelIdx;
    windParams.PerpSize = wP.PerpSize;
    windParams.MinSize = wP.MinSize;
    windParams.MethodName = wP.MethodName;
    windParams.StartContour = wP.StartContour;
    windParams.StartPointPropag = wP.StartPointPropag;
    windParams.ParaSize = wP.ParaSize;
    windParams.ReInt = 61; 
    
    wP.Visual = wP.Visual;
    
else
    
    windParams.OutputDirectory = [MD.outputDirectory_ filesep 'GCALocalVeilVelAnalysis' filesep 'windowing'];
    windParams.ChannelIdx = 1;
    windParams.SegProcessIndex = 1;
    
    load([MD.outputDirectory_ filesep 'SegmentationPackage' filesep 'StepsToReconstruct' filesep ... 
        'III_veilStem_reconstruction' filesep 'Channel_1' filesep 'veilStem.mat']);
    SPIdx = veilStem(1).idxEnterNeurite;
    [SPy,SPx] = ind2sub(MD.imSize_,SPIdx);
    windParams.StartPoint = [SPx SPy];
    
    windParams.PerpSize = 1;
    windParams.MinSize = 500; % min size of the object to window
    windParams.MethodName = 'ConstantNumber';
    windParams.StartContour = 1;
    windParams.StartPointPropag = false;
    windParams.ParaSize = 5;
    windParams.ReInt= 61; 
    
    wP.Visual = true; 
end % ~isempty
parseProcessParams(MD.processes_{end},windParams);
MD.processes_{end}.run;
MD.save
%% Run protrusion sampling process
protSamp = ProtrusionSamplingProcess(MD);
MD.addProcess(protSamp);
protParams = MD.processes_{end}.funParams_;

protParams.OutputDirectory = [MD.outputDirectory_ filesep 'GCALocalVeilVelAnalysis' filesep 'protrusion_samples'];
parseProcessParams(MD.processes_{end},protParams);
MD.processes_{end}.run;

MD.save;
%% Edge velocity quantification
%         %% perform the edge velocity quantification
%         edgeOut = 'EdgeVelocityQuantification';
%         edgeVelocityQuantification(MD,'scale',true,'outLevel',6,'outputPath',edgeOut);
%% Make the Visual
if wP.Visual
    if isempty(wP.OutputDirectoryV)
        outDir = [ MD.outputDirectory_ filesep 'GCALocalVeilVelAnalysis' filesep 'Overlays' ];
    else
        outDir = wP.OuputDirectoryV;
    end
    GCAVisualsMakeLocalVeilDynamicsOverlayMovie(MD,'overlayLongPath',wP.overlayLongPath,.... 
        'OutputDirectory', outDir,'neuriteElongDir',wP.neuriteElongDir,'StartFrame', wP.StartFrame,...
        'EndFrame',wP.EndFrame,'MakeProtMap',wP.MakeProtMap,'SmoothActivityMap',wP.SmoothActivityMap);
end % wP.Visual

end