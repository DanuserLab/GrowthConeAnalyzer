function GCARunVeilStemSingleMovie(MD,varargin)
% GCARunVeilStemSingleMovie

%% Input check
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('segParameterMFile',[]);


ip.parse(varargin{:});
%%
cFunc = str2func( ip.Results.segParameterMFile);
p = cFunc('step2load',0);
nonGCImage = p.NonGCImage; 
 if ~p.NonGCImage
    
    %% STEP I (Orange Block Figure S2)
    
    if isempty(ip.Results.segParameterMFile)
        % if myParamsName empty just run the defaults internally
        GCAgetNeuriteOrientMovie(MD);
    else % load any parameter changes from your personalized file
        cFunc = str2func( ip.Results.segParameterMFile);
        p = cFunc('step2load',1);
        GCAgetNeuriteOrientMovie(MD,p);
    end
    
    %% STEP II (Orange Block Figure S2)
    if isempty(ip.Results.segParameterMFile)
        GCAneuriteOrientConsistencyCheckMovie(MD);
    else
        cFunc = str2func( ip.Results.segParameterMFile);
        p = cFunc('step2load',2);
        GCAneuriteOrientConsistencyCheckMovie(MD,p);
    end
end
%% STEP III (Blue Block Figure S2)

if isempty(ip.Results.segParameterMFile)
    GCAReconstructVeilStemMovie(MD);
else
    cFunc = str2func( ip.Results.segParameterMFile);
    p = cFunc('step2load',3);
    if ~isempty(p.maskDirectory) 
        p.maskDirectory = [MD.outputDirectory_ filesep p.maskDirectory]; 
    end 
    GCAReconstructVeilStemMovie(MD,p);
end

%% Add the veilStem masks as an external mask process to the MD
% this way you can run through the old protrusion/windowing workflow and
% utilize movieViewer if you would like

 gcaMakeVeilStemMaskFolder(MD);

%% STEP IV (Green Block Figure S2)

if ~nonGCImage
    if isempty(ip.Results.segParameterMFile)
        GCAfindVeilStemLongestPathMovie(MD);
    else
        cFunc = str2func(ip.Results.segParameterMFile);
        p = cFunc('step2load',4);
        GCAfindVeilStemLongestPathMovie(MD,p);
    end
    
    %% Trajectory Partitioning
    
    if isempty(ip.Results.segParameterMFile)
        % if myParamsName empty just run the defaults internally
        GCAAnalysisPartitionTimeSeriesMovie(MD);
        GCAVisualsMakeNeuriteLengthOverlayMovie(MD,'colorByOutgrowthState',true);
    else % load any parameter changes from your personalized file
        cFunc = str2func( ip.Results.segParameterMFile);
        
        p = cFunc('step2load','partTraj');
        if p.labelNeuriteOutgrowthState
            GCAAnalysisPartitionTimeSeriesMovie(MD,'threshPause',p.threshPause,...
                'splineParam',p.splineParam,'OutputDirectory',p.OutputDirectoryP,'makePlots',p.makePlots,...
                'title',p.title);
        end % labelNeuriteOutgrowth
        
        if p.visualMovie
            GCAVisualsMakeNeuriteLengthOverlayMovie(MD,'colorByOutgrowthState',p.colorByOutgrowthState,...
                'colorVeil',p.colorVeil,'colorLength',p.colorLength,'Timer',p.Timer );
        end % visualMovie
        
    end % isempty
end
end % function