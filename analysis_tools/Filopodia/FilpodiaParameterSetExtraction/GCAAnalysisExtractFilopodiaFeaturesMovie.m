
function [ output_args ] = GCAAnalysisExtractFilopodiaFeaturesMovie(movieData,varargin)
%GCAAnalysisExtractFilopodiaFeaturesMovie
%   This function makes the default filopodia filter sets and extracts all
%   default filopodia features
%% INPUT
%
%   movieData - The MovieData object describing the movie, as created using
%   movieSelectorGUI.m
%
%   Parameter Structure Field Names:
%
% Generic Fields: (Input/Output Fields Needed for Wrapper)
%       ('OutputDirectory' -> Optional. A character
%       string specifying the directory to save the filopodia measurements to.
%       If not input, the output will be saved in the same directory
%       as the movieData, in a sub-directory called
%       'GrowthConeAnalyzer' filesep 'GCAMeasurementExtraction' filesep
%       'WholeNeurite'
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) on which to perform the filopodia
%       reconstruct
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
%% Check input

% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS

ip.addParameter('InputDirectory',[]);
ip.addParameter('OutputDirectory',[]);

ip.addParameter('ChannelIndex',1);
ip.addParameter('ProcessIndex',0);
ip.addParameter('Rewrite',false);
ip.KeepUnmatched = true;

ip.addParameter('AnalysisType','defaultAnalysis'); %  MainMovie, MainMovieNoEmbed, defaultAnalysis, if empty create default using 



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

% ip.addParameter('MainMovie',false); % flag to make the output for the
% % primary visualizations (all ext filo color coded by length);
% ip.addParameter('MainMovieNoEmbed',false);

ip.addParameter('writeCSV',true); 

ip.addParameter('trunc',false); 
ip.parse(varargin{:});
p = ip.Results;
%% set up directories
if isempty(ip.Results.InputDirectory)
    inDir = [movieData.outputDirectory_ filesep 'SegmentationPackage' ...
        filesep 'StepsToReconstruct' filesep ...
        'VII_filopodiaBranch_fits' filesep 'Channel_' num2str(ip.Results.ChannelIndex)];
else
    inDir = ip.Results.InputDirectory ;
end

if isempty(ip.Results.OutputDirectory)
    
    outDir = [movieData.outputDirectory_ filesep ...
        'GCAFeatureExtraction' ];
else
    outDir = ip.Results.OutputDirectory;
end

if ~isdir(outDir);
    mkdir(outDir);
end
% Define the filopodia filter
pFilt.pixelSizeNm= movieData.pixelSize_;
pFilt.trunc = ip.Results.trunc; 

pFilt.filoTypes = ip.Results.filoTypes;
pFilt.filterByBundleLength = ip.Results.filterByBundleLength; % in microns [minValue, maxValue]
pFilt.saveFiloByLengthAndSig=ip.Results.saveFiloByLengthAndSig;

pFilt.filterByFit= ip.Results.filterByFit;
pFilt.embedFitCriteria= ip.Results.embedFitCriteria;
pFilt.filoFitCriteria=ip.Results.filoFitCriteria;

pFilt.filterBasedOnGroupUpstream= ip.Results.filterBasedOnGroupUpstream;
pFilt.filterBasedOnBranchConnection = ip.Results.filterBasedOnBranchConnection;
pFilt.filterIntNoConnect=ip.Results.filterIntNoConnect;  

%% Init:
nFrames = movieData.nFrames_;
nChan = numel(p.ChannelIndex);
channels = p.ChannelIndex;
%% Loop for each channel
for iCh = 1:nChan 
    filoBranchNew = []; 
    %% Load the segmenatation analysis
    display('Please Be Patient This File Takes a While to Load...');
    %load([movieData.outputDirectory_ filesep 'filopodia_fits' filesep 'Filopodia_Fits_Channel_' num2str(channels(iCh)) filesep 'analInfoTestSave.mat']);
    load([inDir filesep 'filoBranch.mat'])
    filoBranchC = filoBranch; 
    clear filoBranch
    
    display('Finished Loading');
    
    % Add Check to make sure everything run completely
    
    % Get the analysis instructions, if not predefined create default using
    % the filopodia filters above
    if isempty(ip.Results.AnalysisType)
        analInput(1).filterType = [];
        analInput(1).paramFunc{1} = 'filoLength';
        analInput(1).paramName{1} = 'UserDefinedFilter';
        
        x.filoPart = 'Ext_';
        x.outPercent = false;
        x.umPerPixel = movieData.pixelSize_/1000;
        analInput(1).paramInput{1} = x; 
        
    else % load the predefined analysis protocols
        analInput =  GCAloadAnalysisInstructionFile(movieData,ip.Results.AnalysisType);
    end
    
    if ~isempty(ip.Results.AnalysisType)
        if ~isempty(regexp(ip.Results.AnalysisType,'Main','ONCE'))
            
            % save the analInput
            save([outDir filesep 'AnalysisInput.mat'],'analInput');
        else
            newDir = [outDir filesep 'Descriptor' filesep 'Filopodia'];
            if ~isdir(newDir)
                mkdir(newDir);
            end
            
            save([newDir filesep 'AnalysisInput.mat'],'analInput');
        end
    else
        save([outDir filesep 'AnalysisInput.mat'],'analInput');
        
    end
    
    %% Wrap through for each analysis type
    for iAnalType = 1:length(analInput);
        if ~isempty(ip.Results.AnalysisType)
            if isempty(regexp(ip.Results.AnalysisType,'Main','ONCE'))
                newFiloDir = [outDir filesep 'Descriptor'...
                    filesep 'Filopodia' filesep analInput(iAnalType).filterType ];
            else
                newFiloDir = outDir ;
            end
        else
            newFiloDir = [outDir];
        end
        
        % check for the filter directory
        if exist([newFiloDir filesep 'filoFilterSet.mat'])
            display([newFiloDir filesep 'filoFilterSet.mat Found']);
            
            % overwrite filter only if user specifies
            if ip.Results.Rewrite
                display(['Overwriting ' newFiloDir]);
                rmdir(newFiloDir,'s')
                %system(['mkdir -p '  newFiloDir]);
                mkdir(newFiloDir); 
                pFilt.filterType = analInput(iAnalType).filterType;
                
                % get the filopodia filter for analInput
                [filoFilterSet,filterParams] = GCACreateFilopodiaFilterSet(filoBranchC,pFilt);
                display('Re-writing Filopodia Filter Set');
                
                % save the filter set used
                save([newFiloDir filesep 'filoFilterSet.mat'],'filoFilterSet','filterParams');
            else % load the file
                load([newFiloDir filesep 'filoFilterSet.mat']);
                display('Loading Previous Filopodia Filter Set');
                % eventually check the filter so that know that it is from the same analysis
            end % ip.Results.Rewrite
            
        else % make the new directory and make the filter
            pFilt.filterType = analInput(iAnalType).filterType;
%             system(['mkdir -p '  newFiloDir]);

            mkdir(newFiloDir); 
            [filoFilterSet,filterParams] = GCACreateFilopodiaFilterSet(filoBranchC,pFilt);
            
            filterParams.InputDirectory = inDir;
            save([newFiloDir filesep 'filoFilterSet.mat'],'filoFilterSet','filterParams');
        end
        
        % perform the various assocatied feature extractions by calling the
        % associated function
        nParams  = numel(analInput(iAnalType).paramFunc);
        
        for iParamExtract = 1:nParams
            run = true;
            % Make the output directory: Each extraction has its own folder for now so have a place to
            % save associated plots and movies- will collect in a later
            % step (the movies can take a while to run... )
            
            nameParam = ['meas_'  analInput(iAnalType).paramName{iParamExtract}];
            analOutputDir = [newFiloDir filesep  analInput(iAnalType).paramFunc{iParamExtract}];
            
            % check for the .measC file and overwrite if necessary
            if exist([analOutputDir filesep nameParam '.mat' ])
                display([analOutputDir filesep nameParam '.mat file found']);
                if ip.Results.Rewrite
                    display(['Overwriting ' analOutputDir filesep nameParam '.mat']);  %give the user some feedback
                else
                    run = false; % don't run
                    display(['Skipping: ' nameParam]);
                end
            else % need to make the directory
                if ~isdir(analOutputDir)
                    %system(['mkdir -p '  analOutputDir]);
                    mkdir(analOutputDir); 
                end
            end
            
            if run
                % check for the measC file
                
                % Call the extraction function : all of these extraction functions have names starting with
                % GCAAnalysisExtract_ '
                % some are very short however it useful to keep the
                % different extractions as separate functions for organization purposes:
                % It is arguably easier for the user to cleanly add their own customized feature function
                % The output just needs to be a cell array of iFrames holding the distribution of values
                % for that feature for the given frame.  A
                % GCAAnalysisExtract_ template file is included
                paramFuncC = str2func(['GCAAnalysisExtract_' analInput(iAnalType).paramFunc{iParamExtract}]);
                inputC =  analInput(iAnalType).paramInput{iParamExtract};
                display(['Extracting Feature: ' analInput(iAnalType).paramName{iParamExtract}]);
                if ~isempty(inputC)
                    [measC,filoBranchNew] =  paramFuncC(filoBranchC,filoFilterSet, inputC);
                    
                else
                    [measC,filoBranchNew] = paramFuncC(filoBranchC,filoFilterSet);
                end
                if (~isempty(filoBranchNew) && strcmpi(analInput(iAnalType).paramFunc{iParamExtract},'filoCurvature'))
                    
                    filoBranchC = filoBranchNew;
                    filoBranch = filoBranchC;
                    filoBranchNew = [];
                    display('Updating filoBranch with Curvature Information');
                    save([inDir filesep 'filoBranch.mat'],'filoBranch','-v7.3');
                    clear filoBranch
                end
                timeStamp = clock;
                
                if ip.Results.writeCSV
                    dataMat = reformatDataCell(measC);
                    csvwrite([analOutputDir filesep nameParam '.csv'],dataMat);
                end
                
                save([analOutputDir filesep nameParam '.mat' ],'measC','timeStamp','inDir');
            end  % if run
        end % iParamExtract % iterate over all params using the same filter set
    end % iAnal % iterate over all analysis groups with the same filter set
end  % for iCh
end % Function