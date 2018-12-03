function [ output_args ] = gcaMakeVeilStemMaskFolder(movieData,varargin)
%gcaMakeVeilStemMaskFolder

%% Check input
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end

ip = inputParser;

ip.CaseSensitive = false;


ip.addParameter('OutputDirectory',[]);
ip.addParameter('InputDirectory', []);

ip.addParameter('addExternalSegProcess',true);
ip.addParameter('rewriteExtProcess',true); % flag to replace any external segmentation
% processes : if true will delete any other external segmentation processes

ip.addParameter('OutputDirectoryMD',[],@(x) ischar(x)); % Directory where to save the external seg process masks 

ip.addParameter('ChannelIndex',1); % currently only supports one channel 

ip.addParameter('Overlays',true);
ip.addParameter('saveEPS',true); 
ip.parse(varargin{:});
%% Initialize

if isempty(ip.Results.InputDirectory) 
    inputDirectory = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep 'StepsToReconstruct' filesep ...
     'III_veilStem_reconstruction' filesep 'Channel_' num2str(ip.Results.ChannelIndex)]; 
else 
    inputDirectory = ip.Results.InputDirectory; 
end 

if isempty(ip.Results.OutputDirectory) 
    outputDirectory= [inputDirectory filesep 'finalMask']; 
else 
    outputDirectory = ip.Results.OutputDirectory; 
end 

if ~isdir(outputDirectory)
    mkdir(outputDirectory); 
end 

if isempty(ip.Results.OutputDirectoryMD)
    outDirMD = [movieData.outputDirectory_ filesep...
        'GCALocalVeilVelAnalysis' filesep 'masks_VeilStemFinal'];
else
    outDirMD = ip.Results.OutputDirectoryMD;
end

if ~isdir(outDirMD)
    mkdir(outDirMD)
end 

fileToLoad = [inputDirectory filesep 'veilStem.mat'];
load(fileToLoad);
nImTot = numel(veilStem);
fmt = ['%0' num2str(ceil(log10(nImTot))) 'd'];
%%
for iFrame = 1:nImTot
    mask  = veilStem(iFrame).finalMask;
    
    % should have already been done but just to make sure
    % get largest CC
    mask = double(getLargestCC(mask));
    
    % set the border to zero - this just helps make sure that the windows do not
    % error
    dims = size(mask);
    mask(1:dims(1),1) =0;
    mask(1:dims(1),dims(2))=0;
    mask(1,1:dims(2))= 0;
    mask(dims(1),1:dims(2)) =0;
    
    % should have been filled but again just in case
    
    mask = imfill(mask,'holes');
    mask = double(getLargestCC(mask)); % just in case break after remove edge pixels
    
    imwrite(mask,[outputDirectory filesep 'veilStemMask' num2str(iFrame,fmt) '.tif']);
    if ip.Results.Overlays
        
        img =  movieData.channels_(1).loadImage(iFrame);
        img = double(img);
        [ny,nx] = size(img);
        setFigure(nx,ny,'off');
        imshow(-img,[]);
        hold on
        roiYX = bwboundaries(mask);
        text(5,10,'VeilStem Outline', 'Color','g');
        cellfun(@(x) plot(x(:,2),x(:,1),'color','g'),roiYX);
        if ~isdir([outputDirectory filesep 'Overlays']);
            mkdir([outputDirectory filesep 'Overlays']) ;
        end
        saveas(gcf,[outputDirectory filesep 'Overlays' filesep 'VeilStemOverlay' num2str(iFrame,fmt) '.png']);
        if ip.Results.saveEPS
               saveas(gcf,[outputDirectory filesep 'Overlays' filesep 'VeilStemOverlay' num2str(iFrame,fmt) '.eps'],'psc2');               
        end 
        close gcf
    end
    clear mask
end
% add to process via external segmentation
if ip.Results.addExternalSegProcess
    
    if ip.Results.rewriteExtProcess
        idxExt = find(cellfun(@(x) sum(strcmpi(x.name_,'External Segmentation')),movieData.processes_));
        
        while ~isempty(idxExt)
            
            movieData.deleteProcess(idxExt(1));
            movieData.save;
            idxExt = find(cellfun(@(x) sum(strcmpi(x.name_,'External Segmentation')),movieData.processes_));
        end
        
    end % if rewriteExtProcess

    
    
    extProc = ExternalSegmentationProcess(movieData);
    % % % Save the process in the movie object
    movieData.addProcess(extProc);
    params = movieData.processes_{end}.funParams_;
    params.OutputDirectory = outDirMD;
    params.InputData{1} = outputDirectory; 
    parseProcessParams(movieData.processes_{end},params);
    movieData.processes_{end}.run;
    movieData.save;   
end
end  %function
