function [ output_args ] = GCATroubleshootMakeMovieOfReconstructMovie(movieData,varargin)
% GCATroubleshootMakeMovieOfReconstructMovie


%% INPUTPARSER
% for now check movieData separately.
if nargin < 1 %|| ~isa(movieData,'MovieData')
    %error('The first input must be a valid MovieData object!')
    [name,path ]=  uigetfile(pwd,'Please Select a MovieData.mat File');
    load([path filesep name]); 
    movieData = MD; 
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('nonGCImage',[]);

ip.addParameter('OutputDirectory',[],@(x) ischar(x));
ip.addParameter('InputDirectoryFiloBranch',[],@(x) ischar(x));
ip.addParameter('InputDirectoryVeilStem',[],@(x) ischar(x));
ip.addParameter('ChannelIndexFiloBranch',1);
ip.addParameter('ChannelIndexVeilStem',1);


ip.addParameter('frames',1)
ip.addParameter('outDirType','perFrame');

%ip.addParameter('figFiles',true); % need to add
ip.addParameter('epsFiles',false);

ip.addParameter('makeMovies',false); % flag to make movies : Requires ffmpeg

ip.addParameter('writeTitles',true);
ip.addParameter('screen2png',false);
ip.parse(varargin{:});
%% Initiate



imDir  = movieData.getChannelPaths{1};
frames = ip.Results.frames;
%% Wrap

if isempty(ip.Results.OutputDirectory)
    outDir = [movieData.outputDirectory_ ];
else
    outDir = ip.Results.OutputDirectory;
end

if isempty(ip.Results.InputDirectoryFiloBranch)
    % PARAMETERS
    inDirFiloBranch = [movieData.outputDirectory_ filesep...
        'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'VII_filopodiaBranch_fits' ...
        filesep 'Channel_' num2str(ip.Results.ChannelIndexFiloBranch)];
else
    inDirFiloBranch = [ip.Results.FiloBranch filesep 'Channel_' num2str(ip.Results.ChannelIndexFiloBranch)];
end

if isempty(ip.Results.InputDirectoryVeilStem)
    %if ip.Results.nonGCImage
        inDirVeilStem = [movieData.outputDirectory_ filesep ...
            'SegmentationPackage' filesep 'StepsToReconstruct' ...
            filesep 'III_veilStem_reconstruction' filesep 'Channel_' num2str(ip.Results.ChannelIndexVeilStem)];
    %else
%         inDirVeilStem = [movieData.outputDirectory_ filesep ...
%             'SegmentationPackage' filesep 'StepsToReconstruct' ...
%             filesep 'IV_veilStem_length' filesep 'Channel_' num2str(ip.Results.ChannelIndexVeilStem)];
    %end
else
    inDirVeilStem = ip.Results.InputDirectoryVeilStem;
end



for iFrame = 1:length(frames)
    if strcmpi(ip.Results.outDirType,'perFrame');
        saveDir = [outDir filesep 'Reconstruct_Movies' filesep 'Frame_' num2str(frames(iFrame),'%03d')];
    else
        saveDir = [outDir filesep 'Reconstruct_Movie'];
    end
    
    if ~isdir(saveDir)
        mkdir(saveDir)
    end
    
    
    load( [inDirVeilStem filesep 'veilStem.mat']);
    
    
    load([inDirFiloBranch filesep 'filoBranch.mat']);
    
    pixSize_um = movieData.pixelSize_/1000;
    
    inputFrames = frames(iFrame);
    
    [hSet filoFilterSet,filterParams] = GCATroubleShootMakeMovieOfReconstruct(filoBranch,veilStem,inputFrames,pixSize_um,imDir,...
        'writeTitles',ip.Results.writeTitles);
    
    count = 1;
    for i = 1:length(hSet)
        % double it up to slow it down
        %     saveas(hSet(i).h, [saveDir filesep num2str(count,'%03d') '.png']);
        %     count= count+1;
        cFileName = [saveDir filesep num2str(count,'%03d') '.png'];
        
        if ip.Results.screen2png
            helperScreen2png(cFileName,'figureHandle',hSet(i).h);
        else
            saveas(hSet(i).h,cFileName);
        end
        if ip.Results.epsFiles
            saveas(hSet(i).h,[saveDir filesep num2str(count,'%03d') '.eps'],'psc2');
        end
        count = count+1;
        
        % save two more of the last frame because windows media player clips
        % the last two frames
        %
        %     if i == length(hSet)-2;
        %
        %         saveas(hSet(i).h,[saveDir filesep num2str(count,'%03d') '.png']);
        %         saveas(hSet(i).h,[saveDir filesep num2str(count,'%03d') '.png']);
        %         saveas(hSet(i).h,[saveDir filesep num2str(count,'%03d') '.fig']);
        %         saveas(hSet(i).h,[saveDir filesep num2str(count,'%03d') '.eps'],'psc2');
        %     end
        
    end
    
    load([inDirFiloBranch filesep 'params.mat']);
    copyfile([inDirFiloBranch  filesep 'params.mat'],[saveDir filesep 'paramsRec.mat']);
    
    if ip.Results.makeMovies
        %     execute = ['ffmpeg -r 1 -i ' saveDir filesep '%03d.png' ...
        %         ' -b 2000k ' saveDir filesep 'ReconstructMovie' num2str(frames(iFrame),'%03d') '.wmv'];
        %     system(execute);
        
        execute = ['ffmpeg -r 1 -i ' saveDir filesep '%03d.png' ...
            ' -b 2000k ' '-pix_fmt' ' yuv420p ' saveDir filesep 'ReconstructMovie' num2str(frames(iFrame),'%03d') '.mp4'];
        system(execute);
    end
end

end