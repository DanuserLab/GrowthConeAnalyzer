function GCAPrepareDataPerformManualCropping(varargin)
% GCAPrepareDataPerformManualCropping
% Takes a file containing multiple stk (.tif) movies and prepares them for
% input into GCA.
% Also allows the user to manually crop each movie to specify region for
% segmentation.
% Currently only supports a single channel.
%
%
% Synopsis:    GCAPrepareDataPerformManualCropping
%           
%
% Input:
%      path : optional - path to the directory containing STK files input as a string.
%      If not input, the user will be asked to select a folder.
%
%      Optional parameter/value pairs
%           'OutputDirectory' : path where the output will be stored.
%           Default: isempty
%
%
%           'folderName'  :
%
% Francois Aguet, 09/01/2010
% Maria Bagonis  : added scrolling option for cropping,
%% Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addParameter('path', [], @(x) ischar(x) || isempty(x) || iscell(x));
ip.addParameter('OutputDirectory',[],@(x) ischar(x) || isempty(x) || iscell(x));


ip.addParameter('folderName',[]); % empty or character array (ie 'myFolderName')
% default empty : if empty will convert any underscore '_' in filename to a
% filesep '/'
% Ex. Control_20180807_01_mCherry.tif
% creates a Channel and GrowthConeAnalyzer folder in the new directory
% [path filesep 'Control' filesep '20180807'filesep '01']
% the value before the last '_' and .tif) the characters between the last
% '_' and the .tif will be used as the Channel Name
ip.addParameter('endNumberForID',true); % if true will use the last bit before the underscore for the ID. 
% if false will use a renumbered ID directory. 
ip.addParameter('extraFrame',true); % quick fix to add an extra frame for single images


ip.parse(varargin{:});
stkpath = ip.Results.path;
saveDir = ip.Results.OutputDirectory;


if isempty(stkpath)
    stkpath = uigetdir(pwd,'Select directory containing the STK files:');
    if (stkpath == 0)
        return;
    end
end


stkpath = [stkpath filesep];
stkList = [dir([stkpath '*.tif']) ; dir([stkpath '*.tiff']);  dir([stkpath '*.stk']); dir([stkpath '*.TIF'])];

N = length(stkList);
if N==0
    fprintf('No TIFF files found in input directory.\n');
end

if isempty(saveDir)
    %saveDir = upDirectory(stkpath,1); % make the directory above where the stks are stored
    saveDir = uigetdir(pwd,'Please select an Output Directory');
else
    saveDir = ip.Results.OutputDirectory;
end

if ~isdir(saveDir)
    mkdir(saveDir);
end
% if biosensor option on
% make folders with channels

tifNames = vertcat({stkList.name})';
movieNamesAll= cellfun(@(x) upDirectory(strrep(x,'_',filesep),1),tifNames,'uniformoutput',0);

movieNames = unique(movieNamesAll);
dataSet.movieNames = movieNames;

for iMovie = 1:length(movieNames)
    cMovie = movieNames{iMovie};
    
    
    idx = cellfun(@(x) strcmpi(x,cMovie),movieNamesAll);
    
    cTifs = tifNames(idx);
        
    [~,cNames] = cellfun(@(x) upDirectory(strrep(x,'_',filesep),1),cTifs,'uniformoutput',0);
  
    toCrop = char(cTifs{1}); % choose the first channel toCrop 
    
    fprintf('Converting: %s\n', cMovie);
 
    NChannels = numel(cNames);
    
    stackForCrop = stackRead([stkpath filesep toCrop]);
    
    [ny,nx,nFrames] = size(stackForCrop);
    dataSet.nFrames = nFrames;
    % initiate matrix and get other channels if they exist
    stackAllOrig = zeros(ny,nx,NChannels,nFrames);
    
        % get the other two channels
        for iCh = 1:NChannels
            stackC = stackRead([stkpath filesep cTifs{iCh}]);
            stackAllOrig(:,:,iCh,:) = stackC;
        end
 
    countRegions = 1;
  
        skip =0;
        reply2 = 'yes';
        
        imseriesshow(stackForCrop);
        if countRegions == 1
            replySkip = questdlg('Skip Movie?');
            if strcmpi(replySkip,'cancel');
                display('User Terminated Session: Please launch again if you would like to continue cropping');
                break
            end
        else
            replySkip = 'no';
        end
        
        %  hMsg = msgbox({'Use the expandable rectangle to crop the image'; 'Double Click when finished'},'help');
        
        while strcmpi(reply2,'yes');
            
            if strcmpi(replySkip,'yes')
                skip = 1;
                % copy the file over to a new internal directory 'skipped
                % Movies'
                skipDir = [pwd filesep 'Skipped Movies'];
                if ~isdir(skipDir)
                    mkdir(skipDir)
                end
                old = cellfun(@(x) [stkpath filesep x ],cTifs,'uniformoutput',0);
                new = cellfun(@(x) [skipDir filesep x],cTifs,'uniformoutput',0);
               
                arrayfun(@(x) movefile(old{x},new{x}),1:numel(new),'uniformoutput',0);
                               
                display('Movie Moved to a Skipped Movie Directory');
                reply2 = 'no';
                close all
            else % crop
                set(gcf,'Name','Crop the Image; Double Click Rectangle When Finished','NumberTitle','off')
                
                hr = imrect();
                
                pos = round(wait(hr));
                % cropRoi{k} = pos;
                stackCrop = stackForCrop(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3),:);
                %
                close all
                reply = 'Yes';
                while strcmpi(reply,'Yes')
                    imseriesshow(stackCrop)
                    hMsg=  msgbox({'Click through your movie to make sure cropping is correct';'When Finished Click OK'},'help');
                    uiwait(hMsg);
                    reply = questdlg('Re-do the cropping?');
                    if strcmpi(reply,'Yes')
                        close gcf
                        imseriesshow(stackForCrop);
                        hr = imrect();
                        pos = round(wait(hr));
                        
                        stackCrop = stackForCrop(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3),:);
                        
                    end % if strcmpi
                end
                % finish cropping all channels
                x = pos(1);
                x2 = pos(1)+pos(3);
                y = pos(2);
                y2= pos(2)+pos(4);
                deltX = abs(x2-x+1);
                deltY = abs(y2-y+1);
                stackAllCrop = zeros(deltY,deltX,NChannels,nFrames); % NOTE make channel number fliz
                % for each channel put into the multiChannelStack
                % stackAll = y,x,channel,frame
                % Crop All Channels
                for iChannel = 1:NChannels
                    stackAllCrop(:,:,iChannel,:) = stackAllOrig(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3),iChannel,:);
                end
                close gcf
            end % if skip
            if skip ~=1
                % make directories
                example = strrep(cTifs{1},'_',filesep);
                if countRegions >1
                    add = ['_' num2str(countRegions)];
                else
                    add = '';
                end
                if isempty(ip.Results.folderName)
                    
                    newDir = [saveDir filesep upDirectory(example,1) add ];
                else
                    if ip.Results.endNumberForID
                        [~,~,ID]  =  upDirectory(example,2);
                        newDir = [saveDir filesep ip.Results.folderName filesep ID add];
                    else
                        newDir = [saveDir filesep ip.Results.folderName filesep num2str(iMovie,'%03d') add];
                        
                    end
                end
                
                if ~isdir(newDir)
                    mkdir(newDir);
                end
              
                    % record cropping coords
                    cropRoi{countRegions,1} = pos;
                    save([newDir filesep 'cropRegion.mat'],'cropRoi');
               
                
                for iChannel = 1:NChannels
                    cChanDir = [newDir filesep 'Channels' filesep 'C' num2str(iChannel) '_' cNames{iChannel}];
                    
                    if ~isdir(cChanDir)
                        
                        mkdir(cChanDir);
                    end
                    
                    stackAllCrop = uint16(stackAllCrop);
                    
                    for iFrame = 1:nFrames
                        imwrite(stackAllCrop(:,:,iChannel,iFrame),[cChanDir filesep 'C' num2str(iChannel) '_' num2str(iFrame,'%03d') '.tif'],'tif');
                        if iFrame == 1 && ip.Results.extraFrame && size(stackAllCrop,4)==1
                            % write a second frame
                            imwrite(stackAllCrop(:,:,iChannel,iFrame),[cChanDir filesep 'C' num2str(iChannel) '_' num2str(iFrame+1,'%03d') '.tif'],'tif');
                        end
                    end
                end % iChannel
                
                % move the original stack
                newStackDir = [newDir filesep 'OriginalStack'];
                if ~isdir(newStackDir)
                    mkdir(newStackDir);
                end
                
                old = cellfun(@(x) [stkpath filesep x ],cTifs,'uniformoutput',0);
                new = cellfun(@(x) [newStackDir filesep x],cTifs,'uniformoutput',0);
                arrayfun(@(x) copyfile(old{x},new{x}),1:numel(new),'uniformoutput',0);
                close all
                % replot imserries
                imseriesshow(stackForCrop);
                % plot the previously plotted regions
                
                for iR = 1:numel(cropRoi)
                    hold on
                    pos = cropRoi{iR};
                    % plot the crop region
                    x = pos(1);
                    x2 = pos(1)+pos(3);
                    y = pos(2);
                    y2= pos(2)+pos(4);
                    line([x,x],[y,y2],'color','r','Linewidth',1);
                    line([x,x2],[y,y],'color','r','Linewidth',2);
                    line([x,x2],[y2,y2],'color','r','Linewidth',2);
                    line([x2,x2],[y,y2],'color','r','Linewidth',2);
                end % for iR
                
                reply2 = questdlg('Crop Another Region?');
                countRegions = countRegions +1;
                if strcmpi(reply2,'no')
                    arrayfun(@(x) movefile(old{x},new{x}),1:numel(new),'uniformoutput',0);
                    clear cropRoi
                    countRegions =1;
                    close all
                end
                
            end % if skip ~=1
        end % while cropMore (reply2)
    
    close all
end % for iMovie
end % function