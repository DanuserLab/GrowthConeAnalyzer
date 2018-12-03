function [ output_args ] = gcaAddNotes(MD,varargin)
%gcaAddNotes: This is a simple utility function mainly for me as I need to
% parse through the data manually. It allows one to add a flag to exclude
% (which helps when making the groupLists) and make
% comments about the structure and morphology/segmentation. (Might want to
% eventually directly incorporate into movie data itself.
%%
ip = inputParser;
ip.CaseSensitive = false;

ip.addParameter('ChannelIndexVeil',1); 
ip.parse(varargin{:});
% s = {'LastStepChecked','Include_Currently'...
%     ,'Specific_Questions','Task','Morphology','Crop','VeilStem','Filopodia','ExcludeAtStep'};

%
% reply = listdlg('ListString',s);


% look for notes
if exist([MD.outputDirectory_ filesep 'notes.mat'],'file')~=0;
    load([MD.outputDirectory_ filesep 'notes.mat']);
else
    % initiate the structure
    notes.Include = true; 
end



% toDocument = s(reply);
% put last step checked first
% put Last Step Checked First


% for i = 1:numel(toDocument)
%     % if field isn't empty show the previous comments
%
%     if strcmpi(toDocument{i},'LastStepChecked');

checkStepFlag  = questdlg('CheckAStep?');


stepsFolder{1} = [MD.outputDirectory_ filesep 'SegmentationPackage' filesep 'StepsToReconstruct' filesep ... 
    'I_neurite_orientation' filesep 'Channel_' num2str(ip.Results.ChannelIndexVeil) filesep '05_CandSeeds']; 
stepsFolder{2} = [MD.outputDirectory_ filesep 'SegmentationPackage/StepsToReconstruct'...
    filesep 'II_neurite_orientation_refinements' filesep 'Channel_' num2str(ip.Results.ChannelIndexVeil) ];

stepsFolder{3} =  [MD.outputDirectory_ filesep 'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'III_veilStem_reconstruction' filesep ...
    'Channel_' num2str(ip.Results.ChannelIndexVeil) filesep 'finalMask' filesep 'Overlays'];

stepsFolder{4} =  [MD.outputDirectory_ filesep 'GCAMainVisualization' filesep 'filoLength' filesep 'ForMainMovie_Feature_Movie'  ...
    filesep 'Channel1Detect_OverlaidOnChannel1']; 

%stepsFolder{4} =  [MD.outputDirectory_ filesep 'VisualizationOverlays/WholeNeurite/VeilWindows_SignalProcCheck/Channel_1'];
% stepsFolder{4} = [MD.outputDirectory_ filesep 'protrusion_samples_ProtrusionBased_windSize_5ReInit61' filesep 'EdgeVelocityQuantification' ...
%filesep 'VeilTrackCheck_ProtrusionBased_WindSize5' filesep 'Channel_1' ];
 
if (exist(stepsFolder{3})==0 || exist(stepsFolder{4})==0)
    notes.Include = false; 
else 

% test if exist
idxExist = cellfun(@(x) exist(x,'file'),stepsFolder);
idxExist = idxExist>0;
if idxExist(2) == 1
    load([stepsFolder{2} filesep 'backboneInfoFix.mat']);
    if length(backboneInfo)<MD.nFrames_-1
        idxExist(2) = 0 ;
    end
    
end
stepNames{1} = 'StepI';
stepNames{2}= 'StepII';
stepNames{3} = 'StepIII';
stepNames{4} = 'Final';

while strcmpi(checkStepFlag,'yes');
    
    toCheck = stepNames(idxExist);
    replyCheck = listdlg('ListString',toCheck);
    
    
    switch replyCheck
        case 1
            %                 cd([MD.outputDirectory_ filesep 'SegmentationPackage/StepsToReconstruct' ...
            %                     filesep 'I_neurite_orientation' filesep 'Channel_1' filesep 'CandSeeds']);
            %
            test = searchFiles('.png',[],stepsFolder{1},0,'all',1);
            for iFrame = 1:size(test,1);
                forMovie(:,:,:,iFrame) = imread(test{iFrame});
            end
            movieC = immovie(forMovie);
            imSize = MD.imSize_;
            h = setFigure(imSize(2),imSize(1),'on');
            
            imshow(forMovie(:,:,:,1),[]);
            implay(movieC);
            uiwait(h);
            
            
            
            
        case 2
            
            direct =   [MD.outputDirectory_ filesep 'SegmentationPackage/StepsToReconstruct'...
                filesep 'II_neurite_orientation_refinements' filesep 'Channel_1' filesep 'AfterFix'];
            if exist(direct)~=0
                
                
                
                test2 = searchFiles('.tif',[],[MD.outputDirectory_ filesep 'SegmentationPackage/StepsToReconstruct'...
                    filesep 'II_neurite_orientation_refinements' filesep 'Channel_1' filesep 'AfterFix'] ...
                    ,0,'all',1);
                for iFrame = 1:size(test2,1);
                    forMovie2(:,:,:,iFrame) = imread(test2{iFrame});
                end
                movieC = immovie(forMovie2);
                imSize = MD.imSize_;
                h = setFigure(imSize(2),imSize(1),'on');
                
                imshow(forMovie2(:,:,:,1),[]);
                
                implay(movieC);
                uiwait(h);
            else
                msg=   msgbox('No Orientation Outliers Found');
                uiwait(msg);
            end
        case 3
            direct =   stepsFolder{3};
            if exist(direct)~=0
                
                
                
                test2 = searchFiles('.png',[],[stepsFolder{3}] ...
                    ,0,'all',1);
                for iFrame = 1:size(test2,1);
                    forMovie2(:,:,:,iFrame) = imread(test2{iFrame});
                end
                
                movieC = immovie(forMovie2);
                imSize = MD.imSize_;
                h1 = setFigure(imSize(2),imSize(1),'on');
                
                imshow(forMovie2(:,:,:,1),[]);
                
                h = implay(movieC);
               
                
                 uiwait(h1);
                 if notes.Include;
                     x = 'Included'; 
                 else 
                     x = 'Excluded'; 
                 end 
 
                 reply = questdlg(['This movie is currently ' x ': Exclude Movie?']); 
                 
                 switch reply
                     case 'No'
                     notes.Include = true; 
                     case 'Yes' 
                         notes.Include = false; 
                     case 'Cancel'
                         notes.Include = notes.Include;  
                 end 
                 close(h)
            else 
                notes.Include = false; 
            end
            
        case 4
            
            direct= stepsFolder{4};
      
            
            if exist(direct)~=0
                
                
                
                test2 = searchFiles('.png',[],[stepsFolder{4}] ...
                    ,0,'all',1);
                for iFrame = 1:size(test2,1);
                    forMovie2(:,:,:,iFrame) = imread(test2{iFrame});
                end
                
                movieC = immovie(forMovie2);
                                imSize = MD.imSize_;
                                h1 = setFigure(imSize(2),imSize(1),'on');
                
                                imshow(forMovie2(:,:,:,1),[]);
                
                h = implay(movieC);
                
                
                                 uiwait(h1);
                if notes.Include;
                    x = 'Included';
                else
                    x = 'Excluded';
                end
                
                reply = questdlg(['This movie is currently ' x ': Exclude Movie?']);
                
                switch reply
                    case 'No'
                        notes.Include = true;
                    case 'Yes'
                        notes.Include = false;
                    case 'Cancel'
                        notes.Include = notes.Include;
                end
                close(h)
                
                
               
            end
    end
    checkStepFlag  = questdlg('Do you want to check this movie again?');
    
end % while
end % quick fix 
save([MD.outputDirectory_ filesep 'notes.mat'],'notes');
 
end

