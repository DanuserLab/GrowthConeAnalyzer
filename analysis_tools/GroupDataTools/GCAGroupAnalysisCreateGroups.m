function [gcaGroupData] = GCAGroupAnalysisCreateGroups(varargin)
%% GCAGroupAnalysisCreateGroups
% Creates an initial gcaGroupData.mat file defining a list of movie paths 
% grouped by condition ready for final analysis. 
%
%% Input Parser
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('Interactive',true); % 
ip.addParameter('InputDirectory',[]); % if empty and Interactive == true, will ask the user to input 
%  if empty and Interactive == false, will assume the current directory
ip.addParameter('OutputDirectory',[]); % where the gcaGroupData.mat file will be saved 

ip.addParameter('screenByNotes',false); % flag to screen by notes associated with the movie
ip.addParameter('manualProjSelect',false); % flag to manually select the projects
ip.addParameter('searchPhrase',[]); % used KD as my search phrase so I can keep other 
% folders in there and it will ignore them 

ip.parse(varargin{:});
%%
screenByNotes = ip.Results.screenByNotes;

if isempty(ip.Results.InputDirectory)
    if ip.Results.Interactive
        inDir= uigetdir(pwd,'Please select your grouped GCAAnalysis Directory');
    else
        inDir = pwd;
    end  
else 
    inDir = ip.Results.InputDirectory; 
end 

if isempty(ip.Results.OutputDirectory)
    if ip.Results.Interactive
        outDir= uigetdir(pwd,'Please select a directory to store your gcaGroupData.mat file');
    else
        outDir= pwd;
    end  
else 
    outDir = ip.Results.OutputDirectory; 
end


s = dir(inDir);
if ~isempty(ip.Results.searchPhrase)
idxKeep =  arrayfun(@(x) ~isempty(regexp(x.name,ip.Results.searchPhrase,'ONCE')),s);
else 
   idxKeep = true(size(s,1),1); 
   idxKeep(1:2) = 0; % take out the comments 
   idxKeep(vertcat(s(:).isdir)==0) = 0; 
end 
s = s(idxKeep);
groupName = arrayfun(@(x) [x.name],s,'uniformoutput',0);

idxInclude  = listSelectGUI(groupName,[],'move');

groupName = groupName(idxInclude) ;

nGroups = numel(groupName);
%% START
% for each group
for iGroup = 1:nGroups
    
    gcaGroupData.info.names{iGroup} = groupName{iGroup};
    % for each group look for experimental day folders
    directoriesExpDay = dir([inDir filesep groupName{iGroup}]);
    idxDir = vertcat(directoriesExpDay(:).isdir);
    directoriesExpDay = directoriesExpDay(idxDir);
    % keep only those experimental day folders with numbers
    
    idxKeep = arrayfun(@(x) ~isnan(str2double(strrep(x.name,'_',''))),directoriesExpDay);
    directoriesExpDay = directoriesExpDay(idxKeep);
       
    directoriesExpDay =   arrayfun(@(x) [inDir filesep groupName{iGroup} filesep x.name],directoriesExpDay,'uniformoutput',0);
    if exist([inDir filesep groupName{iGroup} filesep 'groupColor.mat'])
        load([inDir filesep groupName{iGroup} filesep 'groupColor.mat'])
    else
        if iGroup ==1 % assume the first group is the control
            
            [mainColor,shades] = gcaHelperGetColorForGroup(iGroup);
            groupColor.single = mainColor;
            groupColor.mult = shades;
        else
            [mainColor,shades] = gcaHelperGetColorForGroup(2);% for now just set all the rest to red.
            groupColor.single = mainColor;
            groupColor.mult = shades;
        end
    end
    % for each day folder look for movie folders (a flag will be in the
    % notes on whether not to currently include or not)
    for iDay = 1:length(directoriesExpDay);
        directoriesMovie = dir(directoriesExpDay{iDay});
        idxDir = vertcat(directoriesMovie(:).isdir);
        directoriesMovie = directoriesMovie(idxDir);
        % keep only the number folders in case I put crap in there.
        idxKeep = arrayfun(@(x) ~isnan(str2double(strrep(x.name,'_',''))),directoriesMovie);
        
        directoriesMovie = directoriesMovie(idxKeep);
        
        forProjList  = arrayfun(@(x) [directoriesExpDay{iDay} filesep  x.name],directoriesMovie,'uniformoutput',0);
        %%
        if ip.Results.manualProjSelect
            idxInclude = listSelectGUI(forProjList,[],'move');
        else
            if screenByNotes
                % test the notes to see if should include
                notesAll = cellfun(@(x) load([x filesep 'GrowthConeAnalyzer' filesep 'notes.mat']), forProjList);
                
                idxInclude = arrayfun(@(x) notesAll(x).notes.Include,1:length(notesAll));
            else
                idxInclude = true(size(forProjList,1),1); % include them all
            end % screenByNotes
        end
        %%
        forProjList = forProjList(idxInclude);
        namesProj = cell(numel(forProjList),1);
        
        for iProj = 1:numel(forProjList)
            
            [~,date] = upDirectory(forProjList{iProj},2,1);
            [~,num] = upDirectory(forProjList{iProj},1,1);
            [~,name] = upDirectory(forProjList{iProj},3,1);
            
            namesProj{iProj,1}= [name '_' date '_' num];
        end
        projListC = [forProjList namesProj];
        projListGroup{iDay} = projListC;
    end % for iDay
    
    gcaGroupData.info.projList{iGroup} = vertcat(projListGroup{:});
    nProjs = size(gcaGroupData.info.projList{iGroup},1);
    grouping{iGroup} = iGroup.*ones(nProjs,1);
    gcaGroupData.info.color{iGroup} = groupColor.single;
    gcaGroupData.info.colorShades{iGroup} = groupColor.mult;
    
    clear projListGroup
end % iGroup
gcaGroupData.info.grouping = vertcat(grouping{:});


save([outDir filesep 'gcaGroupData.mat'],'gcaGroupData');

end % function