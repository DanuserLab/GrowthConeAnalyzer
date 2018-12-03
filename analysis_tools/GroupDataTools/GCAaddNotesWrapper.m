function [ output_args ] = GCAaddNotesMovieList( toPlot,varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ip = inputParser;
ip.CaseSensitive = false;
ip.addParameter('ChannelIndexVeil',1);

ip.parse(varargin{:});

%%
projList = vertcat(toPlot.info.projList{:}); 
projList = projList(:,1); 


for iProj = 1:numel(projList)
    cd([projList{iProj,1} filesep 'GrowthConeAnalyzer']);
    load([projList{iProj} filesep 'GrowthConeAnalyzer' filesep 'movieData.mat']); 
    
    gcaAddNotes(MD,'ChannelIndexVeil',ip.Results.ChannelIndexVeil);
    close all
clear MD
end

