function [ ID,group,date,numID ] = gcaGetNeuriteID(path,varargin)
% gcaGetNeuriteID. Simple function based on 

ip = inputParser; 
ip.CaseSensitive = false;

ip.addParameter('biosensor',false);
ip.parse(varargin{:});

[~,date] = upDirectory(path,3,1);
[~,numID] = upDirectory(path,2,1);
[~,group] = upDirectory(path,4,1);

if ip.Results.biosensor
    [~,biosensor] = upDirectory(path,5,1);
    ID = [biosensor '_' date '_' numID];
else
    ID = [group ' ' date ' ' numID];
end
end

