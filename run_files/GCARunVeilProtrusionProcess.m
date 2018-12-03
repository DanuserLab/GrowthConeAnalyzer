function [ output_args ] = GCARunVeilProtrusionProcess(MD)
% GCARunVeilProtrusionProcess
  %% Protrusion Vector Calc (For Veil/Stem Dynamics and Filopodia Orientation Calcs)
    prot = ProtrusionProcess(MD);
    MD.addProcess(prot);
    protParams = MD.processes_{end}.funParams_;
    
    % Check to make external process run 
    idxExtSeg = find(cellfun(@(x) sum(strcmpi(x.name_,'External Segmentation')),MD.processes_));
    if length(idxExtSeg) >1 
        dis('More than one External Segmentation Process found: Using most recent'); 
        idxExtSeg = idxExtSeg(1); 
    end 
   
    
    protParams.ChannelIndex = 1;
    protParams.SegProcessIndex = idxExtSeg;
    protParams.OutputDirectory = [MD.outputDirectory_ filesep 'GCALocalVeilVelAnalysis' filesep ...
        'protrusion']; 
    if ~isdir(protParams.OutputDirectory)
        mkdir(protParams.OutputDirectory)
    end 
    
    protParams.BatchMode = 1; 
    parseProcessParams(MD.processes_{end},protParams);
    MD.processes_{end}.run
    MD.save;


end

