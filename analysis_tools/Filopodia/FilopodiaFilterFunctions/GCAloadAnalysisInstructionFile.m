function [ analInput] = GCAloadAnalysisInstructionFile(movieData,instructionType)
% GCAloadAnalysisInstructionFlie


% load some predefined analInput files. 
%
% Output: 
% analInput: N x 1 structure where N is the number of filopodia filter
% types you would like to run 
%           with fields 
%           .filterType : String specifying one of the pre-defined filter types
%           .paramFunc : an Mx1 cell array of strings specifying the name of the 
%                        'GCAAnalysisExtract_  .m' function to call
%                         to call GCAAnalysisExtract_filoLength.m : this
%                         .paramFunc{1} = 'filoLength'; 
%           .paramName : an Mx1 cell array of strings specifying the name
%                        of the feature where M is the number of features 
%                        to run using the same filopodia filter. 
%           .paramInput : an Mx1 cell array of strings specifying any input
%                        parameters to the GCAAnalysisExtract_ .m function 
%                        being called. (Some related features may call the 
%                        same function using different input) 
%              
%% Predefined Analysis Input


switch instructionType
   %% Start formatting analysis Input
    % Note need to eventually make this a switch...
    case 'MainMovie'
        analInput(1).filterType = 'Validation';
        analInput(1).paramFunc{1} = 'filoLength';
        analInput(1).paramName{1} = 'ForMainMovie';
        %x.filoPart = 'Ext_';
        x.filoPart = 'Tot';
        x.outPercent = false;
        %x.outPercent = true;
        %x.filterZeroPercent = false;
        x.umPerPixel = movieData.pixelSize_/1000;
        analInput(1).paramInput{1} = x;
        
    case 'MainMovieNoEmbed'
        analInput(1).filterType = 'Validation_NoEmbed';
        analInput(1).paramFunc{1} = 'filoLength';
        analInput(1).paramName{1} = 'ForMainMovieNoEmbed';
        x.filoPart = 'Ext_';
        %x.filoPart = 'Tot';
        x.outPercent = false;
        x.umPerPixel = movieData.pixelSize_/1000;
        analInput(1).paramInput{1} = x;   
    case 'defaultAnalysis'
        %% Features for Filter I 'ConnectToVeil_LengthInt' : Veil Measurements that require good fitting (Length and Intensity)Filopodia/Branch Length (0 and 1st Order) : 
        analInput(1).filterType =   'ConnectToVeil_LengthInt';
            %% Features that call GCAAnalysisExtract_filoLength.m
            analInput(1).paramFunc{1} = 'filoLength'; % % function ID
            analInput(1).paramName{1} = 'filoLengthToVeil'; % featureName
            % Function Input
            x.filoPart = 'Ext_';
            x.outPercent = false;
            analInput(1).paramInput{1} = x;
            clear x

            % Veil Embedded Actin Bundle (Life-Act Only) : Length
            analInput(1).paramFunc{2} = 'filoLength';
            analInput(1).paramName{2} = 'filoLengthEmbedded';
            % Function Input
            x.filoPart = 'Int_';
            x.outPercent = false;
            analInput(1).paramInput{2} = x;
            clear x
        
            % Full Actin Bundle (Life-Act Only) : Length
            analInput(1).paramFunc{3} = 'filoLength';
            analInput(1).paramName{3} = 'filoLengthFullActinBundle';
            % Function Input
            x.filoPart = 'Tot';
            x.outPercent = false;
            analInput(1).paramInput{3} = x ;
            clear x
            
            % Percentage of Each Actin Bundle Embedded
            analInput(1).paramFunc{7} = 'filoLength';
            analInput(1).paramName{7} = 'percentEachActinBundleEmbed';
            x.filoPart = 'Tot';
            x.outPercent = true;
            analInput(1).paramInput{7} = x;

            %%  Features that call GCAAnalysisExtract_filoAvgIntensity.m
            
            analInput(1).paramFunc{4} = 'filoAvgIntensity'; % function ID
            analInput(1).paramName{4} = 'filoIntensityToVeil_Norm'; % paramName for output
            x.filoPart = 'Ext_';
            x.normToVeil = true;
            analInput(1).paramInput{4} = x; % other information for function
            clear x

            % Veil Embedded Actin Bundle
            analInput(1).paramFunc{5} = 'filoAvgIntensity';
            analInput(1).paramName{5} = 'filoIntensityEmbedded_Norm';
            x.filoPart = 'Int_';
            x.normToVeil = true;
            analInput(1).paramInput{5} = x;
            clear x
            
            %%  Features that call GCAAnalysisExtract_percentOfActinBundlesVeilEmbedded.m

            % Percentage of total Actin Bundles Embedded
            analInput(1).paramFunc{6} = 'percentOfActinBundlesVeilEmbedded'; % % function ID
            analInput(1).paramName{6} = 'percentTotalActinBundlesVeilEmbedded'; % paramName-
            analInput(1).paramInput{6} = [];
            %% Features that call GCAAnalysisExtract_filoOrient.m
 
            % Orientation of Filopodia
            analInput(1).paramFunc{8} = 'filoOrient';
            analInput(1).paramName{8} = 'filoOrientation';
            analInput(1).paramInput{8} = [];
            
            %% Features that call GCAAnalysisExtract_filoCurvature.m    
            % Curvature of Filopodia
            analInput(1).paramFunc{9} = 'filoCurvature';
            analInput(1).paramName{9} = 'filoMaxCurvature';
            analInput(1).paramInput{9} = [];
        
        %% Filter II: 'ConnectToVeil_Density': Density Filo Veil: Measurements that do not require accurate fitting and are more sensitive to false negatives. 
        analInput(2).filterType = 'ConnectToVeil_Density'; % note moved orient up so can potentially correlate info

            %% Features that call GCAAnalysisExtract_filoDensityAlongVeil.m   
            analInput(2).paramFunc{1} = 'filoDensityAlongVeil';
            analInput(2).paramName{1} = 'filoDensityAlongVeil';

            load([movieData.outputDirectory_ filesep 'SegmentationPackage' filesep ...
                'StepsToReconstruct' filesep 'III_veilStem_reconstruction' filesep 'Channel_1'...
                filesep 'veilStem.mat']);
            analInput(2).paramInput{1} = veilStem;
        
        %% Filter III: 'Branch2ndOrder_LengthInt'
        % Length/Int 2nd Order Filo/Branches
        % Note my definitions are 0 filo attached to veil no branch
        % 1 filo attached to veil with a branch
        % 2 branch attached to filo type 1
        analInput(3).filterType = 'Branch2ndOrder_LengthInt';
            %%
            analInput(3).paramFunc{1} = 'filoOrient';
            analInput(3).paramName{1} = 'branchOrientation_2ndOrder';
            analInput(3).paramInput{1} = [];
            %%
            analInput(3).paramFunc{2} = 'filoLength';
            analInput(3).paramName{2} = 'branchLength_2ndOrder';
            x.filoPart = 'Ext_';
            x.outPercent = false;
            analInput(3).paramInput{2} = x;          
            clear x
            
            %%
            analInput(3).paramFunc{3} = 'filoCurvature';
            analInput(3).paramName{3} = 'branchMaxCurvature_2ndOrder';
            analInput(3).paramInput{3} = [];

            %%
            analInput(3).paramFunc{4} = 'filoAvgIntensity';
            analInput(3).paramName{4} = 'branchIntensity_2ndOrder';
            %analInput(3).paramInput{4} = 'Ext';
            x.filoPart = 'Ext_';
            x.normToVeil = true;
            analInput(3).paramInput{4} = x;
            clear x        
        %%  Filter IV: 'Branch2ndOrder_Density_NoZero'
        analInput(4).filterType = 'Branch2ndOrder_Density_NoZero';
            %%
            analInput(4).paramFunc{1} = 'distanceToBranch';
            analInput(4).paramName{1} = 'branchDistanceFromVeil'  ;
            analInput(4).paramInput{1} = [];
            %%
            analInput(4).paramFunc{2} =  'filoDensityAlongBranch';
            analInput(4).paramName{2} = 'branchDensity_2ndOrder';
            analInput(4).paramInput{2} = [];        
        %% Filter V: 'Validation' (for branch complexity calc)
        analInput(5).filterType = 'Validation';
            %%
            analInput(5).paramFunc{1} = 'filoBranchComplexity';
            analInput(5).paramName{1} = 'filoBranchComplexity';
            analInput(5).paramInput{1} = [];      
end % switch instruction type