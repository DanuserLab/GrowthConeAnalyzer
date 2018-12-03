function [ filterParams ] = GCAloadDefaultFilopodiaFilters( filterType)
% GCAloadFilopodiaFilterParameters

%% Current predefined filopodia filters: 
% 'ConnectToVeil_LengthInt' 
% 'ConnectToVeil_Density' 
% 'Branch2ndOrder_LengthInt'
% 'Validation'


switch filterType
    case 'ConnectToVeil_LengthInt';
        filterParams.filoTypes = [0,1]; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1; % turn on filtering by filopodia fit.
        filterParams.filterByBundleLength= [0.3,inf]; % for final going to filter by the length of the bundle...
        % this will automatically save short filo with embedded.
        filterParams.saveFiloByLengthAndSig = [];
        filterParams.embedFitCriteria = 95 ; % change threshold to be a percentile of the residuals
        filterParams.filoFitCriteria = 95;
        filterParams.filterBasedOnGroupUpstream = 0;
        filterParams.filterBasedOnBranchConnection =0 ;
        filterParams.filterIntNoConnect=0;
    case 'ConnectToVeil_Density';
        filterParams.filoTypes = [0,1]; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        
        % currently the first row is the min/max cut offs for the Ext Filo
        % the 2nd row is the upper and lower signal percentile
        % Filopodia that fall into both these categories will be saved.
        
        % for density you want to be slightly less rigid when filering out bad fits
        % bad fits can happen because of overlapping filo- want to include
        % these in the final calculation- do the same for orientation.
        
        filterParams.embedFitCriteria = 95 ; % change threshold to be a percentile of the residuals
        filterParams.filoFitCriteria = 95;
        filterParams.filterBasedOnGroupUpstream = 0;
        filterParams.filterBasedOnBranchConnection =0 ;
        filterParams.filterIntNoConnect=0;
        
    case 'Branch2ndOrder_LengthInt';
        filterParams.filoTypes = 2;
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3, inf];
        filterParams.saveFiloByLengthAndSig = [];
        filterParams.embedFitCriteria = 95 ; % Not relavent
        filterParams.filoFitCriteria = 95;
        filterParams.filterBasedOnGroupUpstream = 0;
        filterParams.filterBasedOnBranchConnection =0 ;
        filterParams.filterIntNoConnect=0;
        
    case 'Branch2ndOrder_Density_WithZeroOrder';
        filterParams.filoTypes = [0,1,2]; % 1st order attached to a veil with a branch, 2 branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        %filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        filterParams.saveFiloByLengthAndSig = [];
        filterParams.embedFitCriteria = 0 ;
        filterParams.filoFitCriteria = 95;
        filterParams.filterBasedOnGroupUpstream = 0;
        filterParams.filterBasedOnBranchConnection =0 ;
        filterParams.filterIntNoConnect=0;
        
    case 'Branch2ndOrder_Density_NoZero';
        filterParams.filoTypes = [1,2]; % 1st order attached to a veil with a branch, 2 branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        filterParams.embedFitCriteria = 0 ; %
        filterParams.filoFitCriteria = 95;
        filterParams.filterBasedOnGroupUpstream = 0;
        filterParams.filterBasedOnBranchConnection =0 ;
        filterParams.filterIntNoConnect=0;
        
    case 'Branch3rdOrder_LengthInt';
        filterParams.filoTypes = 3;
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.2, inf];
        filterParams.saveFiloByLengthAndSig = [];
        filterParams.embedFitCriteria = 0 ; %
        filterParams.filterBasedOnGroupUpstream = 0;
        filterParams.filterBasedOnBranchConnection =0 ;
        filterParams.filterIntNoConnect=0;
        
    case 'Branch3rdOrder_DensityOrient';
        filterParams.filoTypes = 3; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.2,inf];
        filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        filterParams.embedFitCriteria = 0 ;
        filterParams.filterBasedOnGroupUpstream = 0;
        filterParams.filterBasedOnBranchConnection =0 ;
        filterParams.filterIntNoConnect=0;
    case 'Validation';
        filterParams.filoTypes = [0 Inf]; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        % filterParams.saveFiloByLengthAndSig = [];
        filterParams.saveFiloByLengthAndSig = [5 inf; 50 100]; 
        filterParams.embedFitCriteria =95 ;
        filterParams.filoFitCriteria = 95;
        filterParams.filterBasedOnGroupUpstream = 0;
        filterParams.filterBasedOnBranchConnection = 0;
        filterParams.filterIntNoConnect=1;
        
    case 'Validation2'
        filterParams.filoTypes = [0 Inf]; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        filterParams.saveFiloByLengthAndSig = [];
        %filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        filterParams.embedFitCriteria =95 ;
        filterParams.filoFitCriteria = 95;
        filterParams.filterBasedOnGroupUpstream = 1;
        filterParams.filterBasedOnBranchConnection =1 ;
        filterParams.filterIntNoConnect=0;
        
    case 'Validation_NoEmbed_NoSave'    
        filterParams.filoTypes = [0 Inf]; % all 
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        %filterParams.filterByBundleLength = [0.15,inf];% for the Gupton mCherry movies
        filterParams.saveFiloByLengthAndSig = [];
        %filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        %filterParams.saveFiloByLengthAndSig = [5 inf; 40 100];
        filterParams.filoFitCriteria = 95;
        filterParams.filterBasedOnGroupUpstream = false;
        
        filterParams.filterBasedOnBranchConnection = false;

    case 'Validation_NoEmbed'
        filterParams.filoTypes = [0 Inf]; % all
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf]; 
        %filterParams.filterByBundleLength = [0.15,inf];% for the Gupton mCherry movies
        filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        %filterParams.saveFiloByLengthAndSig = [5 inf; 40 100];
        filterParams.filoFitCriteria = 95;
        filterParams.filterBasedOnGroupUpstream = false;
        
        filterParams.filterBasedOnBranchConnection = false;
        %% Extra Filts
%         filterParams.filterBasedOnGroupUpstream = true; % set to true for the gupton mCherry movies
%         filterParams.filterBasedOnBranchConnection = true; % set to true the gupton mCherry movies
%         % Original Defaults
%         filterParams.minPercentFiloAboveBack = 80; % default is 75 percent if empty will not filter
%         % minimum percentage of the filopodia that needs to be above
%         filterParams.nSTDParticleDetector = 3;
        %         filterParams.minPercentFiloAboveBack = 90;
        % filterParams.nSTDParticleDetector = 1;
        
    case 'curvVsLength'
        filterParams.filoTypes = 0; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        
        % filterParams.saveFiloByLengthAndSig = [];
        filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        
        filterParams.embedFitCriteria = 95 ; % change threshold to be a percentile of the residuals
        filterParams.filoFitCriteria = 95;
        filterParams.filterBasedOnGroupUpstream = 0;
        filterParams.filterBasedOnBranchConnection =0 ;
        filterParams.filterIntNoConnect=0;
    case 'ConnectToVeil_LengthInt_Biosensors'
        filterParams.filoTypes = 0; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1; % turn on filtering by filopodia fit.
        filterParams.filterByBundleLength= [0.3,inf]; % for final going to filter by the length of the bundle...
        % this will automatically save short filo with embedded.
        filterParams.saveFiloByLengthAndSig = [];
        filterParams.embedFitCriteria = 95 ; % change threshold to be a percentile of the residuals
        filterParams.filoFitCriteria = 95;
        filterParams.filterBasedOnGroupUpstream = 0;
        filterParams.filterBasedOnBranchConnection =0 ;
        filterParams.filterIntNoConnect=0;
        
        % filter based on local threshold around filopodia in FRET/Donor channel
        filterParams.sigmaNumerator = 3; % multiple times the standard deviation of the background estimation
        % for the numerator of a FRET ratio - typically FRET channel
        
        filterParams.sigmaDenominator = 1; % multiple time the standard deviation of the background estimation
        % for the denominator of the FRET ratio = typically the Donor
        % channel - as this value should decrease upon signal activity
        % this should likely not be quite as stringent
        
        filterParams.backSampleSize = 20; % the sample size for the background has to be at least this large
        filterParams.excludeCrosses = false;  %need to implement this...

end

