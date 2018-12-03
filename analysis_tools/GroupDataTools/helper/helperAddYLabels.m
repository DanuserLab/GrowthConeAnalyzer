function [ toPlot ] = helperAddYLabels( toPlot )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

meas = fieldnames(toPlot); 

for imeas = 1:numel(meas)
    measC = meas{imeas};
   % if ~isfield(toPlot.(meas{imeas}),'yLabel');
        
        switch measC
            %% Veil
            case 'VeilStemThickness'
                
                toPlot.(meas{imeas}).yLabel  = {'Veil/Stem Thickness ' ; '(um)'};
                toPlot.(meas{imeas}).yGroup = 'Veil';
            case 'VeilStemArea' 
                toPlot.(meas{imeas}).yLabel = {'Veil Area' ; 'um'}; 
                toPlot.(meas{imeas}).yGroup = 'Veil'; 
            case 'VeilStemAspectRatio' 
                toPlot.(meas{imeas}).yLabel = {'Veil Aspect Ratio'}; 
                toPlot.(meas{imeas}).yGroup = 'Veil'; 
                
                %%  Actin Bundle
            case 'filoLengthToVeil'
                toPlot.(meas{imeas}).yLabel = {'Filopodia Length to Veil' ; '(um)'};
                toPlot.(meas{imeas}).yGroup = 'ActinBundle';
            case 'filoLengthEmbedded' 
                toPlot.(meas{imeas}).yLabel = {'Filopodia Length Embedded in Veil' ; '(um)' ; 'Non-Zero Only'};
                toPlot.(meas{imeas}).yGroup = 'ActinBundle'; 
            case 'filoLengthFullActinBundle'
                toPlot.(meas{imeas}).yLabel = {'Full Actin Bundle Length' ; '(um)'};
                toPlot.(meas{imeas}).yGroup = 'ActinBundle';
            case 'filoMaxCurvature' 
                toPlot.(meas{imeas}).yLabel = {'Filopodia Maximum Curvature'; '(1/pixels)'}; 
                toPlot.(meas{imeas}).yGroup = 'ActinBundle'; 
            case 'filoDensityAlongVeil' 
                toPlot.(meas{imeas}).yLabel = {'Filopodia Density Along Veil' ; '(Number of Filopodia per 10 um)'};
                toPlot.(meas{imeas}).yGroup = 'ActinBundle'; 
            case 'filoOrientation' 
                toPlot.(meas{imeas}).yLabel = {'Filopodia Orientation Relative to Veil/Stem' ; '(Degrees)'}; 
                toPlot.(meas{imeas}).yGroup = 'ActinBundle'; 
            case 'filoIntensityEmbedded_Norm'
                toPlot.(meas{imeas}).yLabel = {'Veil Embedded Actin Bundle'; 'LifeAct Fluorescence Intensity' ; '(Normalized to Non-Bundle Veil Intensity)'}; 
                toPlot.(meas{imeas}).yGroup = 'ActinBundle'; 
                
            case 'filoIntensityToVeil_Norm' 
                toPlot.(meas{imeas}).yLabel = {'Filopodia Bundle to Veil' ; 'LifeAct Fluorescence Intensity' ; '(Normalized to Non-Bundle Veil Intensity)'}; 
                toPlot.(meas{imeas}).yGroup = 'ActinBundle'; 
                
            case 'percentEachActinBundleEmbed'
                toPlot.(meas{imeas}).yLabel = {'Percent Length of Each' ; 'Actin Bundle Embedded in Veil'};
                toPlot.(meas{imeas}).yGroup = 'ActinBundle'; 
                
            case 'percentTotalActinBundlesVeilEmbedded'
                toPlot.(meas{imeas}).yLabel = {'Percentage of the Total' ; 'Actin Bundles Veil Embedded'}; 
                toPlot.(meas{imeas}).yGroup =  'ActinBundle'; 
                
             case 'filoIntensityToVeil' 
                toPlot.(meas{imeas}).yLabel = {'Filopodia Bundle to Veil' ; 'LifeAct Fluorescence Intensity' ; '(Non-Bundle Veil Intensity)'}; 
                toPlot.(meas{imeas}).yGroup = 'ActinBundle';  
              case 'filoIntensityEmbedded'
                toPlot.(meas{imeas}).yLabel = {'Veil Embedded Actin Bundle'; 'LifeAct Fluorescence Intensity' ; '(Non-Bundle Veil Intensity)'}; 
                toPlot.(meas{imeas}).yGroup = 'ActinBundle';   
                
                %% Branch Measurements
            case 'branchDistanceFromVeil'
                toPlot.(meas{imeas}).yLabel = {'Branch (2nd Order) Distance from Veil ' ; '(um)'};
                toPlot.(meas{imeas}).yGroup = 'Branch';
            case 'branchDensity_2ndOrder' 
                toPlot.(meas{imeas}).yLabel = {'Branch (2nd Order) Density ' ; '(Actin Branch per 10 um)' };  
                toPlot.(meas{imeas}).yGroup = 'Branch'; 
            case 'branchIntensity_2ndOrder' 
                toPlot.(meas{imeas}).yLabel = {'Branch (2nd Order) LifeAct Fluorescence Intensity' ; '(Normalized to Non-Bundle Veil Intensity)' };  
                toPlot.(meas{imeas}).yGroup = 'Branch';
            case 'branchMaxCurvature_2ndOrder' 
                toPlot.(meas{imeas}).yLabel = {'Branch (2nd Order) Maximum Curvature' ; '(1/pixels)'}; 
                toPlot.(meas{imeas}).yGroup = 'Branch'; 
            case 'branchLength_2ndOrder'
                toPlot.(meas{imeas}).yLabel = {'Branch (2nd Order) Length ' ; '(um)'}; 
                toPlot.(meas{imeas}).yGroup = 'Branch'; 
            case 'branchOrientation_2ndOrder'
                toPlot.(meas{imeas}).yLabel = {'Branch (2nd Order) Orientation' ; '(Degrees)'};
                toPlot.(meas{imeas}).yGroup = 'Branch'; 
            case 'filoBranchComplexity'
                toPlot.(meas{imeas}).yLabel = {'Filo/Branch Network Complexity' ; '(Number of branches per 10 um)'}; 
                toPlot.(meas{imeas}).yGroup = 'Branch'; 
                
            case 'retractionAnalysis_Veloc_' 
                toPlot.(meas{imeas}).yGroup = 'Veil'; 
            case 'protrusionAnalysis_Veloc_' 
                toPlot.(meas{imeas}).yGroup = 'Veil'; 
               
            case 'retractionAnalysis_persTime_greaterThan75th_Perc'
                toPlot.(meas{imeas}).yGroup= 'Veil'; 
            case  'protrusionAnalysis_persTime_greaterThan75th_Perc'
               toPlot.(meas{imeas}).yGroup = 'Veil'; 
            case 'info' 
                toPlot.(meas{imeas}).yGroup = 'A'; 
            case 'ExpressionNormalizationBS'; 
                toPlot.(meas{imeas}).yGroup = 'Veil'; 
                toPlot.(meas{imeas}).yLabel = {'Veil F-actin Intensity (AU)' ; 'No Bundles :  Background Subtracted'}; 
            case 'ExpressionNormalization'
                toPlot.(meas{imeas}).yGroup = 'Veil'; 
                toPlot.(meas{imeas}).yLabel = {'Veil F-actin Intensity (AU)' ;  'No Bundles'}; 
                
        end   % switch measC
    %end
end


end

