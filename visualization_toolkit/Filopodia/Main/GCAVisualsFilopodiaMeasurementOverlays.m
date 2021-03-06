function [mask] = GCAVisualsFilopodiaMeasurementOverlays(filoInfo,imgSize,varargin)
%GCAVisualsFilopodiaOverlaysFromFilterSet
% INPUT:
% filoInfo
% imgSize
% filoFilterSet:  set of
% justExt
%% INPUTPARSER
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
ip.addRequired('filoInfo');
ip.addRequired('imgSize');


ip.addParameter('filoFilterSet',[]); % note keep the filter set here as it makes it easier to plot original IDs if user desires
ip.addParameter('plotValues',[]); % values per filopodia for plotting and or color/coding values : if set to 'IDs' will plot the
% filo's order in the data structure

ip.addParameter('justExt',1); % flag for internal filo versus ext filo plotting

% flag to plotSigma as well as mean
ip.addParameter('plotSigma',[]); % if empty don't plot default 1/2 sigma (as plotting positive side only)

% FiloColor Options
ip.addParameter('colorEmbed',[0 0 1]);
ip.addParameter('colorFiloBranch',[0 0 0]);
ip.addParameter('colorBranch',[0 0 0]); % only applicable in branch mode

% cMap Parameters: Options related to color coding by value
ip.addParameter('colorByValue',false);
ip.addParameter('extraColor',[0,0,0]); % add an extra color to the cMap
ip.addParameter('cMapLimits',[]); % [min max] limits only applied if colorByValue is true
ip.addParameter('colorMap',[]);

ip.addParameter('forceFiloText',0); %
ip.addParameter('branchMode',false); % branch mode will by default color code by type
ip.addParameter('plotText',true,@(x) islogical(x));

ip.addParameter('plotTextAtBranches',false,@(x) islogical(x)); % just make an explicit input for now instead of testing
ip.addParameter('plotFiltered',false);

% Options related to plotting
ip.addParameter('LineWidth',1);
ip.addParameter('MarkerSize',10);
ip.addParameter('LineStyle','-');

% Options related to smoothing
ip.addParameter('UseSmoothedCoords',false);

% output mask pixels.
ip.addParameter('createMask',false);

ip.parse(filoInfo,imgSize,varargin{:});
%% Initiate

% Filters
if isempty(ip.Results.filoFilterSet);
    filoFilterSet = true(length(filoInfo),1);
else
    filoFilterSet = ip.Results.filoFilterSet;
    filoFilterSet = filoFilterSet(:,1);
end

if strcmpi(ip.Results.plotValues,'IDs')
    plotValues = 1:length(filoInfo); % plot the ID number
    plotValues = plotValues';
    plotValues = plotValues(filoFilterSet);
else
    plotValues = ip.Results.plotValues;
end

IDsCurrentSet = find(filoFilterSet);

% filter filoInfo
filoInfoFilt = filoInfo(filoFilterSet(:,1));

% cumbersome historical flags for plotting embedded vs extra-veil filo
switch ip.Results.justExt
    case 1
        typeEnd =1;
        typeStart =1;
    case 2
        typeEnd = 2;
        typeStart = 2;
    case 3
        typeStart  = 1;
        typeEnd = 2;
end

toAdd{1} = 'Ext_';
toAdd{2} = 'Int_';

% initiate the mask
mask = zeros(imgSize);

% get the mapper if applicable
if ip.Results.colorByValue && ~isempty(plotValues)
    
    if ~isempty(ip.Results.cMapLimits);
        minValue = ip.Results.cMapLimits(1);
        maxValue = ip.Results.cMapLimits(2);
    else
        
        minValue = min(plotValues);
        maxValue = max(plotValues);
    end
    
    if isempty(ip.Results.colorMap)
        cMapLength=128; % cMap=jet(cMapLength);
        cMap = jet(128);
        %         cMap = brewermap(128,'spectral');
        %         cMap = flip(cMap,1);
    else
        cMap = ip.Results.colorMap;
        cMapLength = size(cMap,1);
    end
    mapper=linspace(minValue,maxValue,cMapLength)';
    D=createDistanceMatrix(plotValues,mapper);
    [sD,idxCMap]=sort(abs(D),2);
    if ~isempty(ip.Results.extraColor);
        cMap = [ip.Results.extraColor; cMap]; % assumes adding the color to the first value which is blue in jet
    end
    % make a cmap for sigma this brightened with respect to the previous
    % color map
    if ~isempty(ip.Results.plotSigma)
        cMapSigma = brighten(cMap,0.7);
    end
end

c{1} = [ip.Results.colorFiloBranch];
c{2} = [ip.Results.colorEmbed];

% parition data if in branch mode
if ip.Results.branchMode % will use ip.Results.colorFiloBranch as the stem
    % find the smallest type
    types = vertcat(filoInfoFilt(:).type);
    types(types==0) = 1; % treat any no branch filo as stem
    if length(unique(types))>2;
        display(['You are plotting in branch mode, but' ...
            'but the filter has more than two filoTypes: Check FilterSet']);
    end
    
    stemType = min(unique(types));
    
    filoInfoBranch = filoInfoFilt(types~=stemType);
    filoInfoFilt= filoInfoFilt(types==stemType); % set the filter info equal to the stem
    % test if the values equal the number of branches or the nubmer of
end % if ip.Results.branchMode
%% START Loop over the embedded/extra veil filo if user desires.
for iType = typeStart:typeEnd
    colorC = c{iType};
    
    if ~isempty(filoInfoFilt)
        
        for i = 1:numel(filoInfoFilt)
            value = filoInfoFilt(i).([toAdd{iType} 'endpointCoordFitXY']);
            if ~isnan(value);
                xycoordsEndFit(i,:) = filoInfoFilt(i).([toAdd{iType} 'endpointCoordFitXY']);
            else
                xycoordsEndFit(i,:) = [NaN,NaN];
            end
        end
        
        if ip.Results.UseSmoothedCoords
            fieldname = '_SplineFit';
        else
            fieldname = '';
        end
        
        % make sure to take out any NaN
        toRemove1 = arrayfun(@(x) numel(x.([toAdd{iType} 'coordsXY' fieldname])),filoInfoFilt);
        toRemove2= arrayfun(@(x) isempty(x.([toAdd{iType} 'coordsXY' fieldname ])),filoInfoFilt);
        toKeep = ~(toRemove1==1 | toRemove2 ==1 );
        
        xCoordsBase= arrayfun(@(x) x.([toAdd{iType} 'coordsXY' fieldname])(1,1),filoInfoFilt(toKeep));
        yCoordsBase = arrayfun(@(x) x.([toAdd{iType} 'coordsXY' fieldname])(1,2),filoInfoFilt(toKeep));
        xycoordsBase = [xCoordsBase' yCoordsBase'];
        xycoordsEndFit =  xycoordsEndFit(toKeep,:);
        
        % currently it is easiest now just to plot separately...
        filoIDs = 1:length(filoInfoFilt);
        filoIDs = filoIDs(toKeep);
        if ip.Results.colorByValue && ~isempty(plotValues)
            for k=1:cMapLength
                filoToPlot = find(idxCMap(:,1)==k);
                filoToPlot = intersect(filoIDs,filoToPlot);
                
                if ~isempty(filoToPlot)
                    for ifilo = 1:size(filoToPlot,1)
                        test = find(filoIDs == filoToPlot(ifilo));
                        if ip.Results.UseSmoothedCoords
                            xy=  filoInfoFilt(filoToPlot(ifilo)).([toAdd{iType} 'coordsXY' fieldname]);
                            xC = xy(:,1);
                            yC = xy(:,2);
                            
                        else
                            pixIndices = filoInfoFilt(filoToPlot(ifilo)).([toAdd{iType} 'pixIndices']);
                            idxEnd = find(pixIndices == filoInfoFilt(filoToPlot(ifilo)).([toAdd{iType} 'endpointCoordFitPix']));
                            pixIndicesPlot = pixIndices(1:idxEnd);
                            if ip.Results.createMask
                                mask(pixIndicesPlot) = 1;
                            end
                            [yC,xC] = ind2sub( imgSize  ,pixIndicesPlot);
                            
                            if ~isempty(ip.Results.plotSigma)
                                % note eventually might want to do this
                                % calc and save in filo branch structure
                                sigma =  filoInfoFilt(filoToPlot(ifilo)).([toAdd{iType} 'params'])(2) ...
                                    + filoInfoFilt(filoToPlot(ifilo)).([toAdd{iType} 'params'])(3)*ip.Results.plotSigma;
                                dist = filoInfoFilt(filoToPlot(ifilo)).([toAdd{iType} 'distFilo']);
                                delt = (abs(dist-sigma));
                                idxS = find(delt==min(delt));
                                pixIndicesSig = pixIndices(idxEnd:idxS);
                                [yS,xS] = ind2sub(imgSize,pixIndicesSig);
                                
                            end
                            
                        end
                        
                        plot(xC,yC,'color',cMap(k,:),'Linewidth',ip.Results.LineWidth,'LineStyle',ip.Results.LineStyle);
                        if ~isempty(ip.Results.plotSigma)
                            plot(xS,yS,'color',cMapSigma(k,:), 'Linewidth', ip.Results.LineWidth,'LineStyle',ip.Results.LineStyle);
                        end
                        
                        scatter(xycoordsEndFit(test,1),xycoordsEndFit(test,2),ip.Results.MarkerSize,cMap(k,:),'filled');
                        scatter(xycoordsBase(test,1),xycoordsBase(test,2),ip.Results.MarkerSize, cMap(k,:),'filled');
                        
                        
                        if ~isempty(ip.Results.plotSigma)
                            scatter(xS(end),yS(end),ip.Results.MarkerSize,cMapSigma(k,:),'filled');
                        end
                        
                    end % ifilo
                end
            end % for k
        else % don't plot by colorValues
            
            for i = 1:numel(filoInfoFilt)
                if ip.Results.UseSmoothedCoords
                    xy=  filoInfoFilt(filoToPlot(i)).([toAdd{iType} 'coordsXY' fieldname]);
                    xC = xy(:,1);
                    yC = xy(:,2);
                    
                else
                    pixIndices = filoInfoFilt(i).([toAdd{iType} 'pixIndices']);
                    idxEnd = find(pixIndices == filoInfoFilt(i).([toAdd{iType} 'endpointCoordFitPix']));
                    pixIndicesPlot = pixIndices(1:idxEnd);
                    if ip.Results.createMask
                        mask(pixIndicesPlot) = 1;
                    end
                    [yC,xC] = ind2sub( imgSize  ,pixIndicesPlot);
                end
                plot(xC,yC,'color',colorC,'Linewidth',ip.Results.LineWidth,'LineStyle',ip.Results.LineStyle);
            end
            
            scatter(xycoordsEndFit(:,1),xycoordsEndFit(:,2),ip.Results.MarkerSize,colorC,'filled');
            scatter(xycoordsBase(:,1),xycoordsBase(:,2),ip.Results.MarkerSize, colorC,'filled');
            clear  idxEnd pixIndicesPlot xycoordsEndFit
        end
    end % isempty
end    % iType
%% Plot Branches if applicable
% Ok the quickest way to do this right now is to repeat the above for branches - can maybe re-code if have time
%  the other option would be just to loop over the entire set with one
% with one color and to only turn on the embed flag if user wants it .
% could just make below a subfunction
if ip.Results.branchMode
    
    for i = 1:numel(filoInfoBranch)
        value = filoInfoBranch(i).('Ext_endpointCoordFitXY');
        if ~isnan(value);
            xycoordsEndFit(i,:) = filoInfoBranch(i).('Ext_endpointCoordFitXY');
        else
            xycoordsEndFit(i,:) = [NaN,NaN];
        end
        
    end
    
    xCoordsBase= arrayfun(@(x) x.('Ext_coordsXY')(1,1),filoInfoBranch);
    yCoordsBase = arrayfun(@(x) x.('Ext_coordsXY')(1,2),filoInfoBranch);
    xycoordsBase = [xCoordsBase' yCoordsBase'];
    
    % plot the other pixels associated with the filo
    for i = 1:numel(filoInfoBranch)
        pixIndices = filoInfoBranch(i).( 'Ext_pixIndices');
        idxEnd = find(pixIndices == filoInfoBranch(i).( 'Ext_endpointCoordFitPix'));
        pixIndicesPlot = pixIndices(1:idxEnd);
        [yC,xC] = ind2sub( imgSize  ,pixIndicesPlot);
        plot(xC,yC,'color',ip.Results.colorBranch,'Linewidth',ip.Results.LineWidth);
    end
    scatter(xycoordsEndFit(:,1),xycoordsEndFit(:,2),ip.Results.MarkerSize,ip.Results.colorBranch,'filled');
    scatter(xycoordsBase(:,1),xycoordsBase(:,2),ip.Results.MarkerSize, ip.Results.colorBranch,'filled');
    clear  idxEnd pixIndicesPlot xycoordsEndFit
end % branchMode
%% plot the filtered filo segments if user desires
if ip.Results.plotFiltered
    
    if   ~isempty(filoInfo)
        [toPlotY,toPlotX]  =  arrayfun(@(x) ind2sub(imgSize,x.([toAdd{iType} 'pixIndices'])),filoInfo,'uniformoutput',0);
        
        for i = 1:numel(toPlotY)
            plot(toPlotX{i}(:),toPlotY{i}(:),'color','k','Linewidth',ip.Results.LineWidth);
        end
    end
end
%% Plot Text if applicable
if ~isempty(plotValues) && ip.Results.plotText
    plotValues = plotValues';
    if (length(plotValues) > 1 || ip.Results.forceFiloText)
        if ip.Results.plotTextAtBranches
            % arrayfun(@(i,j) text(i.Ext_coordsXY(1,1),i.Ext_coordsXY(1,2),num2str(j,3),'color','k'),filoInfoBranch,plotValues);
            % you would have to filter again the by the conIdx
            [~,idxKeep] = arrayfun(@(x) intersect(filoInfoFilt(x).conIdx,IDsCurrentSet),1:length(filoInfoFilt),'uniformoutput',0);
            coordsAndD = arrayfun(@(x) filoInfoFilt(x).conXYCoords(idxKeep{x},:),1:length(idxKeep),'uniformoutput',0);
            toPlot = vertcat(coordsAndD{:});
            arrayfun(@(x) text(toPlot(x,1),toPlot(x,2),num2str(toPlot(x,3).*0.216,3),'color','k','FontSize',10),1:size(toPlot,1));
        else
            
            arrayfun(@(i,j) text(i.Ext_coordsXY(1,1),i.Ext_coordsXY(1,2),num2str(j,3),'color','k','FontSize',6),filoInfoFilt,plotValues);
        end
    else
        text(10,5,num2str(plotValues,3)); % single value parameter not per filopodia
    end
end % isempty(plotValues)
mask = logical(mask);
end % function