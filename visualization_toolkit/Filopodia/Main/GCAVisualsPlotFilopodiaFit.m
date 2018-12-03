function [ output_args ] = GCAVisualsPlotFilopodiaFit( img,filoInfo,varargin)
% GCAVisualsPlotFilopodiaFit
%
%
%%Input check
ip = inputParser;
ip.CaseSensitive = false;
ip.addParameter('OutputDirectory',[]);
ip.addParameter('filoPart','Ext_');
ip.addParameter('idxToPlot',[]);
ip.addParameter('localLims',true);
ip.addParameter('saveOutput',false);

ip.parse(varargin{:});

%% Set up
if ip.Results.saveOutput
    if  isempty(ip.Results.OutputDirectory)
        outDir = [pwd filesep 'FilopodiaFits'];
    else
        outDir = ip.Results.OutputDirectory;
    end
    if ~isdir(outDir)
        mkdir(outDir)
    end
end

if isempty(ip.Results.idxToPlot)
    idxToPlot = [1:length(filoInfo)]; % plot them all
else
    idxToPlot = ip.Results.idxToPlot;
end
purple = [0.5020, 0.4510, 0.6745];
cyan = [ 0.0039  ,  0.4264 ,   0.3848];
dred = [ 0.6471 ,        0 ,   0.1490  ];
std3 =  [0.7804    0.9176    0.8980];
%std3 =  [0.5020    0.8039    0.7569];
std2 =  [0.3529    0.7059    0.6745];


%% START
for ifilo = 1:length(idxToPlot)
    idxCurrent = idxToPlot(ifilo);
    maskIndices = filoInfo(idxCurrent).([ip.Results.filoPart 'maskIndices']);
    maskIndices = maskIndices(:,[5,3,1,2,4]);
    maskIndices = maskIndices'; % switch for plotting
    % load the distance coords
    distFilo =  filoInfo(idxCurrent).([ip.Results.filoPart 'distFilo']);
    
    yData =  filoInfo(idxCurrent).([ip.Results.filoPart 'weightedAvg']);
    
    
    params = filoInfo(idxCurrent).([ip.Results.filoPart 'params']);
    if ~isnan(params)
        yFit = params(1)*(1-0.5*(1+erf((distFilo-params(2))/(params(3)*sqrt(2)))))+params(4);
    end
    resid = filoInfo(idxCurrent).([ip.Results.filoPart 'resid']);
    
    % note in this output I dont' GCAVisualsPlotFilopodiaFit( img,filoInfo,varargin)currently save the slopeMaxNeg and the
    % slopeMaxPos I also dont' save the start and end fits
    
    fsFigure(0.75,'visible','on');
    
    %% First Subplot: View Local Detection of Ridge Signal
    if sum(isnan(maskIndices(:)))==0; % if no-mask NaNs (ie not at the border);
        subplot(3,2,(1:2));
        imgFilo = img(maskIndices);
        
        if ip.Results.localLims
            lims = [-max(yData(:)),-min(yData(:))];
        else
            lims = [-max(img(:)),-min(img(:))];
        end
        
        imagesc(-imgFilo,lims);
        
        colormap('gray')
        % mark pixels for potential background
        % estimation
        forSearch =filoInfo(idxCurrent).([ip.Results.filoPart 'pixIndicesFor']);
        nSearch = length(forSearch(~isnan(forSearch)));
        backMask = zeros(size(imgFilo));
        backMask([1:2,4:5],end-nSearch:end) = 1;
        idx = find(backMask);
        
        if isfield(filoInfo,[ip.Results.filoPart 'indicesRemoveBackground']);
            indicesRemove = filoInfo(idxCurrent).([ip.Results.filoPart 'indicesRemoveBackground']);
            
            if ~isempty(indicesRemove);
                idxRemove = arrayfun(@(x) find(maskIndices==indicesRemove(x)),1:length(indicesRemove),'uniformoutput',0);
                idxRemove = vertcat(idxRemove{:});
                idx = setdiff(idx,idxRemove);
            end
        end % isfield
        [y,x] = ind2sub(size(backMask),idx);
        hold on
        scatter(x,y,100, [0.5020, 0.4510, 0.6745],'x'); % color in purple for now
        
        ylabel('Filopodia Width (Pixels)');
        axis([0.5,size(maskIndices,2),0.5,5.5])
        title(['Filopodia Detection ' num2str(idxCurrent) ' : Purple Crosses Mark Local Background'],'FontSize',10,'FontName','Arial');
    end
    
    %% Second Subplot: Intensity Values Along Automated Linescan Etc
    subplot(3,2,(3:4));
    hData = scatter(distFilo,yData,50,'k');
    hold on
    %%
    % mark the estimated local background value
    backEstMean = nanmean(filoInfo(idxCurrent).([ip.Results.filoPart 'localBackValues']));
    hBack = line([0 distFilo(end)], [backEstMean,backEstMean],'color','k','Linestyle','--');
    dataFit = filoInfo(idxCurrent).([ip.Results.filoPart 'coordsForFit']);
    % mark the data used for the fit
    hDataFit =  scatter(dataFit(:,1),dataFit(:,2),'k','filled'); % color in the data specifically used for the fitting.
    %%
    % plot the fit
    hFit =  plot(distFilo,yFit,'r');
    
    hMean = line([params(2),params(2)],[min(yData),max(yData)],'Color','k','Linewidth',2);
    %% plot the sigma
    
    t = get(gcf,'Children');
    gca = t(1);
    ylim = get(gca,'YLim');
    s = filoInfo(idxCurrent).([ip.Results.filoPart 'params'])(3);
    m = filoInfo(idxCurrent).([ip.Results.filoPart 'params'])(2);
    a = filoInfo(idxCurrent).([ip.Results.filoPart 'params'])(1)+filoInfo(idxCurrent).([ip.Results.filoPart 'params'])(4);
    a2 = filoInfo(idxCurrent).([ip.Results.filoPart 'params'])(1)/2;
    yForVar = a-a2;
    params = filoInfo(idxCurrent).([ip.Results.filoPart 'params']);
    startVar3 = params(2) + (3*params(3))/2;
    stopVar3 = params(2) - (3*params(3))/2;
    yForVar = params(1)/2 + params(4); % half the amplitude + the background estimation
    hVar3 = line([stopVar3 startVar3],[yForVar,yForVar],'Color',std3,'Linewidth',3);
    
    startVar2 = params(2) + (params(3)*2)/2;
    stopVar2 = params(2) - (params(3)*2)/2;
    
    hVar2 = line([stopVar2 startVar2],[yForVar,yForVar],'Color',std2,'Linewidth',3);
    
    startVar1 = params(2) + params(3)/2;
    stopVar1 = params(2) - params(3)/2;
    
    hVar1 = line([stopVar1 startVar1],[yForVar,yForVar],'Color',cyan,'Linewidth',3);
    %line(gca,[m+3*(s/2),m-3*(s/2)],[a-a2,a-a2],'color','k')
    %line(gca,[m,m],[ylim(1),ylim(2)],'color','k');
    
    %%
    % plot sites of potential decay and positive slopes
    % for troubleshooting
    negPlot = filoInfo(idxCurrent).([ip.Results.filoPart 'slopeMaxNeg']);
    posPlot = filoInfo(idxCurrent).([ip.Results.filoPart 'slopeMaxPos']);
    
    hDec = scatter(negPlot(:,1),negPlot(:,2),50,[ 0.0039  ,  0.4264 ,   0.3848],'filled'); % cyan
    
    hInc = scatter(posPlot(:,1),posPlot(:,2),50,[ 0.6471 ,        0 ,   0.1490  ] ,'filled');
    %%
    warning('off',('MATLAB:legend:IgnoringExtraEntries'));
    %      legend([hData,hDataFit,hFit,hDec,hBack,hMean,hInc],'Raw Data', 'Data For Fit','Fit','Potential Signal Decay',...
    %          'Background Estimate','Tip Position','Potential Signal Rise','location','northeast');
    %      legend('boxoff');
    xlabel('Distance Along Filopodia (Pixels)','FontSize',10,'FontName','Arial');
    ylabel('Fluorescence Intensity (AU)','FontSize',10,'FontName','Arial');
    %      if exitFlag>1
    %          title({['Length (Mean Gaussian Survival) = ' num2str(params(2).*0.216,3) ' um'] , ['Intensity Filament (Amplitude) = ' num2str(params(1,1),3)] })
    %      end
    axis( [0,max(distFilo),min(yData),max(yData)]);
    %% Third Subplot: Filopodia Overlay
    subplot(3,2,5);
    imshow(-img,[]) ;
    hold on
    filoInfoC = filoInfo(idxCurrent);
    text(5,5,['Filopodia ' num2str(idxCurrent)],'FontSize',6,'FontName','Arial');
    
    pixelsF = zeros(size(img));
    pixF= filoInfoC.([ip.Results.filoPart 'pixIndicesFor']);
    pixF = pixF(~isnan(pixF));
    pixelsF(pixF)=1;
    % pixelsB = zeros(size(img));
    % pixelsB(filoInfoC.([toAdd 'pixIndicesBack']))=1;
    [yb,xb] = ind2sub(size(img), filoInfoC.([ip.Results.filoPart 'pixIndicesBack']));
    scatter(xb,yb,10,[ 0.0039  ,  0.4264 ,   0.3848],'filled');
    %                         spy(pixelsF,'g',10);
    %                         spy(pixelsB,'r');
    [yf,xf] = ind2sub(size(img),filoInfoC.([ip.Results.filoPart 'pixIndicesFor']));
    scatter(xf,yf,10,[0.5020, 0.4510, 0.6745],'filled');  % purple
    %% 4th Subplot the residuals compared to the amplitude
    subplot(3,2,6);
    boxplot(resid);
    hold on
    thresh = prctile(resid,95);
    if params(1) > thresh;
        add = '*';
    else
        add = [];
    end
    
    if ip.Results.localLims
        yLimHi = max(yData);
    else
        yLimHi = max(img(:));
    end
    yLimLo = min(resid)-2;
    axis([0.5 1.5 yLimLo yLimHi]);
    
    scatter(1,params(1),'r','filled');
    text(1,params(1),['Intensity Filament ' add],'FontSize',6,'FontName','Arial');
    ylabel({'Residuals' ; 'From Local Fit'},'FontSize',6,'FontName','Arial');
    
    if ip.Results.saveOutput
        
        filename{1} = ['Filopodia_' num2str(idxCurrent,'%03d') ip.Results.filoPart '.fig' ];
        filename{2} = ['Filopodia_' num2str(idxCurrent,'%03d') ip.Results.filoPart '.png'];
        for i = 1:2
            saveas(gcf, [ outDir  filesep filename{i}]);
        end
        warning('off','MATLAB:print:FigureTooLargeForPage')
        saveas(gcf,[outDir filesep 'Filopodia_' num2str(idxCurrent,'%03d') ip.Results.filoPart '.psc2'],'psc2');
        close gcf
    end % if saveOutput
end

end