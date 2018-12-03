function [ h ] = GCAVisualsProtrusionMap(imgData, varargin)
%GCAVisualsProtrusionMap

%%Input check
ip = inputParser;

ip.CaseSensitive = false;

ip.addRequired('imgData'); % what windowing method to use - currently assumes you will use one of the
% windowing methods associated with the MD.

ip.addParameter('frameInSec',5);

ip.addParameter('umPerMin',false); % default assumed in nm/sec

ip.addParameter('CLim',[-100,100]); % currently in nm/s

% ip.addParameter('blackOutHi',false);
% 
% ip.addParameter('blackOutLo',false);

ip.addParameter('colorMap','jet');

ip.addParameter('FontSize',20);

ip.addParameter('ShowColorBar',true);
ip.addParameter('ColorBarForRatio',false);
ip.addParameter('ScaleFactor',1000);

ip.addParameter('alphaMaskNaN',false);

ip.addParameter('title',[]);
ip.addParameter('visible','on');

ip.parse(imgData,varargin{:});
%%
frameInSec = ip.Results.frameInSec;
h = setupFigure( 'DisplayMode', 'screen');
%h = setAxis(ip.Results.visible,0.95,ip.Results.FontSize);
x = ip.Results.colorMap;
% if ip.Results.blackOutLo
%     % colormap set to go from lowest to highest
%     x(1,:) = [1,1,1]; % set the original data to white
% end
% if ip.Results.blackOutHi
%     x(end,:) = [0,0,0]; % set the outlier data to black
% end
colormap(x)

h = imagesc(imgData,ip.Results.CLim);

axis xy
% white out NaN values
%                              alphamask =true(size(imgData));
%                              alphamask(isnan(imgData))=false;
%                              set(h,'AlphaData',alphamask,'AlphaDataMapping','none');
ylabel('Window Number');
%
if ip.Results.umPerMin
    forLabel = 'min';
else   
    forLabel = 'sec';
end

xlabel(['Time ' '( ' forLabel ' )'] )

if ip.Results.ShowColorBar
    c = colorbar;
    if ip.Results.ColorBarForRatio
        c.Label.String = 'FRET/Donor';
        v = c.TickLabels;
        new = cellfun(@(x) (str2double(x)./ip.Results.ScaleFactor),v,'uniformoutput',0);
        c.TickLabels = new;
    end
end
[ny,nx] = size(imgData);
axis([1,nx-1,1,ny]); 
frameNum = get(gca,'XTickLabel');
converted = str2double(frameNum).*frameInSec; % SET THE TIME INTERVAL
num2str(converted);
set(gca,'XTickLabel',converted);

if ip.Results.alphaMaskNaN
    % mask out NaNs
    alphaMask = false(size(imgData));
    alphaMask(isnan(imgData)) = true;
    set(gcf,'AlphaMap', [1,0.8]);
    % first number is will set the transparancy of
    % of the true values in the alphaMask while the second number will set the
    % transparancy of the false number of the alpha mask
    set(h, 'alphaData', alphaMask,'alphaDataMapping','direct');
end

if ~isempty(ip.Results.title)
    title(ip.Results.title);
end
%% This was a quick and dirty to reset the axis
%    h1 = get(gcf,'CurrentAxes');
%             yLim = h1.YLim; % use whatever they used
%             axis([0.5 29.5 yLim(1) yLim(2)]);
%%
end % function 