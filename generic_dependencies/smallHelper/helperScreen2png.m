
function helperScreen2png(filename,varargin)
%SCREEN2JPEG Generate a JPEG file of the current figure with
%   dimensions consistent with the figure's screen dimensions.
%
%   SCREEN2JPEG('filename') saves the current figure to the
%   JPEG file "filename".
%
%    Sean P. McCarthy
%    Copyright (c) 1984-98 by MathWorks, Inc. All Rights Reserved
% From: https://www.mathworks.com/matlabcentral/answers/102382-how-do-i-specify-the-output-sizes-of-jpeg-png-and-tiff-images-when-using-the-print-function-in-mat
%% add optional 
if nargin < 1
     error('Not enough input arguments!')
end

ip = inputParser;
ip.CaseSensitive = false;
ip.addParameter('figureHandle',gcf); 
ip.parse(varargin{:});

%%
oldscreenunits = get(ip.Results.figureHandle,'Units');
oldpaperunits = get(ip.Results.figureHandle,'PaperUnits');
oldpaperpos = get(ip.Results.figureHandle,'PaperPosition');
set(ip.Results.figureHandle,'Units','pixels');
scrpos = get(ip.Results.figureHandle,'Position');
newpos = scrpos/100;
set(ip.Results.figureHandle,'PaperUnits','inches',...
     'PaperPosition',newpos)
print(ip.Results.figureHandle,'-dpng', filename, '-r100');
drawnow
set(ip.Results.figureHandle,'Units',oldscreenunits,...
     'PaperUnits',oldpaperunits,...
     'PaperPosition',oldpaperpos)

end

