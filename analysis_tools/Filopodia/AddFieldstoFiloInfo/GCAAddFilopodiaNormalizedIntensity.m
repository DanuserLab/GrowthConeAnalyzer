function [ filoInfo, avgVeilIntensity ] = GCAAddFilopodiaNormalizedIntensity(img,veilMask,filoInfo,varargin)
% GCAAddFilopodiaNormalizedIntensity: 
% Adds fields to filoInfo.mat corresponding to average filopodia
% intensities (absolute and normalized to the mean veil intensity)
% 
% INPUT: 
% img: rxc matrix corresponding to the image from which to extract the
% intensity values. 
% 
% veilMask: rxc matrix corresponding to the veilMask for a single frame as 
%           output by Step III of GCA. 
%
% filoInfo: rx1 structure output by Step VII of GCA: GCAfitFilopodiaMovie
%           where r is the number of filopodia segments detected for a
%           given frame 
%% Check Input 
ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('PSFSigma',0.43,@(x) isnumeric(x)) ;
ip.CaseSensitive = false;
ip.parse(varargin{:});
%% Calculate the veil intensity without internal filopodia as a normalization factor
% for the actin content intensity.
[ny,nx] =size(img);
% filter the image based on the psf of the microscope.
H = fspecial('gaussian',3,ip.Results.PSFSigma); 
imgFilt = imfilter(img,H); % for weighted averaging

if isfield(filoInfo,'Int_endpointCoordFitPix');
    pixIndicesInt = cell(length(filoInfo),1);
    % get the internal filo to subtract out from the veil normalization
    % intensity
    for i = 1:length(filoInfo)
        pixIndices = filoInfo(i).('Int_pixIndices');
        idxEnd = find(pixIndices == filoInfo(i).('Int_endpointCoordFitPix')(1,1));
        x = pixIndices(1:idxEnd);
        if ~isempty(x) % having trouble with concatenation with the empty cells - in future won't allow
            
            pixIndicesInt{i,1} = pixIndices(1:idxEnd);
        end
        %[yC,xC] = ind2sub( imgSize  ,pixIndicesPlot);
        % plot(xC,yC,'color',colorC,'Linewidth',3);
    end
    maskIntFilo = zeros(size(veilMask));
    maskIntFilo(vertcat(pixIndicesInt{:})) = 1;
    maskIntFilo = imdilate(maskIntFilo,(strel('disk',2)));
    
    veilMask = (veilMask-maskIntFilo);
end % isfield

veilIntensityValues = imgFilt(logical(veilMask));
avgVeilIntensity = mean(veilIntensityValues(:)); % normalization factor
%%
sanityCheck = 0;
if sanityCheck == 1
    setFigure(nx,ny,'on');
    imshow(imgFilt,[]);
    hold on
    spy(~veilMask,'g');
    text(5,5,['Mean Fluorescence Veil/Stem ' num2str(avgVeilIntensity,4) 'AU'],'color','y');
end
%%
% run through all the filopodia collect the intensity to the
% endpoint- calc an normInt Ext, normInt Int, normInt Total.
for iFilo = 1:length(filoInfo)
    
    % actin content of non-embedded filopodia first.
    pixIndices = filoInfo(iFilo).('Ext_pixIndices');
    idxEnd = find(pixIndices == filoInfo(iFilo).Ext_endpointCoordFitPix);
    pixIndicesExt = pixIndices(1:idxEnd);
    
    % Extra-veil filo
    avgFiloIntensityExt =  mean(imgFilt(pixIndicesExt));
    normIntensityExt = (avgFiloIntensityExt/avgVeilIntensity);
    filoInfo(iFilo).Ext_IntensityNormToVeil = normIntensityExt;
    filoInfo(iFilo).Ext_Intensity = avgFiloIntensityExt;
    
    if isfield(filoInfo,'Int_endpointCoordFitPix');
        % Intra-veil filo
        pixIndices = filoInfo(iFilo).('Int_pixIndices');
        idxEnd = find(pixIndices == filoInfo(iFilo).('Int_endpointCoordFitPix')(1,1));
        pixIndicesInt = pixIndices(1:idxEnd);
        
        avgFiloIntensityInt = mean(imgFilt(pixIndicesInt));
        normIntensityInt = (avgFiloIntensityInt/avgVeilIntensity);
        
        filoInfo(iFilo).Int_IntensityNormToVeil = normIntensityInt;
        filoInfo(iFilo).Int_Intensity = avgFiloIntensityInt;
        
        % Total Actin Bundle
        avgFiloIntensityTot = mean(imgFilt([pixIndicesExt;pixIndicesInt]));
        filoInfo(iFilo).Tot_IntensityNormToVeil = avgFiloIntensityTot/avgVeilIntensity;
        filoInfo(iFilo).Tot_Intensity = avgFiloIntensityTot;
    end % isfield
end % iFilo
end % function