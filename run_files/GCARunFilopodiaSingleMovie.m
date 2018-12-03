function GCARunFilopodiaSingleMovie(MD,varargin)
% GCARunFilopodiaSingleMovie
% Run the filopodia workflow 

%% Input 
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('segParameterMFile',[]);

ip.parse(varargin{:});
%% STEP VI: (Purple Block in Figure S2)
   
    if isempty(ip.Results.segParameterMFile)
        % if myParamsName empty just run the defaults internally
        GCAReconstructFilopodiaMovie(MD);
    else % load any parameter changes from your personalized file
        cFunc = str2func([ip.Results.segParameterMFile]);
        p = cFunc('step2load',6);
        GCAReconstructFilopodiaMovie(MD,p);
    end
%% STEP VII (Pink Block in Figure S2)
    
    if isempty(ip.Results.segParameterMFile)
        GCAfitFilopodiaMovie(MD);
    else
        cFunc = str2func([ip.Results.segParameterMFile]);
        p = cFunc('step2load',7);
        GCAfitFilopodiaMovie(MD,p)
    end % 
end % function
    
    