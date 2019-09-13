function start(SavedConfig)
	% Entry to the cgrisim GUI system
	% Pass an object of type CGRISimGui to open with a previously saved
	% configuration 
    
    Version = 1.0;
    
    fprintf('Welcome to OAIQSim Version %1.2f\n',Version);
	
	close all;
	%addpath([pwd,'/gui'],[pwd,'/util'],[pwd,'/imaging'],[pwd,'/objects'],[pwd,'/observers'],[pwd,'/recon'])
	addpath(genpath([pwd,'/gui']),genpath([pwd,'/imaging']),genpath([pwd,'/util']),...
            genpath([pwd,'/objects']),genpath([pwd,'/test']),genpath([pwd,'/geom']));
    
    
    % Set default figure formatting
    set(0,'defaultAxesFontSize',16)
    set(groot,'defaulttextinterpreter','latex');  
    set(groot, 'defaultAxesTickLabelInterpreter','latex');  
    set(groot, 'defaultLegendInterpreter','latex');  
    
    
	%if(nargin>0 && isa(SavedConfig,'CGRISimGui'))
		
	%else
	%	CGRISimGui;
	%end

end