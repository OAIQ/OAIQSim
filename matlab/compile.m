% ---------------------------------------------------------------------------
% **********************CGRI SIMULATION TOOLBOX***************************
% File: compile.m
% Purpose: Matlab script to compile the CGRI Simulation Toolbox mex files
% Notes: 
% Author:  Nick Henscheid
% Date:    9-2016
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without technical
% support, and with no warranty, express or implied, as to its usefulness for
% any purpose.
% ---------------------------------------------------------------------------

%Set home directory.
str = {pwd};  
if(~strcmp(str{1}(end-5:end),'matlab'))
    errordlg('Must run script from /matlab directory!');
    return;
end
CGRIMATLAB = str{1};
CGRIHOME = str{1}(1:end-7);  %Cute trick to extract base dir

% Directories
CGRISRC  = strcat(CGRIHOME,'/src');
CUDALIB  = '';

% Verbose mode 
VERBOSE = 0;

% Get architecture
arch=computer;
mac=strcmp(arch,'MACI64') || strcmp(arch,'MACI') || strcmp(arch,'MAC');
linux=strcmp(arch,'GLNXA64') || strcmp(arch,'GLNX86');
pc= strcmp(arch,'PCWIN64');
GPU = sign(gpuDeviceCount); % To compile with gpu or not
if(GPU)
    CUDALIB = '/usr/local/cuda/lib64';
    if(usejava('desktop'))
        choice = questdlg('Nvidia GPU Found!  Compile CUDA files?','Compile with CUDA?','Yes','No','Cancel Compilation','Yes');
        switch choice
            case 'Yes'
                questdlg(['Cuda library currently set to ',CUDALIB,', is this correct?'],'CUDA Library correct?','Yes','No','Cancel Compilation','Yes'); 
                GPU = 1;
            case 'No'
                GPU = 0;
            otherwise
                msgbox('Canceling compilation');
                return;
        end
    else
        GPU = input('Nvidia GPU Found!  Compile CUDA files? (1 = yes, 0 = no) ');
        input(['Compiling with CUDA.  CUDA library currently set to ',CUDALIB,', if this is incorrect modify compile.m line 34! (Press any key to continue)']);
    end
end
if mac
    mexext = '.mexmaci64'; % Assuming 64 bit.
end
if linux
    mexext = '.mexa64'; % Assuming 64 bit.
end
if pc
    mexext = '.win64'; % Assuming 64 bit.
end

MEXOPTS = struct('CGRISRC',CGRISRC,'mexext',mexext,'CUDALIB',CUDALIB,'GPU',GPU,'VERBOSE',VERBOSE);

%folders = {'objects/diffusion'}
%folders = {'objects/lumpy','objects/diffusion','objects/anisolumpy','imaging/rte/RTELumpy'};
folders  = {'objects/lumpy'};
%folders = {'imaging/rte/RTELumpy'};
%folders = {'imaging/rte/RTEFourier'};
for iFolder = 1:length(folders)
    addpath([CGRIMATLAB,'/',folders{iFolder}]);
    cd(folders{iFolder})
    if(exist('./mexsrc','dir')==7)
        disp('mex folder exists!  Compiling...')
        cd mexsrc
        try
            compile_mex(MEXOPTS);
        catch me
            warning(['Something went wrong with mex compilation in foler ',pwd,'!']);
            disp(me.message)
            cd(CGRIMATLAB)
        end
            
    end
    cd(CGRIMATLAB)
end

clear all;