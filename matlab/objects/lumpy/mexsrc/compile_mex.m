% ---------------------------------------------------------------------------
% **********************CGRI SIMULATION TOOLBOX***************************
% File: compile_lumpy.m
% Purpose: Matlab script to compile the lumpy background mex files
% Notes: 
% Author:  Nick Henscheid
% Date:    9-2016
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without technical
% support, and with no warranty, express or implied, as to its usefulness for
% any purpose.
% ---------------------------------------------------------------------------
function compile_mex(mexopts)
disp('Compiling Lumpy Object library');
GPU     = mexopts.GPU;
CGRISRC = mexopts.CGRISRC
mexext  = mexopts.mexext;
verb    = mexopts.VERBOSE;
if(GPU)
    CUDALIB = mexopts.CUDALIB;
    if(verb)
        GPUFLAGS = sprintf('-v -L"%s" -lcudart -I"%s"',CUDALIB,CGRISRC);
    else
        GPUFLAGS = sprintf('-L"%s" -lcudart -I"%s"',CUDALIB,CGRISRC);
    end
        
end

if(verb)
    CPUFLAGS = sprintf('-v -I"%s"',CGRISRC);
else
    CPUFLAGS = sprintf('-I"%s"',CGRISRC);
end

%!!NOTE!!: Must leave a space at the beginning of each file name!

if(~exist('../mexbin','dir'))
    mkdir('../mexbin')
end

switch GPU
    case 1
        % Compile binaries for CPU/GPU
        disp('Compiling for CPU/GPU');
        str = [GPUFLAGS,' lumpy_mex_gpu.cu'];
        args = regexp(str,'\s+','split');
        mexcuda(args{:})    % for R2016 and above
        movefile(strcat('lumpy_mex_gpu',mexext),strcat('../mexbin/lumpy_mex_gpu',mexext));
        str = [CPUFLAGS,' lumpy_mex_cpu.cc'];
        args = regexp(str,'\s+','split');
        mex(args{:})
        movefile(strcat('lumpy_mex_cpu',mexext),strcat('../mexbin/lumpy_mex_cpu',mexext));
        cd ../
        addpath(strcat(pwd,'/mexbin'));
    case 0
        % Compile binaries for CPU only 
        disp('Compiling for CPU only');
        str = [CPUFLAGS,' lumpy_mex_cpu.cc'];
        args = regexp(str,'\s+','split');
        mex(args{:})
        movefile(strcat('lumpy_mex_cpu',mexext),strcat('../mexbin/lumpy_mex_cpu',mexext));
        cd ../
        addpath(strcat(pwd,'/mexbin'));
end