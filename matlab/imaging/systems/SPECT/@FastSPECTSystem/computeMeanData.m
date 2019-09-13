function gbar = computeMeanData(obj,f)
    % Computes the mean image data for the system, given the object f
    % \bar{g} = Hf (c-to-d)  or \bar{g} = Lf (listmode/particle-processing)
    
    
    % First, need to decide if system is list-mode/PP or integrating/binned
    % If the system is list-mode, gbar will be an object of type
    % ListModeMeanData, which is effectively a function-type class
    % If the system is integrated/binned, gbar will be an object of type
    % BinnedModeMeanData, which is a vector-type class 
    
    
    
    % Then, need to set up the appropriate evaluation points (r,s)\in
    % \partial_+ \Gamma
    
    
    % Then, call the RTE solver to obtain the mean data 
    


end