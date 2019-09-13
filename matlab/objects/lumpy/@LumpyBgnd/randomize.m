function randomize(obj)
    %  Should start calling the object's support function unifRand function
    %  instead - the support set might not be rectangular 
    L = obj.L;  % Get bounding box;  L = [a,b;c,d] for [a,b]x[c,d] or L = [a,b;c,d;e,f] for 
    %                                                                      [a,b]x[c,d]x[e,f].
    rng('shuffle');
    padx = obj.pad_factor*sqrt(max(obj.cov(1,:))/2);  %padding to improve stationarity
    pady = obj.pad_factor*sqrt(max(obj.cov(2,:))/2);
   
    minx = L(1,1)-padx;  
    maxx = L(1,2)+padx;
    miny = L(2,1)-pady;
    maxy = L(2,2)+pady;
    if(obj.dim==3)
        padz = obj.pad_factor*sqrt(max(obj.cov(3,:))/2);
        minz = L(3,1)-padz;
        maxz = L(3,2)+padz;
    end
    if(obj.random_number_of_lumps)
        k = poissrnd(obj.Kbar);
    else 
        k = obj.Kbar;
    end
    
    obj.centers = zeros(k,obj.dim);
    obj.centers(:,1) = minx + (maxx-minx).*rand(k,1);
    obj.centers(:,2) = miny + (maxy-miny).*rand(k,1);
    if(obj.dim==3)
        obj.centers(:,3) = minz + (maxz-minz).*rand(k,1);
    end
end