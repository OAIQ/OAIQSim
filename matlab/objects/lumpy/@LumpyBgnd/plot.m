function varargout = plot(obj,varargin)
    % Eventually want to make this a GUI (some pieces already exist but need to
    % be updated...
    
    
    p = inputParser; 
    addParameter(p,'axes',0);
    addParameter(p,'alpha',0.5);
    addParameter(p,'grid',0);
    addParameter(p,'centers',0);
   
    p.parse(varargin{:});
    axes_for_plot        = p.Results.axes;
    alpha                = p.Results.alpha;
    input_params.grid    = p.Results.grid;
    input_params.axes    = p.Results.axes;
    %input_params.centers = p.Results.centers;
       
    %axes(axes_for_plot);

    n = obj.N;

    if(obj.dim == 2)
         plot(obj.eval,input_params);  % Pass through to the VoxelData plotting object
%         X = obj.support.UnifGrid([n,n]);
%         xplot = X(1,:,1);
%         imagesc(xplot,xplot,obj.eval);set(gca,'YDir','normal');axis image;
        if(p.Results.centers)
            hold on;
            plot(obj.centers(:,1),obj.centers(:,2),'k.');
            hold off;
        end
    elseif(obj.dim == 3)
        plot(obj.support,axes_for_plot,alpha)
        fprintf('Plotting 3D Isosurface...this may take a moment...\n');
        X = obj.support.UnifGrid([n,n,n]);
        Xeval = [reshape(X(:,:,:,1),[n^3,1]),reshape(X(:,:,:,2),[n^3,1]),reshape(X(:,:,:,3),[n^3,1])];
        u  = obj.eval(Xeval,[n,n,n]);
        if(nargin<2)
            alpha = mean(u);
        end
        s = isosurface(X(:,:,:,1),X(:,:,:,2),X(:,:,:,3),u.data_array,alpha);
        p = patch(s);
        p.FaceColor = [1,0,0];
        p.EdgeAlpha = 0.2;
        title(sprintf('Isosurface plot for alpha = %f',alpha));
    end
    
    
        
   
    if(nargout==1)
        varargout{1} = gcf;
    elseif(nargout==2)
        varargout{1} = gcf;
        varargout{2} = p;
    end
    
end