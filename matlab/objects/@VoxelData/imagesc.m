 function imagesc(obj,varargin)
 
    delete(obj.imfig);  % Delete if already exists
    if(obj.dim == 2)
       nparam   = obj.N(end);
       dx       = (obj.L(1,2)-obj.L(1,1))/obj.N(1);
       dy       = (obj.L(2,2)-obj.L(2,1))/obj.N(2);
       Xplot{1} = (obj.L(1,1)+dx/2):dx:(obj.L(1,2)-dx/2);
       Xplot{2} = (obj.L(2,1)+dy/2):dy:(obj.L(2,2)-dy/2);
       obj.imfig = figure('Position',[500,500,600,600]);
       obj.imfig.CloseRequestFcn = @(src,data) obj.imfigCloseFcn(src,data);
       
       uicontrol('Parent',obj.imfig,'Style','text','FontSize',12,...
           'String','Show Grid:','Position',[20,20,100,20]);
       uicontrol('Parent',obj.imfig,'Style','checkbox',...
           'FontSize',12,'Position',[130,20,20,20],...
           'Callback',@(src,evt)obj.gridCallback(src,evt),...
           'Value',0);
       uicontrol('Parent',obj.imfig,'Style','text','FontSize',12,...
           'String','Parameter:','Position',[170,20,100,20]);
       uicontrol('Parent',obj.imfig,'Style','popupmen',...
           'String',{1:nparam},'FontSize',12,'Position',[280,20,80,20],...
           'Callback',@(src,evt)obj.paramCallback(src,evt),'Value',1);
       
       if(isreal(obj.Y))
           imagesc(Xplot{1},Xplot{2},obj.Y(:,:,1));
           set(gca,'YDir','normal');
           axis image;
       else

           for j=1:obj.datadim 
               T1{j} = sprintf('Real part, layer %i',j);
               T2{j} = sprintf('Imaginary part, layer %i',j);
               T3{j} = sprintf('Modulus, layer %i',j);
           end
           %DisplayImages([1,obj.datadim],real(obj.Y),'title',T1);
           %DisplayImages([1,obj.datadim],imag(obj.Y),'title',T2);
           %DisplayImages([1,obj.datadim],abs(obj.Y),'title',T3);
       end
       title(sprintf('%i-dimensional voxel data with %i parameter(s)',obj.dim,obj.N(end)));
    end
 end
 