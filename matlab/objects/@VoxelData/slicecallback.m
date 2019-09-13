function slicecallback(obj,source,event,panel,handles)
            x0 = handles.slx.Value;set(handles.xtxt,'String',['x= ',num2str(x0)]);
            y0 = handles.sly.Value;set(handles.ytxt,'String',['y= ',num2str(y0)]);
            z0 = handles.slz.Value;set(handles.ztxt,'String',['z= ',num2str(z0)]);
            nx = handles.slnx.Value;set(handles.nxtxt,'String',['nx= ',num2str(nx)]);
            ny = handles.slny.Value;set(handles.nytxt,'String',['ny= ',num2str(ny)]);
            nz = handles.slnz.Value;set(handles.nztxt,'String',['nz= ',num2str(nz)]);
            
            l  = norm([nx,ny,nz]);
            nplot = str2double(handles.nplot.String);
            cmin = str2double(handles.cmin.String);
            cmax = str2double(handles.cmax.String);
            [az,el] = view;
            ax = gca;
            n = length(ax.Children);
            i = 1;
            while i<=n
                if(isa(ax.Children(i),'matlab.graphics.primitive.Surface')||isa(ax.Children(i),'matlab.graphics.chart.primitive.Surface')||isa(ax.Children(i),'matlab.graphics.chart.primitive.Quiver'))
                    delete(ax.Children(i));
                    i = 1;
                    n = n-1;
                else
                    i = i+1;
                end
            end
            if(nz~=0)
                fprintf('Rendering slice...')
                [xsurf,ysurf] = meshgrid(linspace(0,1,nplot));
                zsurf = z0 - (nx/nz)*(xsurf-x0) - (ny/nz)*(ysurf-y0);
                hold on;
                h = slice(obj.X(:,:,:,1),obj.X(:,:,:,2),obj.X(:,:,:,3),obj.Y,xsurf,ysurf,zsurf,'cubic');              
                view(az,el);
            elseif(nx~=0)
                fprintf('Rendering slice...')
                [ysurf,zsurf] = meshgrid(linspace(0,1,nplot));
                xsurf = x0 - (ny/nx)*(ysurf-y0) - (nz/nx)*(zsurf-z0);    
                hold on;
                h = slice(obj.X(:,:,:,1),obj.X(:,:,:,2),obj.X(:,:,:,3),obj.Y,xsurf,ysurf,zsurf,'cubic');              
                view(az,el);
            elseif(ny~=0)
                fprintf('Rendering slice...')
                [xsurf,zsurf] = meshgrid(linspace(0,1,nplot));
                ysurf = y0 - (nx/ny)*(xsurf-x0) - (nz/ny)*(zsurf-z0);
                hold on;
                h = slice(obj.X(:,:,:,1),obj.X(:,:,:,2),obj.X(:,:,:,3),obj.Y,xsurf,ysurf,zsurf,'cubic');    
                view(az,el);
            else     
                warning('Zero normal vector detected.  Using [0,0,1].')
                fprintf('Rendering slice...')
                [xsurf,ysurf] = meshgrid(linspace(0,1,nplot));
                zsurf = z0*ones(256);
                hold on;
                h = slice(obj.X(:,:,:,1),obj.X(:,:,:,2),obj.X(:,:,:,3),obj.Y,xsurf,ysurf,zsurf,'cubic');     
                view(az,el);
            end
            fprintf('done!\n')
            xlim([0,1]);ylim([0,1]);zlim([0,1]);  
            caxis([cmin,cmax]);
            %plot3(x0,y0,z0,'k.','MarkerSize',12);
            quiver3(x0,y0,z0,0.2*nx/l,0.2*ny/l,0.2*nz/l,'LineWidth',1.5,'Marker','.','MarkerSize',12,'Color','k');
            hold off;
           
           set(gca,'Parent',panel);
           set(h,'EdgeColor','none');      
end
