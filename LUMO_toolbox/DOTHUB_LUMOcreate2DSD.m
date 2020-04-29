function SD = DOTHUB_LUMOcreate2DSD(nTiles,outname)

% Quick and dirty gui to create 2D layouts
% Requires entry of tile centres, and angles.  Grid arrangement is fixed,
% but can be rotated to the preferred orientation.

% Create grid with x-spacing of tiles of 35 mm, with offset rows
% and y-spacing of 30.31 mm.

if ~exist('nTiles','var')
    tmp = inputdlg('Enter number of LUMO tiles...','Tile number');
    nTiles = str2num(tmp{1});
end
xspace = 35;
yspace = 30.31;
ng = 6;

%buildgrid
xlin = -ng*xspace:xspace:+ng*xspace;
ylin = -ng*yspace:yspace:+ng*yspace;
[tmp1,tmp2] = meshgrid(xlin,ylin);
tmp1(2:2:end,:) = tmp1(2:2:end,:)+xspace/2;%Offset alternate rows
grid(:,1) = tmp1(:);
grid(:,2) = tmp2(:);
grid(:,3) = 0;

centers_zero = zeros(nTiles,3);
centers_zero(:,1) = 0:xspace:(nTiles-1)*xspace;

if ~exist('initialCent_Orient','var')
    centers = centers_zero;
    orientations = zeros(nTiles,1);
else
    centers(:,1:2) = initialCent_Orient(:,1:2);
    centers(:,3) = 0;
    centers_zero = centers;
    orientations = initialCent_Orient(:,3);
end


%% Now enter editable display;
ss = get(0,'screensize');ss = ss(3:4);
fig = uifigure('Position',[0 0 ss(1) ss(2)]);
set(fig,'Name','Define 2D layout');
uit = uitable(fig);
uit.Data = zeros(nTiles,3);
uit.Data = [centers(:,1:2) orientations];
uit.ColumnName = {'X','Y','Angle'};
uit.ColumnWidth = {60,60,60};
uit.ColumnEditable = true;
uit.CellEditCallback = @updatePlot;
uit.Position = [0.05*ss(1) 0.05*ss(2) 0.15*ss(1) 0.75*ss(2)];
ax = uiaxes(fig);
ax.Position = [0.25*ss(1) 0.25*ss(2) 0.7*ss(1) 0.60*ss(2)];
bt = uibutton(fig,'Text','FINISHED');
bt.Position = [0.85*ss(1) 0.05*ss(2) 0.1*ss(1) 0.075*ss(2)];
bt.ButtonPushedFcn = @SaveandQuit;

bt2 = uibutton(fig,'Text','Rotate Grid & re-fit');
bt2.Position = [0.25*ss(1) 0.05*ss(2) 0.1*ss(1) 0.075*ss(2)];
bt2.ButtonPushedFcn = @RotateGrid;

[source_array, detector_array] = getTileOptodePos_Callback;


%% Nested functions

    function [nnode, ind, offset] = Nearest_node_RJC(point,nodes)
        
        %This function outputs the node (1x3) present in the list nodes (Mx3) that is closest
        %to the 3D location specified by point (1x3) and the offset in whatever
        %dimensions are provided.
        
        %Inputs (point,nodes) : the single point and the list of nodes
        %Outputs [nnode, offset]: the index of the specific node in the list
        %'nodes' which is nearest to the input point, and the offset (euclidean
        %error).
        
        if size(point,2)~=3
            point = point';
            if size(point,2)~=3
                error('Dimensions of input point not 1x3');
            end
        end
        
        if size(nodes,2)~=3;
            nodes = nodes';
            if size(nodes,2)~=3;
                error('Dimensions of input nodes not Mx3');
            end
        end
        
        M = size(nodes,1);
        
        if size(point,1)>1
            for i = 1:size(point,1)
                tmpPoint = point(i,:);
                dists = sqrt(sum((nodes - repmat(tmpPoint,M,1)).^2,2));
                [offset(i),ind(i)] = min(dists);
            end
        else
            dists = sqrt(sum((nodes - repmat(point,M,1)).^2,2));
            [offset,ind] = min(dists);
        end
        
        nnode = nodes(ind,:);
    end

    function updatePlot(ssrc,event)
        tmp = uit.Data;
        centers(:,1:2) = tmp(:,1:2);
        orientations = tmp(:,3);
        %Force to grid
        [centers] = Nearest_node_RJC(centers,grid);
        orientations = round((orientations)/30)*30;
        
        uit.Data(:,1:2) = centers(:,1:2);
        uit.Data(:,3) = orientations;
        getTileOptodePos_Callback
    end

    function [source_array, detector_array] = getTileOptodePos_Callback
        
        %DO NOT EDIT THIS #########################################################
        sources(1,:) = [0 10.81 0];
        sources(2,:) = [-10.81*cosd(30) -10.81*sind(30) 0];
        sources(3,:) = [10.81*cosd(30) -10.81*sind(30) 0];
        detectors(1,:) = [0 -8.88 0];
        detectors(2,:) = [-8.88*cosd(30) 8.88*sind(30) 0];
        detectors(3,:) = [0 0 0];
        detectors(4,:) = [8.88*cosd(30) 8.88*sind(30) 0];
        % #########################################################
        
        source_array = [];
        detector_array = [];
        for i = 1:size(centers,1)
            srcTmp = sources;
            detTmp = detectors;
            
            srcTmpRot = rotz(srcTmp,orientations(i)) + repmat(centers(i,:),3,1);
            detTmpRot = rotz(detTmp,orientations(i)) + repmat(centers(i,:),4,1);
            
            source_array = [source_array;srcTmpRot];
            detector_array = [detector_array;detTmpRot];
        end
        
        cla(ax);
        hold(ax,'on');
        
        %plot grid
        plot(ax,grid(:,1),grid(:,2),'k.','markersize',10);hold on;
        %plot array
        plot(ax,source_array(:,1),source_array(:,2),'r.','markersize',30);
        plot(ax,detector_array(:,1),detector_array(:,2),'b.','markersize',30);
        for i = 1:size(source_array,1)
            text(ax,source_array(i,1)+3,source_array(i,2),['S' num2str(i)],'color','r','FontSize',8);
        end
        for i = 1:size(detector_array,1)
            text(ax,detector_array(i,1)+3,detector_array(i,2),['D' num2str(i)],'color','b','FontSize',8);
        end
        linend = source_array(1:3:end,[1 2]);
        linstart = detector_array(3:4:end,[1 2]);
        for i = 1:length(linend);
            line(ax,[linstart(i,1) linend(i,1)],[linstart(i,2) linend(i,2)],'color','g','LineWidth',2);hold on;
            plot(ax,linstart(i,1),linstart(i,2),'g.','MarkerSize',20);
        end
        
        tmp = [source_array(:,1); detector_array(:,1)];
        xl = [min(tmp(:))-0.1*abs(min(tmp(:))) max(tmp(:))+0.1*abs(max(tmp(:)))];
        ax.XLim = xl;
        tmp = [source_array(:,2); detector_array(:,2)];
        yl = [min(tmp(:))-0.1*abs(min(tmp(:))) max(tmp(:))+0.1*abs(max(tmp(:)))];
        ax.YLim = yl;
        ax.DataAspectRatio = [1 1 1];
        
        figHandles = findobj('Type', 'figure');
        for i = 1:length(figHandles)
            if isempty(figHandles(i).Name)
                close(figHandles(i));
            end
        end
        %remove erroneous figures generated I don't know where...
        
    end

    function RotateGrid(~,~)
        grid = rotz(grid,90,[0 0 0]);
        orientations = round((orientations)/30)*30 + 90;
        orientations(orientations>360) = orientations(orientations>360)-360;
        %Force to grid
        [centers] = Nearest_node_RJC(centers_zero,grid);
        
        uit.Data(:,1:2) = centers(:,1:2);
        uit.Data(:,3) = orientations;
        
        getTileOptodePos_Callback;
    end

    function SaveandQuit(~,~)
        
        SD_2Dtmp.SrcPos = source_array;
        SD_2Dtmp.DetPos = detector_array;
        SD_2Dtmp.nSrcs = size(source_array,1);
        SD_2Dtmp.nDets = size(detector_array,1);
        SD_2Dtmp.Lambda = [735 850];
        SD_2Dtmp.SpatialUnit = 'mm';
        
        count = 1;
        for s = 1:SD_2Dtmp.nSrcs
            for d = 1:SD_2Dtmp.nDets
                SD_2Dtmp.MeasList(count,:) = [s d 0 1];
            end
        end
        SD_2Dtmp.MeasList = [SD_2Dtmp.MeasList; SD_2Dtmp.MeasList];
        SD_2Dtmp.MeasList(end/2+1:end,4) = 2;
             
        SD = SD_2Dtmp;
        
        if ~exist('outname','var')
            tmp = inputdlg('Enter save name...','Filename');
            outname = tmp{1};
        end
        
        if ~strcmpi(outname(end-2:end),'.SD')
            outname = [outname '.SD'];
        end
        
        save(outname,'SD');
        close all;
        delete(fig);
        
    end

    function out = rotz(in,ang,about)
        
        %in = Nx3 points in 3D
        %ang = angle in degrees
        %about = centre point of rotation (1x3)
        R = [cosd(ang) -sind(ang) 0; sind(ang) cosd(ang) 0; 0 0 1];
        
        if ~exist('about','var')
            out = (R*in')';
        else
            tmpIn = in - repmat(about,size(in,1),1);
            tmpRot = (R*tmpIn')';
            out = tmpRot + repmat(about,size(in,1),1);
        end
    end

    function out = rotx(in,ang,about)
        
        %in = Nx3 points in 3D
        %ang = angle in degrees
        %about = centre point of rotation (1x3)
        R = [1 0 0;0 cosd(ang) -sind(ang); 0 sind(ang) cosd(ang)];
        
        if ~exist('about','var')
            out = (R*in')';
        else
            tmpIn = in - repmat(about,size(in,1),1);
            tmpRot = (R*tmpIn')';
            out = tmpRot + repmat(about,size(in,1),1);
        end
    end
end
