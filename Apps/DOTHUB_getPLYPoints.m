classdef DOTHUB_getPLYPoints < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                    matlab.ui.Figure
        GridLayout                  matlab.ui.container.GridLayout
        LeftPanel                   matlab.ui.container.Panel
        UITable                     matlab.ui.control.Table
        RightPanel                  matlab.ui.container.Panel
        UIAxes                      matlab.ui.control.UIAxes
        LoadPLYButton               matlab.ui.control.Button
        CollectselectedpointButton  matlab.ui.control.Button
        AligntolandmarksButton      matlab.ui.control.Button
        SaveResultsButton           matlab.ui.control.Button
        LoadlabelslistButton        matlab.ui.control.Button
        LoadpointsButton            matlab.ui.control.Button
        ScalepositionstotileedgelengthButton  matlab.ui.control.Button
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    
    properties (Access = private)
        selectedCell = []; %
        mesh
        landmarks
        positions
        plyFileName
        landmarksLineHandle
        positionsLineHandle
        positionLabelHandles
    end
    
    methods (Access = private)
        
        function hm = UIPlotMesh(app)
            axis(app.UIAxes);
            hold(app.UIAxes,'off');
            if ~isempty(app.mesh.node)
                hm = trimesh(app.mesh.face,app.mesh.node(:,1),app.mesh.node(:,2),app.mesh.node(:,3),'Parent',app.UIAxes);
                hm.EdgeColor = 'none';
                hm.FaceColor = 'flat';
                hm.FaceVertexCData = app.mesh.CData;
                %Also plot points
                UIPlotPoints(app)
            else
                hm = [];
            end
        end
        
        function UIPlotPoints(app)
            axis(app.UIAxes);
            if isgraphics(app.landmarksLineHandle)
                delete(app.landmarksLineHandle) %Clear points before replot
            end
            if ~isempty(app.landmarks)
                hold(app.UIAxes,'on');
                app.landmarksLineHandle = plot3(app.UIAxes,app.landmarks(:,1),app.landmarks(:,2),app.landmarks(:,3),'g.','MarkerSize',30);
            end
            
            if isgraphics(app.positionsLineHandle)
                delete(app.positionsLineHandle) %Clear points before replot
            end
            if isgraphics(app.positionLabelHandles)
                delete(app.positionLabelHandles)
            end
            if ~isempty(app.positions)
                hold(app.UIAxes,'on');
                app.positionsLineHandle = plot3(app.UIAxes,app.positions(:,1),app.positions(:,2),app.positions(:,3),'r.','MarkerSize',30);
                for ii = 1:size(app.positions,1)
                    app.positionLabelHandles(ii) = text(app.UIAxes,app.positions(ii,1)*1.1,app.positions(ii,2)*1.1,app.positions(ii,3)*1.1,app.UITable.Data{ii+5,1},'color','r');
                end
            end
            axis(app.UIAxes,'equal');
            hold(app.UIAxes,'off');
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.UITable.Data = cell(100,4);
            app.UITable.ColumnName = {'Label';'x';'y';'z'};
            app.UITable.Data(1,1:4) = {'Nasion',0,0,0};
            app.UITable.Data(2,1:4) = {'Inion',0,0,0};
            app.UITable.Data(3,1:4) = {'Ar',0,0,0};
            app.UITable.Data(4,1:4) = {'Al',0,0,0};
            app.UITable.Data(5,1:4) = {'Cz',0,0,0};
            app.UIAxes.DataAspectRatio = [1 1 1];
            %For debugging###:
            app.mesh.node = [];
            app.mesh.face = [];
            app.mesh.CData = [];
            % #####
            app.landmarks = [];
            app.positions = [];
            app.plyFileName = '';
            app.landmarksLineHandle = [];
            app.positionsLineHandle = [];
            hm = UIPlotMesh(app);
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {558, 558};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {256, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end

        % Cell edit callback: UITable
        function UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            
            %Make landmark labels uneditable
            if indices(1)<=5 && indices(2)==1
                tmp = {'Nasion' 'Inion' 'Ar' 'Al' 'Cz'};
                app.UITable.Data{indices(1),indices(2)} = tmp(indices(1));
                return
            else
                app.UITable.Data{indices(1),indices(1)} = newData;
            end       
        end

        % Button pushed function: LoadPLYButton
        function LoadPLYButtonPushed(app, event)
            [fname, pname] = uigetfile('*.ply','Selet .PLY or .mshs to load');
            if fname==0
                return
            end
            app.plyFileName = fullfile(pname,fname);
            app.mesh = DOTHUB_ply2Mesh(app.plyFileName);
            hm = UIPlotMesh(app);
        end

        % Button pushed function: CollectselectedpointButton
        function CollectselectedpointButtonPushed(app, event)
            dcm = datacursormode(app.UIFigure);
            dcm.Enable = 'on';
            info = getCursorInfo(dcm);
            %Determine currently selected row
            currCellInd = app.selectedCell;
            if ~isempty(currCellInd) && isstruct(info)
                app.UITable.Data(currCellInd(1),2) = {info.Position(1)};
                app.UITable.Data(currCellInd(1),3) = {info.Position(2)};
                app.UITable.Data(currCellInd(1),4) = {info.Position(3)};
                
                if currCellInd(1)<=5 %Landmark update
                    app.landmarks(currCellInd(1),1) = info.Position(1);
                    app.landmarks(currCellInd(1),2) = info.Position(2);
                    app.landmarks(currCellInd(1),3) = info.Position(3);
                else %position update
                    app.positions(currCellInd(1)-5,1) = info.Position(1);
                    app.positions(currCellInd(1)-5,2) = info.Position(2);
                    app.positions(currCellInd(1)-5,3) = info.Position(3);
                end
                UIPlotPoints(app);
            else
                return
            end
            
        end

        % Cell selection callback: UITable
        function UITableCellSelection(app, event)
            indices = event.Indices;
            app.selectedCell = event.Indices;
        end

        % Button pushed function: SaveResultsButton
        function SaveResultsButtonPushed(app, event)
            %Get data from table
            tmp = app.UITable.Data;
            popInd = ~cellfun(@isempty,tmp(:,1));
            tmpData = tmp(popInd,1:4);
            
            %Check scaling
            answer = questdlg('Do you want to scale these measurements (LUMO TILES)?','Scale?','Yes','No','Yes');
            if strcmpi(answer,'Yes') %Determine mean euc distance between points on the same tile, and scale
                dataTmp = cell2mat(app.UITable.Data(6:end,2:4));
                if mod(length(dataTmp),3)~=0
                    error('Number of positions is not a multiple of 3...')
                end
                nTile = length(dataTmp)/3;
                count = 1;
                for ii = 1:nTile
                    A = dataTmp(1+(ii-1)*3,:);
                    B = dataTmp(2+(ii-1)*3,:);
                    C = dataTmp(3+(ii-1)*3,:);
                    triDist(count) = sqrt(sum((A-B).^2));
                    count = count+1;
                    triDist(count) = sqrt(sum((A-C).^2));
                    count = count+1;        
                    triDist(count) = sqrt(sum((B-C).^2));
                    count = count+1;                          
                end
                %Determine scaling factor;
                scaleFact = 18/mean(triDist);
                figure;
                hist(triDist*scaleFact);
                xlabel('Extracted triangle edge distance (should be 18) (mm)');
                ylabel('Count');
                
                for ii = 1:size(tmpData,1)
                    for jj = 2:4
                        tmpData{ii,jj} = tmpData{ii,jj}.*scaleFact;
                    end
                end
                    
            end
            if ~isempty(app.plyFileName)
                [pathtmp,nametmp,~] = fileparts(app.plyFileName);
                [file,path,~] = uiputfile([pathtmp '/' nametmp '_Positions.csv']);
            else
                [file,path,~] = uiputfile([pwd '/' 'Positions.csv']);
            end
           
            writecell(tmpData,fullfile(path,file));
        end

        % Button pushed function: AligntolandmarksButton
        function AligntolandmarksButtonPushed(app, event)
            
            landmarks = cell2mat(app.UITable.Data(1:5,2:4));
            meshPoints = app.mesh.node(:,1:3);
            
            %Shift so inion is origin
            Iz = landmarks(2,:);
            landmarks = landmarks - repmat(Iz,size(landmarks,1),1);
            meshPoints = meshPoints - repmat(Iz,size(meshPoints,1),1);
                        
            %Rotate Cz to be in +ve z;
            Cz = landmarks(5,:);            
            [azimuth,elevation,~] = cart2sph(Cz(1),Cz(2),Cz(3));
            landmarks = rotz(landmarks,90-rad2deg(azimuth),[0 0 0]);
            meshPoints = rotz(meshPoints,90-rad2deg(azimuth),[0 0 0]);
            landmarks = rotx(landmarks,-rad2deg(elevation),[0 0 0]);
            meshPoints = rotx(meshPoints,-rad2deg(elevation),[0 0 0]);
            
            %Rotate so Nz is at (0 y 0)
            Nz = landmarks(1,:);
            Iz = landmarks(2,:);
            Cz = landmarks(5,:);
            [azimuth,elevation,~] = cart2sph(Nz(1)-Iz(1),Nz(2)-Iz(2),Nz(3)-Iz(3));
            landmarks = rotz(landmarks,90-rad2deg(azimuth),[0 0 0]);
            meshPoints = rotz(meshPoints,90-rad2deg(azimuth),[0 0 0]);
            landmarks = rotx(landmarks,-rad2deg(elevation),[0 0 0]);
            meshPoints = rotx(meshPoints,-rad2deg(elevation),[0 0 0]);
            
            %Rotate so Ar and Al are perpendicular to z y plane
            Ar = landmarks(3,:);
            Al = landmarks(4,:);
            [~,elevation,~] = cart2sph(Ar(1)-Al(1),Ar(2)-Al(2),Ar(3)-Al(3));
            landmarks = roty(landmarks,rad2deg(elevation),[0 0 0]);
            meshPoints = roty(meshPoints,rad2deg(elevation),[0 0 0]);
            
            %Reset landmarks in table
            app.UITable.Data(1:5,2:4) = num2cell(landmarks);
            
            %Reset landmarks in appdata
            app.landmarks = landmarks;
            
            %Replot mesh
            app.mesh.node = meshPoints;
            hm = UIPlotMesh(app);
 
        end

        % Button pushed function: LoadlabelslistButton
        function LoadlabelslistButtonPushed(app, event)
            [fname, pname] = uigetfile({'*.txt';'*.csv'},'Selet label list to load');
            if fname==0
                return
            end
            inputList = readtable(fullfile(pname,fname),'ReadVariableNames',0);
            inputList = table2cell(inputList(:,1));
            inputList = inputList(~cellfun(@isempty,inputList));
            app.UITable.Data(6:end,:) = [];
            app.UITable.Data(6:6+length(inputList)-1,1) = inputList;
            %Replot points
            UIPlotPoints(app);
        end

        % Button pushed function: LoadpointsButton
        function LoadpointsButtonPushed(app, event)
            [fname, pname] = uigetfile({'*.csv'},'Selet points to load');
            if fname==0
                return
            end
            inputList = readtable(fullfile(pname,fname),'ReadVariableNames',0);
            inputList = table2cell(inputList(:,1:4));
            app.UITable.Data = inputList;
            %Replot points
            UIPlotPoints(app);
        end

        % Button pushed function: 
        % ScalepositionstotileedgelengthButton
        function ScalepositionstotileedgelengthButtonPushed(app, event)
            tmp = app.UITable.Data;
            popInd = ~cellfun(@isempty,tmp(:,1));
            tmpData = tmp(popInd,1:4);
            
            dataTmp = cell2mat(app.UITable.Data(6:end,2:4));
            if mod(length(dataTmp),3)~=0
                error('Number of positions is not a multiple of 3...')
            end
            nTile = length(dataTmp)/3;
            count = 1;
            for ii = 1:nTile
                A = dataTmp(1+(ii-1)*3,:);
                B = dataTmp(2+(ii-1)*3,:);
                C = dataTmp(3+(ii-1)*3,:);
                triDist(count) = sqrt(sum((A-B).^2));
                count = count+1;
                triDist(count) = sqrt(sum((A-C).^2));
                count = count+1;
                triDist(count) = sqrt(sum((B-C).^2));
                count = count+1;
            end
            %Determine scaling factor;
            scaleFact = 18/mean(triDist);
            f1 = figure;
            hist(triDist*scaleFact);
            xlabel('Extracted triangle edge distance (should be 18) (mm)');
            ylabel('Count');
            
            for ii = 1:size(tmpData,1)
                for jj = 2:4
                    tmpData{ii,jj} = tmpData{ii,jj}.*scaleFact;
                end
            end
            
            answer = questdlg('DO you want to proceed with scaling?','Scale?','Yes','No','Yes');
            close(f1)
            if strcmpi(answer,'Yes')
                for ii = 1:size(tmpData,1)
                    for jj = 2:4
                        app.UITable.Data{ii,jj} = tmpData{ii,jj};
                    end
                end
                app.landmarks = cell2mat(app.UITable.Data(1:5,2:4));
                app.positions = cell2mat(app.UITable.Data(6:end,2:4));
                app.mesh.node(:,1:3) = app.mesh.node(:,1:3).*scaleFact;
                UIPlotMesh(app);
            end
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 883 558];
            app.UIFigure.Name = 'UI Figure';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {256, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create UITable
            app.UITable = uitable(app.LeftPanel);
            app.UITable.ColumnName = {'Label'; 'x'; 'y'; 'z'};
            app.UITable.ColumnWidth = {80, 40, 40, 40};
            app.UITable.ColumnEditable = [true false false false];
            app.UITable.CellEditCallback = createCallbackFcn(app, @UITableCellEdit, true);
            app.UITable.CellSelectionCallback = createCallbackFcn(app, @UITableCellSelection, true);
            app.UITable.Position = [6 19 244 533];

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create UIAxes
            app.UIAxes = uiaxes(app.RightPanel);
            title(app.UIAxes, 'PLY Mesh')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            app.UIAxes.Position = [15 19 592 394];

            % Create LoadPLYButton
            app.LoadPLYButton = uibutton(app.RightPanel, 'push');
            app.LoadPLYButton.ButtonPushedFcn = createCallbackFcn(app, @LoadPLYButtonPushed, true);
            app.LoadPLYButton.Position = [15 509 119 40];
            app.LoadPLYButton.Text = 'Load PLY';

            % Create CollectselectedpointButton
            app.CollectselectedpointButton = uibutton(app.RightPanel, 'push');
            app.CollectselectedpointButton.ButtonPushedFcn = createCallbackFcn(app, @CollectselectedpointButtonPushed, true);
            app.CollectselectedpointButton.Position = [184 461 127 88];
            app.CollectselectedpointButton.Text = {'Collect selected'; 'point'; ''};

            % Create AligntolandmarksButton
            app.AligntolandmarksButton = uibutton(app.RightPanel, 'push');
            app.AligntolandmarksButton.ButtonPushedFcn = createCallbackFcn(app, @AligntolandmarksButtonPushed, true);
            app.AligntolandmarksButton.Position = [321 509 127 40];
            app.AligntolandmarksButton.Text = {'Align to'; 'landmarks'};

            % Create SaveResultsButton
            app.SaveResultsButton = uibutton(app.RightPanel, 'push');
            app.SaveResultsButton.ButtonPushedFcn = createCallbackFcn(app, @SaveResultsButtonPushed, true);
            app.SaveResultsButton.Position = [507 513 100 36];
            app.SaveResultsButton.Text = 'Save Results';

            % Create LoadlabelslistButton
            app.LoadlabelslistButton = uibutton(app.RightPanel, 'push');
            app.LoadlabelslistButton.ButtonPushedFcn = createCallbackFcn(app, @LoadlabelslistButtonPushed, true);
            app.LoadlabelslistButton.Position = [15 464 119 40];
            app.LoadlabelslistButton.Text = 'Load labels list';

            % Create LoadpointsButton
            app.LoadpointsButton = uibutton(app.RightPanel, 'push');
            app.LoadpointsButton.ButtonPushedFcn = createCallbackFcn(app, @LoadpointsButtonPushed, true);
            app.LoadpointsButton.Position = [15 418 119 40];
            app.LoadpointsButton.Text = 'Load points';

            % Create ScalepositionstotileedgelengthButton
            app.ScalepositionstotileedgelengthButton = uibutton(app.RightPanel, 'push');
            app.ScalepositionstotileedgelengthButton.ButtonPushedFcn = createCallbackFcn(app, @ScalepositionstotileedgelengthButtonPushed, true);
            app.ScalepositionstotileedgelengthButton.Position = [321 461 127 40];
            app.ScalepositionstotileedgelengthButton.Text = {'Scale positions to '; 'tile edge length'};

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = DOTHUB_getPLYPoints

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end