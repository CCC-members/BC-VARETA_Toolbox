classdef Dataset < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        DatasetsUIFigure  matlab.ui.Figure
        Toolbar           matlab.ui.container.Toolbar
        LoadPushTool      matlab.ui.container.toolbar.PushTool
        PushTool2         matlab.ui.container.toolbar.PushTool
        DeletePushTool    matlab.ui.container.toolbar.PushTool
        PushTool          matlab.ui.container.toolbar.PushTool
        CloseButton       matlab.ui.control.Button
        UITable           matlab.ui.control.Table
    end

    
    properties (Access = public)
        BCVdir;
        Datasets_file;
        Datasets;
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.UITable.SelectionType = 'row';
        end

        % Close request function: DatasetsUIFigure
        function DatasetsUIFigureCloseRequest(app, event)
            delete(app)            
        end

        % Button pushed function: CloseButton
        function CloseButtonPushed(app, event)
            delete(app);
        end

        % Cell edit callback: UITable
        function UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            uDataset = app.Datasets(indices(indices(1)));
            switch indices(2)
                case 1
                    uDataset.Dataset_name = newData;
                case 2
                    uDataset.Description = newData;
            end
            app.Datasets(indices(indices(1))) = uDataset;
            saveJSON(app.Datasets,app.Datasets_file);
            msg = msgbox({'The Dataset was updated successfully!!!'},'Success',"custom",icondata,iconcmap);
            movegui(msg,'north');
        end

        % Callback function: LoadPushTool
        function LoadPushToolClicked(app, event)
            folder = uigetdir("Select the Dataset folder");
            dataset_file = fullfile(folder,'BC_VARETA.mat');
            if(isfile(dataset_file))
                BC_VARETA_info = load(dataset_file);
                BC_VARETA_info.Path = folder;
                if(~isempty(app.Datasets) && contains(BC_VARETA_info.Dataset_name,{app.Datasets.Dataset_name}))
                    opts.Interpreter = 'tex';
                    % Include the desired Default answer
                    opts.Default = 'Yes';
                    % Use the TeX interpreter to format the question
                    quest = {strcat("The dataset: ", BC_VARETA_info.Dataset_name, " is already loaded"), ...
                        'Would you like import this the dataset anyway?'};
                    answer = questdlg(quest,'Apply changes',...
                        'Add','Replace','Cancel',opts);
                    % Handle response
                    switch answer
                        case 'Add'
                            app.Datasets(end + 1) = BC_VARETA_info;
                            tableData.Name = {app.Datasets.Dataset_name}';
                            tableData.Description = {app.Datasets.Description}';
                            tableData.Path_Location = {app.Datasets.Path}';
                            app.UITable.Data = struct2table(tableData);
                            disp("-->> Dataset loaded");
                            [icondata,iconcmap] = imread("guide/images/success.png");
                            msg = msgbox({'The Dataset was loaded successfully!!!'},'Success',"custom",icondata,iconcmap);
                            movegui(msg,'north');
                            saveJSON(app.Datasets,app.Datasets_file);
                        case 'Replace'
                            index = find(contains(BC_VARETA_info.Dataset_name,{app.Datasets.Name}),1);
                            app.Datasets(index) = BC_VARETA_info;
                            disp("-->> Dataset replaced");
                            [icondata,iconcmap] = imread("guide/images/success.png");
                            msg = msgbox({'The Dataset was replaced successfully!!!'},'Success',"custom",icondata,iconcmap);
                            movegui(msg,'north');
                            saveJSON(app.Datasets,app.Datasets_file);
                        case 'Cancel'
                    end
                else
                    app.Datasets = BC_VARETA_info;
                    [icondata,iconcmap] = imread("guide/images/success.png");
                    msg = msgbox({'The Dataset was loaded successfully!!!'},'Success',"custom",icondata,iconcmap);
                    movegui(msg,'north');
                    disp("-->> Dataset added");
                    saveJSON(app.Datasets,app.Datasets_file);
                end
            else
                msgbox({'The selected folder do not contains the dataset file.',...
                    ' Please select a correct folder.'},'Error',"error","modal");
            end
        end

        % Callback function: DeletePushTool
        function DeletePushToolClicked(app, event)
            row = app.UITable.Selection;
            if(isempty(row))
                msgbox({'Please select a Dataset to delete.'},'Error',"error","modal");
            else
                opts.Interpreter = 'tex';
                % Include the desired Default answer
                opts.Default = 'Yes';
                % Use the TeX interpreter to format the question
                quest = {strcat("Would you like to delete the Dataset: ", BC_VARETA_info.Dataset_name)};
                answer = questdlg(quest,'Delete',...
                    'Delete','Cancel',opts);
                % Handle response
                switch answer
                    case 'Delete'
                        app.Datasets(row) = [];
                        [icondata,iconcmap] = imread("guide/images/success.png");
                        msg = msgbox({'The Dataset was deleted successfully!!!'},'Success',"custom",icondata,iconcmap);
                        movegui(msg,'north');
                        tableData.Name = {app.Datasets.Dataset_name}';
                        tableData.Description = {app.Datasets.Description}';
                        tableData.Path_Location = {app.Datasets.Path}';
                        app.UITable.Data = struct2table(tableData);
                        disp("-->> Dataset deleted");
                        saveJSON(app.Datasets,app.Datasets_file);
                    case 'Cancel'
                end
            end
        end

        % Callback function: PushTool2
        function PushTool2Clicked(app, event)
            row = app.UITable.Selection;
            if(isempty(row))
                msg = msgbox({'Please select a Dataset to delete.'},'Error',"error","modal");
                movegui(msg,'north');
            else
                msg = msgbox({'Explore dataset.'},'Info',"help","modal");
                movegui(msg,'north');
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create DatasetsUIFigure and hide until all components are created
            app.DatasetsUIFigure = uifigure('Visible', 'off');
            app.DatasetsUIFigure.Position = [100 100 980 587];
            app.DatasetsUIFigure.Name = 'Datasets';
            app.DatasetsUIFigure.Resize = 'off';
            app.DatasetsUIFigure.CloseRequestFcn = createCallbackFcn(app, @DatasetsUIFigureCloseRequest, true);
            app.DatasetsUIFigure.WindowStyle = 'alwaysontop';

            % Create Toolbar
            app.Toolbar = uitoolbar(app.DatasetsUIFigure);

            % Create LoadPushTool
            app.LoadPushTool = uipushtool(app.Toolbar);
            app.LoadPushTool.Tooltip = {'Load Dataset'};
            app.LoadPushTool.ClickedCallback = createCallbackFcn(app, @LoadPushToolClicked, true);
            app.LoadPushTool.Icon = fullfile(pathToMLAPP, 'images', 'upload.png');

            % Create PushTool2
            app.PushTool2 = uipushtool(app.Toolbar);
            app.PushTool2.Tooltip = {'Explore Dataset'};
            app.PushTool2.ClickedCallback = createCallbackFcn(app, @PushTool2Clicked, true);
            app.PushTool2.Icon = fullfile(pathToMLAPP, 'images', 'explore.png');

            % Create DeletePushTool
            app.DeletePushTool = uipushtool(app.Toolbar);
            app.DeletePushTool.Tooltip = {'Delete Dataset'};
            app.DeletePushTool.ClickedCallback = createCallbackFcn(app, @DeletePushToolClicked, true);
            app.DeletePushTool.Icon = fullfile(pathToMLAPP, 'images', 'delete.png');

            % Create PushTool
            app.PushTool = uipushtool(app.Toolbar);
            app.PushTool.Tooltip = {'Help'};
            app.PushTool.Icon = fullfile(pathToMLAPP, 'images', 'help.png');

            % Create UITable
            app.UITable = uitable(app.DatasetsUIFigure);
            app.UITable.ColumnName = {'Name'; 'Description'; 'Location'};
            app.UITable.RowName = {};
            app.UITable.ColumnEditable = [true true false];
            app.UITable.CellEditCallback = createCallbackFcn(app, @UITableCellEdit, true);
            app.UITable.FontSize = 14;
            app.UITable.Position = [10 49 962 534];

            % Create CloseButton
            app.CloseButton = uibutton(app.DatasetsUIFigure, 'push');
            app.CloseButton.ButtonPushedFcn = createCallbackFcn(app, @CloseButtonPushed, true);
            app.CloseButton.Position = [904 16 68 23];
            app.CloseButton.Text = 'Close';

            % Show the figure after all components are created
            app.DatasetsUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Dataset

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.DatasetsUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.DatasetsUIFigure)
        end
    end
end