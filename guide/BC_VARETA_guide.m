classdef BC_VARETA_guide < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        BCVARETAToolboxv10UIFigure  matlab.ui.Figure
        FileMenu                    matlab.ui.container.Menu
        DownloadtestdataMenu        matlab.ui.container.Menu
        ImportdataMenu              matlab.ui.container.Menu
        ExitMenu                    matlab.ui.container.Menu
        ToolsMenu                   matlab.ui.container.Menu
        CreateDataStructureMenu     matlab.ui.container.Menu
        LeadFieldComputationMenu    matlab.ui.container.Menu
        SingleSubjectMenu_LF        matlab.ui.container.Menu
        BatchProcessingMenu_LF      matlab.ui.container.Menu
        MEEGAnalysisMenu            matlab.ui.container.Menu
        SingleSubjectMenu_A         matlab.ui.container.Menu
        BatchProcessingMenu_A       matlab.ui.container.Menu
        ViewMenu                    matlab.ui.container.Menu
        FigureMenu                  matlab.ui.container.Menu
        SubjectsresultMenu          matlab.ui.container.Menu
        SubjectsconnectivityMenu    matlab.ui.container.Menu
        SubjectsactivityMenu        matlab.ui.container.Menu
        RealEEGMenu                 matlab.ui.container.Menu
        HelpMenu                    matlab.ui.container.Menu
        TextArea                    matlab.ui.control.TextArea
    end

    
    properties (Access = private)
        Property % Description
    end
    
    properties (Access = public)
        single_subject % Description
    end
    
    methods (Access = private)
        
        function setPromptFcn(app,jTextArea,eventData,newPrompt)
            % Prevent overlapping reentry due to prompt replacement
            persistent inProgress
            if isempty(inProgress)
                inProgress = 1;  %#ok unused
            else
                return;
            end
            
            try
                % *** Prompt modification code goes here ***
                cwText = char(jTextArea.getText);
                app.TextArea.Value = cwText;
                % force prompt-change callback to fizzle-out...
                pause(0.02);
            catch
                % Never mind - ignore errors...
            end
            
            % Enable new callbacks now that the prompt has been modified
            inProgress = [];
            
        end  % setPromptFcn
        
    end
    

    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            clc;
            processes = jsondecode(fileread(strcat('app_processes.json')));
            for i = 1: length(processes)
                process = processes(i);
                if(process.active)
                    addpath(process.root_folder);
                end
            end
            try
                jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
                jCmdWin = jDesktop.getClient('Command Window');
                jTextArea = jCmdWin.getComponent(0).getViewport.getView;
                set(jTextArea,'CaretUpdateCallback',@app.setPromptFcn)
            catch
                warndlg('fatal error');
            end
        end

        % Menu selected function: CreateDataStructureMenu
        function CreateDataStructureMenuSelected(app, event)
            folder = uigetdir('tittle','Select the Source Folder');
            if(folder==0)
                return;
            end
            create_data_structure(folder);
            msgbox('Completed operation!!!','Info');
        end

        % Callback function
        function ButtonPushed(app, event)
            jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
            jCmdWin = jDesktop.getClient('Command Window');
            jTextArea = jCmdWin.getComponent(0).getViewport.getView;
            cwText = char(jTextArea.getText);
            
            set(jTextArea,'CaretUpdateCallback',@myUpdateFcn)
            
        end

        % Menu selected function: ExitMenu
        function ExitMenuSelected(app, event)
            delete(app);
        end

        % Menu selected function: DownloadtestdataMenu
        function DownloadtestdataMenuSelected(app, event)
            folder = uigetdir('tittle','Select the Source Folder');
            if(folder==0)
                return;
            end
            
            f = dialog('Position',[300 300 250 80]);
            
            iconsClassName = 'com.mathworks.widgets.BusyAffordance$AffordanceSize';
            iconsSizeEnums = javaMethod('values',iconsClassName);
            SIZE_32x32 = iconsSizeEnums(2);  % (1) = 16x16,  (2) = 32x32
            jObj = com.mathworks.widgets.BusyAffordance(SIZE_32x32, 'Downloading test data...');  % icon, label
            
            jObj.setPaintsWhenStopped(true);  % default = false
            jObj.useWhiteDots(false);         % default = false (true is good for dark backgrounds)
            javacomponent(jObj.getComponent, [50,10,150,80], f);
            jObj.start;
            pause(1);
            
            app_properties= jsondecode(fileread(strcat('app_properties.json')));
            url = app_properties.generals.test_data_url;
            filename = strcat(folder,filesep,'BC_VARETA_test_data.zip');
            matlab.net.http.HTTPOptions.VerifyServerName = false;
            options = weboptions('Timeout',Inf,'RequestMethod','get');
            
            try
                disp('Downloding test data....');
                outfilename = websave(filename,url,options);
            catch
                delete(f);
                errordlg('Download error!!!','Error');
                return;
            end
            try
                disp('Unpacking test data....');
                exampleFiles = unzip(filename,folder);
            catch
                delete(f);
                errordlg('Unpacking error!!!','Error');
                return;
            end
            jObj.stop;
            jObj.setBusyText('All done!');
            disp('All done....');
            pause(2);
            delete(f);
            msgbox('Download complete','Info');
        end

        % Menu selected function: FigureMenu
        function FigureMenuSelected(app, event)
            
            [file,path] = uigetfile('*.fig');
            if isequal(file,0)
                disp('User selected Cancel');
                return;
            end
            openfig(strcat(path,filesep,file));
        end

        % Menu selected function: SubjectsresultMenu
        function SubjectsresultMenuSelected(app, event)
            folder = uigetdir('tittle','Select the Source Folder');
            if(folder==0)
                return;
            end
            files = dir(folder);
            ext='.fig';
            for j=1:size(files,1)
                file_name = files(j).name;
                file_path = strcat(folder,filesep, file_name);
                [~,name,ex]=fileparts(file_name);
                %% ----------Searching de data files ------------------------------------
                if(~isfolder(file_path) & strcmpi(strtrim(ex),ext) )
                    openfig(strcat(file_path));
                end
            end
        end

        % Menu selected function: RealEEGMenu
        function RealEEGMenuSelected(app, event)
            %             [file,path] = uigetfile('*.mat');
            %             if isequal(file,0)
            %                 disp('User selected Cancel');
            %                 return;
            %             end
            %             real_EEG=load(strcat(path,filesep,file));
            %
            %             figure_EEG = figure('Color','k','Name',file,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
            %             subplot(1,3,1); plot(real_EEG.Sjj);
            %             %             ylabel('generators')
            %             %     xlabel('generators')
            %             %     title('simulated partial correlations')
            %             subplot(1,3,2); plot(real_EEG.Thetajj);
            %             %             ylabel('generators')
            %             %     xlabel('generators')
            %             %     title('simulated partial correlations')
            %             subplot(1,3,3); plot(real_EEG.indms);
            
        end

        % Menu selected function: SingleSubjectMenu_A
        function SingleSubjectMenu_ASelected(app, event)
            addpath('functions');
            bcv_properties = jsondecode(fileread(fullfile('bcv_properties','bcv_properties.json')));
            bcv_properties.run_single_subject.value = true;
            saveJSON(bcv_properties,fullfile('bcv_properties','bcv_properties.json'));
            BC_VARETA_bash;
            msgbox('Completed operation!!!','Info');
        end

        % Menu selected function: BatchProcessingMenu_A
        function BatchProcessingMenu_ASelected(app, event)
            bcv_properties = jsondecode(fileread(fullfile('bcv_properties','bcv_properties.json')));
            bcv_properties.run_single_subject.value = false;
            saveJSON(bcv_properties,fullfile('bcv_properties','bcv_properties.json'));
            BC_VARETA_bash;
            msgbox('Completed operation!!!','Info');
        end

        % Menu selected function: SingleSubjectMenu_LF
        function SingleSubjectMenu_LFSelected(app, event)
%             addpath('bst_lf_ppl');
%             bs_lf_ppl;
%             msgbox('Completed operation!!!','Info');
        end

        % Menu selected function: BatchProcessingMenu_LF
        function BatchProcessingMenu_LFSelected(app, event)
%             addpath('bst_lf_ppl');
%             bs_lf_ppl;
%             msgbox('Completed operation!!!','Info');
        end

        % Menu selected function: SubjectsconnectivityMenu
        function SubjectsconnectivityMenuSelected(app, event)
            folder = uigetdir('tittle','Select the Subject''s result folder');
            if(folder==0)
                return;
            end
            files = dir(folder);
            ext='.fig';
            for j=1:size(files,1)
                file_name = files(j).name;
                file_path = strcat(folder,filesep, file_name);
                [~,name,ex]=fileparts(file_name);
                %% ----------Searching de data files ------------------------------------
                if(~isfolder(file_path) & strcmpi(strtrim(ex),ext) & contains(file_name,'roi_conn'))
                    openfig(strcat(file_path));
                end
            end
        end

        % Menu selected function: SubjectsactivityMenu
        function SubjectsactivityMenuSelected(app, event)
            folder = uigetdir('tittle','Select the Subject''s result folder');
            if(folder==0)
                return;
            end
            files = dir(folder);
            ext='.fig';
            for j=1:size(files,1)
                file_name = files(j).name;
                file_path = strcat(folder,filesep, file_name);
                [~,name,ex]=fileparts(file_name);
                %% ----------Searching de data files ------------------------------------
                if(~isfolder(file_path) & strcmpi(strtrim(ex),ext) & contains(file_name,'activity') )
                    openfig(strcat(file_path));
                end
            end
        end

        % Menu selected function: ImportdataMenu
        function ImportdataMenuSelected(app, event)
            import_data_structure;
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create BCVARETAToolboxv10UIFigure
            app.BCVARETAToolboxv10UIFigure = uifigure;
            app.BCVARETAToolboxv10UIFigure.Color = [0.9412 0.9412 0.9412];
            app.BCVARETAToolboxv10UIFigure.Colormap = [0.2431 0.149 0.6588;0.251 0.1647 0.7059;0.2588 0.1804 0.7529;0.2627 0.1961 0.7961;0.2706 0.2157 0.8353;0.2745 0.2353 0.8706;0.2784 0.2549 0.898;0.2784 0.2784 0.9216;0.2824 0.302 0.9412;0.2824 0.3216 0.9569;0.2784 0.3451 0.9725;0.2745 0.3686 0.9843;0.2706 0.3882 0.9922;0.2588 0.4118 0.9961;0.2431 0.4353 1;0.2196 0.4588 0.9961;0.1961 0.4863 0.9882;0.1843 0.5059 0.9804;0.1804 0.5294 0.9686;0.1765 0.549 0.9529;0.1686 0.5686 0.9373;0.1529 0.5922 0.9216;0.1451 0.6078 0.9098;0.1373 0.6275 0.898;0.1255 0.6471 0.8902;0.1098 0.6627 0.8745;0.0941 0.6784 0.8588;0.0706 0.6941 0.8392;0.0314 0.7098 0.8157;0.0039 0.7216 0.7922;0.0078 0.7294 0.7647;0.0431 0.7412 0.7412;0.098 0.749 0.7137;0.1412 0.7569 0.6824;0.1725 0.7686 0.6549;0.1922 0.7765 0.6235;0.2157 0.7843 0.5922;0.2471 0.7922 0.5569;0.2902 0.7961 0.5176;0.3412 0.8 0.4784;0.3922 0.8039 0.4353;0.4471 0.8039 0.3922;0.5059 0.8 0.349;0.5608 0.7961 0.3059;0.6157 0.7882 0.2627;0.6706 0.7804 0.2235;0.7255 0.7686 0.1922;0.7725 0.7608 0.1647;0.8196 0.749 0.1529;0.8627 0.7412 0.1608;0.902 0.7333 0.1765;0.9412 0.7294 0.2118;0.9725 0.7294 0.2392;0.9961 0.7451 0.2353;0.9961 0.7647 0.2196;0.9961 0.7882 0.2039;0.9882 0.8118 0.1882;0.9804 0.8392 0.1765;0.9686 0.8627 0.1647;0.9608 0.8902 0.1529;0.9608 0.9137 0.1412;0.9647 0.9373 0.1255;0.9686 0.9608 0.1059;0.9765 0.9843 0.0824];
            app.BCVARETAToolboxv10UIFigure.Position = [100 100 679 459];
            app.BCVARETAToolboxv10UIFigure.Name = 'BC-VARETA Toolbox v1.0';

            % Create FileMenu
            app.FileMenu = uimenu(app.BCVARETAToolboxv10UIFigure);
            app.FileMenu.Text = 'File';

            % Create DownloadtestdataMenu
            app.DownloadtestdataMenu = uimenu(app.FileMenu);
            app.DownloadtestdataMenu.MenuSelectedFcn = createCallbackFcn(app, @DownloadtestdataMenuSelected, true);
            app.DownloadtestdataMenu.Text = 'Download test data';

            % Create ImportdataMenu
            app.ImportdataMenu = uimenu(app.FileMenu);
            app.ImportdataMenu.MenuSelectedFcn = createCallbackFcn(app, @ImportdataMenuSelected, true);
            app.ImportdataMenu.Text = 'Import data';

            % Create ExitMenu
            app.ExitMenu = uimenu(app.FileMenu);
            app.ExitMenu.MenuSelectedFcn = createCallbackFcn(app, @ExitMenuSelected, true);
            app.ExitMenu.Text = 'Exit';

            % Create ToolsMenu
            app.ToolsMenu = uimenu(app.BCVARETAToolboxv10UIFigure);
            app.ToolsMenu.Text = 'Tools';

            % Create CreateDataStructureMenu
            app.CreateDataStructureMenu = uimenu(app.ToolsMenu);
            app.CreateDataStructureMenu.MenuSelectedFcn = createCallbackFcn(app, @CreateDataStructureMenuSelected, true);
            app.CreateDataStructureMenu.Text = 'Create Data Structure';

            % Create LeadFieldComputationMenu
            app.LeadFieldComputationMenu = uimenu(app.ToolsMenu);
            app.LeadFieldComputationMenu.Text = 'Lead Field Computation';

            % Create SingleSubjectMenu_LF
            app.SingleSubjectMenu_LF = uimenu(app.LeadFieldComputationMenu);
            app.SingleSubjectMenu_LF.MenuSelectedFcn = createCallbackFcn(app, @SingleSubjectMenu_LFSelected, true);
            app.SingleSubjectMenu_LF.Text = 'Single Subject';

            % Create BatchProcessingMenu_LF
            app.BatchProcessingMenu_LF = uimenu(app.LeadFieldComputationMenu);
            app.BatchProcessingMenu_LF.MenuSelectedFcn = createCallbackFcn(app, @BatchProcessingMenu_LFSelected, true);
            app.BatchProcessingMenu_LF.Text = 'Batch Processing';

            % Create MEEGAnalysisMenu
            app.MEEGAnalysisMenu = uimenu(app.ToolsMenu);
            app.MEEGAnalysisMenu.Text = 'MEEG Analysis';

            % Create SingleSubjectMenu_A
            app.SingleSubjectMenu_A = uimenu(app.MEEGAnalysisMenu);
            app.SingleSubjectMenu_A.MenuSelectedFcn = createCallbackFcn(app, @SingleSubjectMenu_ASelected, true);
            app.SingleSubjectMenu_A.Text = 'Single Subject';

            % Create BatchProcessingMenu_A
            app.BatchProcessingMenu_A = uimenu(app.MEEGAnalysisMenu);
            app.BatchProcessingMenu_A.MenuSelectedFcn = createCallbackFcn(app, @BatchProcessingMenu_ASelected, true);
            app.BatchProcessingMenu_A.Text = 'Batch Processing';

            % Create ViewMenu
            app.ViewMenu = uimenu(app.BCVARETAToolboxv10UIFigure);
            app.ViewMenu.Text = 'View';

            % Create FigureMenu
            app.FigureMenu = uimenu(app.ViewMenu);
            app.FigureMenu.MenuSelectedFcn = createCallbackFcn(app, @FigureMenuSelected, true);
            app.FigureMenu.Text = 'Figure';

            % Create SubjectsresultMenu
            app.SubjectsresultMenu = uimenu(app.ViewMenu);
            app.SubjectsresultMenu.MenuSelectedFcn = createCallbackFcn(app, @SubjectsresultMenuSelected, true);
            app.SubjectsresultMenu.Text = 'Subject''s result';

            % Create SubjectsconnectivityMenu
            app.SubjectsconnectivityMenu = uimenu(app.ViewMenu);
            app.SubjectsconnectivityMenu.MenuSelectedFcn = createCallbackFcn(app, @SubjectsconnectivityMenuSelected, true);
            app.SubjectsconnectivityMenu.Text = 'Subject''s connectivity';

            % Create SubjectsactivityMenu
            app.SubjectsactivityMenu = uimenu(app.ViewMenu);
            app.SubjectsactivityMenu.MenuSelectedFcn = createCallbackFcn(app, @SubjectsactivityMenuSelected, true);
            app.SubjectsactivityMenu.Text = 'Subject''s activity';

            % Create RealEEGMenu
            app.RealEEGMenu = uimenu(app.ViewMenu);
            app.RealEEGMenu.MenuSelectedFcn = createCallbackFcn(app, @RealEEGMenuSelected, true);
            app.RealEEGMenu.Text = 'Real EEG ';

            % Create HelpMenu
            app.HelpMenu = uimenu(app.BCVARETAToolboxv10UIFigure);
            app.HelpMenu.Text = 'Help';

            % Create TextArea
            app.TextArea = uitextarea(app.BCVARETAToolboxv10UIFigure);
            app.TextArea.Position = [29 21 625 391];
        end
    end

    methods (Access = public)

        % Construct app
        function app = BC_VARETA_guide

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.BCVARETAToolboxv10UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.BCVARETAToolboxv10UIFigure)
        end
    end
end                                                                                                                                                                                                                                              ��h��p�� �'H�?��'�I�M|M����n`*�m ����4�q���0�=M8��h�K9�L8�i��+ydY�(�:B��=���VIMF�c��,`wC� o4��3��� 6�G(�0�T���؃1�m�3lys�XBݑ9��ܐ `c<̌k��x��x�ɥxl������+c���Dą)l6�&�)���G��2�c�9@1��9ca]S�7��A
l��h�
Zcn�1(�u �g����rR�D#0w��X�t!w�E�XzD�[�a,�2��<��¡C��1�΃h`�@�E�D�IÊƁ�"�?�q��#s�s��)y<bba�0�X�L_�2N�J[��XY����m'"R�-U�fio��1d�ʃc�8�_�����"�
�"�2�E��)T��D0g
�RI0ӱH9��`x�L��c���R��T�<:H-e�b����\��H̅��D��3�"v�"ZL�ݺX���X�4	����h8ݍ������<�$�ȍe�P�5U�P�Uv8��5��o�.�"a�m^6�J�q��"�x�$�b� O�Q,V����>2����i�U���/��T�Ȗ����V�T_�M�@ն)��`��b��
>JѼY���;E3;���
a5���Z�h+հ�U��Z�N~<��i��Ȅ�-��tm�Ut�N��$͘7�B`tl4nuOge ���/�h5X�I���#[���B&9A
 ��o|ӫ�66,��&V�M�v����έ�n�,&f�Y�,��X�j9[�Pk5���^��W�Z����L��V���1�'�:���p���';D�l�����ʰU�;*�g��ٷ��y���j�F*_��Fj�s|�c�5��y㇏ۊD&6l(�k{����4�Ac���}4�o9l5�~>Xnd7��SDNZo�ou�$m�؞�u~Y�V��w��u�ʭv�|W�T	<�Hm뇝���n0X��l�K:U��ט��:S5�&u���0��e����I�G�5Fg4�DlU��Ai�m�V=�o��e?^��[�������_�'�d
��R/�n0T�/��SҷȨ?$��Ux�T[hiX>�?zP\=Q�.Ȱ�X�?cs�oei��\�kM�/xm���.�'
b����o�d��^R�-�b	 "�:L\ew&�sn��w�D'g������q�
]oA;^�)D�	1kkx�O�d��J>TN��j���G�c�GPg�����U�����D+�3�6ڷ�~��O��//%�3���?�7�?
���PK    ڬ�N�_q<  R     metadata/appScreenshot.png�XkTSWNi�U�	�h!�TH��Q�F%T�AD���T���[T@�g��"HD�B$b$F!B&7���Z3�Ƶ�c�����|��}��Ϲ7u7�u��������6}�@���o�u#8S�w����M�}b�����f#�@�Sk꿸�7��#�w�۶u>X��a����IR|��A�y��}I;L�6�%U<�'�՗��\D2��o���H����smò��'����Ϝ��u>��'�6�_�tqq�3W� faNNN*N��nnn�/D����lu�(W(B���c���W5�1����sMG��Il܁WU��룈.�&.�`�w��[F�RR�@�䯵	��<�b��t����H�`6��h[�<ĚK�꾖��Y�JZP�ƭ�
{ ����<j��&*�I83����<�&�J��=&̔�eBܩ��"�� �u�����U�MTMk��'����ٿ��2���6�9x_�這�<H?S���zB���Ũt~�H�-�A�fM��0!��g��b����r��n9�2�����y�UnQ��R�l��J�������s���B\�\��5�	X絽�"��@Z����8N%c��~���ܜ���5=QCGW�2��̬��]`	��;���5;��G#v���ʡV����#��F���K�c	kW���U>3E��,�'�q��nq�gP�f �����a�V@8�:݆���Q�::�Q!yZ�05i����"��\H�bW�����|�64(1�1��)�Ͷ���G$fL.�*�"��U0��:ہu�[ED~��DO1��R�9 Q������X������߰e!��W �������k{����pi,�Ec���0��R�n#Y�\�[~+9�R�m��_�P��|���
.�t�����'Y��+����xwrY�����]'3����'���ת���q{�g������ڑ��z���?ف[e{�?�sI;�-+oh�n١߆G�Bʢ�UZ�=iT�^�f��������/���q5���9�D\�PQ]�H5,�>7�Y�'=���0����j8�َ����	�"���]Z�nWx*nT%�Yj�J[.uV��d����0HS��~)�(��=m~wq[s�Y����o����3ẵ�ؿ�=6X7(c!-�(�5�Ġ�ha63s���� �=�#)���?����w���L���g��X����b��qjy!��Z�5G�F�p6i�i���H�<yp+\���5�)˨'�VǫNUō�|7����� n�8̅V�ZҦ�0Wp�ޫ����Mte�3�f}��ݏ��ݼ��^�R��~���-yqg��ο �k�cB'u��3�yTy]3Ͳ�@��k<a�>��f�\�E��Ə�����T����3�nطBj���e���?����[
ĩ{��!E�3����k�moޤ}5���k�Е7��tʦ�<ک�ꩆw�݋�>Ҽ97J�e�S\i���k�����K��tmї�1^v��W�_	lHk$<)%�n��uľ�h�E��6HD��4\�R �^̻ɔ�!v�1��������{��Vr��
+�'�@�KX!�̝o%A:y*�u�<'���:�3�[�+,��w7�&����I�H�V6r���J<�����Q��������*_3J�Q)�Ǹ���#3�����j�03��&Q6��M$Sh��;�R?���m��q��K\�p!| ��eb�����n�p��F�_L�%���i�(ǐ/��DA|�Nd�{YTƟ����@yMeCk��N�b��L�df�wW4�g� �]k���Ml��˳k&c�E�2Oo�Đ�7.^����Y�JxOΰ�{娻��l��ͻ��r�7%b8����os�e ;�]��t�b�C���JA�)#��Gz
�ѧ=�x;)'��w0�v^X��-�����$X���gJ.������H�yS�s����B�g�"�f��hI��AES����b�ڗ��^�I/2[�j��رԢ�����b*� :Τ�0Ӑ�#�2�ĐY�F�)���/���u,�&�@��u�Z"iYN�b�EQ��r�[o+�+UG0�h�)*S���q�"�������
�|��`] �|<��k���7���]5�:T8Y�>O��Cu"�%����Ͻ�f��X�1�Q�
t9����)�Z���,z�� �^��ٖA#�%��"�\o�nxHb��z_{�Wcp%E�8��o�%�g�+\A3�ػ��
��֔x��m ���A�RZB�g��M�j���b�1Dk?��-����o��芤ᭉ��W�E��F�!u)���O� Ǿ�����"uH�V)$������\lg��ˤd�;�<��������;)�߭
����Lª��-i8|�<��A_G�������h��ޱ���f`��i�I"��Sظ���+^���4�Ǻ�n5��W�>#`��(�d=q�zNGc���\=�;8�v)�1B�����ʗ��S,�0W�-�8����Bs�d9~�N_`���/���l'�h!�q��SD{���V�Ǩ��ؕv����Ĭ[�����d>)�u��8�d�����_2{��g�s"�Vi�BZ"KxG�Քv9�T!��{��^�%�����\�b+���Wz��wj��O-y��-_��q�����Qa�K?�%Q��POm���s�_��~�8T��QL�W�SF�``2���8�j���RbFX&�~م�
e�#R`����8J{�����K�bV���P	"]R�Ռ��+�؇��/���l���[��������K���#H�vn����H`I� ���������e3��F��h�-��_eف�������緔�Ԅ�����=�nh�ِn9��=��{�[���j�l�^�������,̓u���Q��	�ѧ<��� =b�A�y��� ��Qq�~�����Q�-w�HdB�U=>T�s�ݞ2�m�����qŲCN��Me	�PW�׺���>J�#��̀�J�G�а���r�s�T�to_�$_��}y����"��5A󆬸��WK=v}ߠ��0�zi�?<��ep�)��G$i��֯�5��Z:�kX����#�)i�L��E�BK������^N�����ݽN�����>�Lh+�3���%���'/�q�����"�6]�P�]�Xh�?ߍ`w���d�_�d�ğ��W>PF>���o�_�m�VL�<�Q�P�z'��R�cͥ��9S,\�?-[e��![5����@��VY"��m�*DSۑ��=�6ִ�B81����:��Y�@t�0p9�]�����9"_Ӱa-����\JZ�/nq��P��R�C����i�zpXl��z�}��E&{���l�����!���,*ΐ��{�����^���m'S��!ty�G�D��4��I��`�$������@ ����?`�����9���?�G�8�B}�>�K�a��� &��&�M��`�'��&6����-�g��%��%��$�_�}��	�Y�>��B�?��O��˦v6�������k����W��U������_�r�w|��i��]~��7Y�w
����I��$��Ĕ,���
*4A�~�7��7:�OPK      qVzN~��   g                  metadata/appMetadata.xmlPK      qVzN
�0D  \              $  metadata/coreProperties.xmlPK      qVzNɨw��   ?              �  metadata/mwcoreProperties.xmlPK      qVzNj�Æ�   �   &            �  metadata/mwcorePropertiesExtension.xmlPK      qVzN(#GqH  �              {  [Content_Types].xmlPK      qVzN>ھJV  &              �  _rels/.relsPK      ��$P�2]�O9  �9               s  appdesigner/appModel.matPK      qVzN�z���  XJ              �@  matlab/document.xmlPK      ڬ�N�_q<  R               Q  metadata/appScreenshot.pngPK    	 	 w  �_                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            p] = chol(XtX);
if ~p,
	invXtX = RXtX \ (RXtX' \ eye(nRegressors));
	pinvX  = RXtX \ (RXtX' \ designMat');
	hasBadEVs    = false;
	badContrasts = false(nContrasts, 1);
else
	% design matrix was rank deficient
	% is that because we have missing information for certain trial types?
	badEVs    = all(0 == designMat);
	hasBadEVs = any(badEVs);
	if hasBadEVs,
		warning([mfilename ':MissingTrials'],                   ...
			    '%s: file %s is missing trials for %d EVs. \n', ...
				mfilename, fileName, sum(badEVs));
		badContrasts = logical(cellfun(@(C) any(C(badEVs)), useContrasts));
		invXtX = pinv(XtX);
		pinvX  = invXtX * designMat';
	else
		error([mfilename ':RankDeficientDesign'],                     ...
			  ['%s: the design matrix is rank deficient. ',           ...
			   'Check that you''ve specified your EVs sensibly. \n'], ...
			  mfilename);
	end%if
end%if

% declare memory
[rho, prho, prhoReg] = deal(zeros(nModes, nModes, nContrasts));

% run GLM on each edge
for i = 1:nModes,
    for j = i+1:nModes,
        rho(i,j,:) = glm_fast_for_meg(squeeze(CorrMats.envCorrelation_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useContrasts, 0);
        prho(i,j,:) = glm_fast_for_meg(squeeze(CorrMats.envPartialCorrelation_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useContrasts, 0);
								  
	    % fill in uninformative values with NaN.
		if hasBadEVs,
			rho(i,j,badContrasts)  = NaN;
			prho(i,j,badContrasts) = NaN;
		end%if
        if isfield(CorrMats, 'envPartialCorrelationRegularized_z'),
        prhoReg(i,j,1:nContrasts) = glm_fast_for_meg(squeeze(CorrMats.envPartialCorrelationRegularized_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useContrasts, 0); 
        prhoReg(i,j,badContrasts) = NaN;
        else
            prhoReg(i,j) = 0;
        end%if
    end%for
end%for

% symmetrise and reformat
for iContrast = nContrasts:-1:1,
    FirstLevel(iContrast).cope.correlation                   = rho(:,:,iContrast) + rho(:,:,iContrast)';
    FirstLevel(iContrast).cope.partialCorrelation            = prho(:,:,iContrast) + prho(:,:,iContrast)';
    FirstLevel(iContrast).cope.partialCorrelationRegularized = prhoReg(:,:,iContrast) + prhoReg(:,:,iContrast)';
end%for

end%run_first_level_glm






%--------------------------------------------------------------------------
function [designMat, goodTrials, trialID, nConditions] = set_up_first_level(D, Settings, excludeTrials)
%SET_UP_FIRST_LEVEL creates the design matrix and trial identifiers

nConditions = length(Settings.SubjectLevel.conditionLabel);

% hold the relevant indices for trials in each condition
for iCondition = nConditions:-1:1,
    tI = D.indtrial(Settings.SubjectLevel.conditionLabel{iCondition}, ...
                    'GOOD');
    trialInds{iCondition} = ROInets.setdiff_pos_int(tI, excludeTrials);
    % check we've found something                               
    if isempty(trialInds{iCondition}), 
        warning([mfilename ':EmptyCondition'], ...
                'No good trials found in %s for condition %s. \n', ...
                D.fname, Settings.SubjectLevel.conditionLabel{iCondition});
    end%if
end%for

% extract a list of all good, relevant trials
goodTrials = sort([trialInds{:}]);

% generate a set of IDs linking trials to condition number
trialID = zeros(length(goodTrials), 1);
for iCondition = nConditions:-1:1,
    trialID(ismember(goodTrials,trialInds{iCondition})) = iCondition;
end%for

% check design matrix size
assert(all(cellfun(@length, Settings.SubjectLevel.designSummary) == nConditions), ...
       [mfilename ':DesignMatrixSizeFault'],                                      ...
       ['The design matrix summary must, in each cell, contain a vector of the ', ...
        'same length as the number of conditions. \n']);
% use an OSL function to generate the subject-specific design matrix
designMat = oat_setup_designmatrix(struct('Xsummary', {Settings.SubjectLevel.designSummary}, ...
                                          'trialtypes', trialID));
end%set_up_first_level
    
%--------------------------------------------------------------------------
function [t, tI, tR] = time_range(time, timeRange, iSession)
%TIME_RANGE selects time range for analysis of each session
% TIME is a vector of times
% TIMERANGE is either a cell array of two-component vectors, a single
% two-component vector, or a null vector

if isempty(timeRange),
    % use the whole time range
    t = time;
    tI = true(size(time));
    tR = [];
else
    % subselect time range
    if iscell(timeRange),
        tR = timeRange{iSession};
    else
        tR = timeRange;
    end%if
    validateattributes(tR, {'numeric'}, ...
                       {'vector', 'numel', 2, 'nondecreasing'}, ... % can have negative times in task data. 
                       'time_range', 'timeRange', 2);
                   
    tI = (time <= tR(2)) & (time >= tR(1));
    t  = time(tI);
end%if
end%time_range



%--------------------------------------------------------------------------
% [EOF]
                                                                                                                                                                                                                                                                                                                                                                                                        (std(voxelData, [], 2), eps);
        voxelDataScaled = ROInets.demean(voxelData, 2);
        clear voxelData
        
        % pre-allocate PCA weightings for each parcel
        voxelWeightings = zeros(size(spatialBasis));
        
        % perform PCA on each parcel and select 1st PC scores to represent
        % parcel
        for iParcel = nParcels:-1:1,
%             progress = nParcels - iParcel + 1;
%             ft_progress(progress / nParcels, ...
%                         [mfilename ...
%                          ':    Finding PCA time course for ROI %d out of %d'], ...
%                         iParcel, nParcels);
                
            thisMask = spatialBasis(:, iParcel);
            if any(thisMask), % non-zero
                parcelData = voxelDataScaled(thisMask, :);
                
                [U, S, V]  = ROInets.fast_svds(parcelData, 1);
                PCAscores  = S * V';
                
                % restore sign and scaling of parcel time-series
                % U indicates the weight with which each voxel in the
                % parcel contributes to the 1st PC
                TSsign          = sign(mean(U));
                relVoxelWeights = abs(U) ./ sum(abs(U)); % normalise the linear combination
                % weight the temporal STDs from the ROI by the proportion used in 1st PC
                TSscale         = dot(relVoxelWeights, temporalSTD(thisMask)); 
                nodeTS          = TSsign .*                               ...
                                  (TSscale / max(std(PCAscores), eps)) .* ... 
                                  PCAscores;
                      
                % return the linear operator which is applied to the data
                % to retrieve the nodeTS
                voxelWeightings(thisMask, iParcel) = TSsign .* ...
                                                     (TSscale / max(std(PCAscores), eps)) ...
                                                     .* U';
                
            else
                warning([mfilename ':EmptySpatialComponentMask'],          ...
                        ['%s: When calculating ROI time-courses, ',        ...
                         'an empty spatial component mask was found for ', ...
                         'component %d. \n',                               ...
                         'The ROI will have a flat zero time-course. \n',  ...
                         'Check this does not cause further problems ',    ...
                         'with the analysis. \n'],                         ...
                         mfilename, iParcel);
                     
                nodeTS = zeros(1, ROInets.cols(voxelDataScaled));
            end%if
            
            nodeData(iParcel,:) = nodeTS;
        end%for
        
        clear parcelData voxelDataScaled

    case 'peakvoxel'
        if any(spatialBasis(:)~=0 & spatialBasis(:)~=1),
            warning([mfilename ':NonBinaryParcelMask'],    ...
                    ['Input parcellation is not binary. ', ...
                     'It will be binarised. \n']);
        end%if
        spatialBasis = logical(spatialBasis);
        
        % find rms power in each voxel
        voxelPower = sqrt(ROInets.row_sum(voxelData.^2) ./ ...
                          ROInets.cols(voxelData));
                      
        % pre-allocate weightings for each parcel
        voxelWeightings = zeros(size(spatialBasis));
                      
        % take peak voxel in each parcel
        for iParcel = nParcels:-1:1,
%             progress = nParcels - iParcel + 1;
%             ft_progress(progress / nParcels, ...
%                         [mfilename ...
%                          ':    Finding peak voxel time course for ROI %d out of %d'], ...
%                         iParcel, nParcels);
            
            thisMask = spatialBasis(:, iParcel);
            
            if any(thisMask), % non-zero
                % find index of voxel with max power
                thisParcPower            = voxelPower;
                thisParcPower(~thisMask) = 0;
                [~, maxPowerInd]         = max(thisParcPower);
                
                % select voxel timecourse
                nodeData(iParcel,:) = voxelData(maxPowerInd,:);
                
                % save which voxel was used
                voxelWeightings(maxPowerInd, iParcel) = 1;
                
            else
                warning([mfilename ':EmptySpatialComponentMask'],          ...
                        ['%s: When calculating ROI time-courses, ',        ...
                         'an empty spatial component mask was found for ', ...
                         'component %d. \n',                               ...
                         'The ROI will have a flat zero time-course. \n',  ...
                         'Check this does not cause further problems ',    ...
                         'with the analysis. \n'],                         ...
                         mfilename, iParcel);
                     
                nodeData(iParcel,:) = zeros(1, ROInets.cols(voxelData));
            end%if
        end%loop over parcels
        
        clear voxelData parcelData
        
    case 'spatialbasis'
        % scale group maps so all have a positive peak of height 1
        % in case there is a very noisy outlier, choose the sign from the
        % top 5% of magnitudes
        top5pcInd = abs(spatialBasis) >=                        ...
                         repmat(prctile(abs(spatialBasis), 95), ...
                                [ROInets.rows(spatialBasis), 1]);
        for iParcel = nParcels:-1:1,
            mapSign(iParcel) = sign(mean(...
                              spatialBasis(top5pcInd(:,iParcel), iParcel)));
        end%for
        scaledSpatialMaps = ROInets.scale_cols(spatialBasis, ...
                                   mapSign ./                ...
                                   max(max(abs(spatialBasis), [], 1), eps));
        
        % find time-course for each spatial basis map
        for iParcel = nParcels:-1:1, % allocate memory on the fly
%             progress = nParcels - iParcel + 1;
%             ft_progress(progress / nParcels, ...
%                         [' ' mfilename ...
%                          ':    Finding spatial basis time course for ROI %d out of %d'], ...
%                         iParcel, nParcels);
            
            % extract the spatial map of interest
            thisMap     = scaledSpatialMaps(:, iParcel);
            parcelMask  = logical(thisMap);
            
            % estimate temporal-STD for normalisation
            temporalSTD = max(std(voxelData, [], 2), eps);
            
            % variance-normalise all voxels to remove influence of
            % outliers. - remove this step 20 May 2014 for MEG as data are
            % smooth and little risk of high-power outliers. Also, power is
            % a good indicator of sensible signal. 
            % Weight all voxels by the spatial map in question
            weightedTS  = ROInets.scale_rows(voxelData, thisMap);
            
            % perform svd and take scores of 1st PC as the node time-series
            % U is nVoxels by nComponents - the basis transformation
            % S*V holds nComponents by time sets of PCA scores - the 
            % timeseries data in the new basis
            [U, S, V]   = ROInets.fast_svds(weightedTS(parcelMask,:), 1);
            clear weightedTS
            
            PCAscores   = S * V';
            maskThresh  = 0.5; % 0.5 is a decent arbitrary threshold chosen by Steve Smith and MJ after playing with various maps.
            thisMask    = thisMap(parcelMask) > maskThresh;   
            
            if any(thisMask), % the mask is non-zero
                % U is the basis by which voxels in the mask are weighted
                % to form the scores of the 1st PC
                relativeWeighting = abs(U(thisMask)) ./ ...
                                    sum(abs(U(thisMask)));
                
                TSsign  = sign(mean(U(thisMask)));
                TSscale = dot(relativeWeighting, temporalSTD(thisMask));       
                nodeTS  = TSsign .*                               ...
                          (TSscale / max(std(PCAscores), eps)) .* ...      
                          PCAscores;
                      
                % for Mark: this is the linear operator which is applied to
                % the voxel data to get nodeTS.
                voxelWeightings(parcelMask,iParcel) = TSsign .* ...
                                             (TSscale / max(std(PCAscores), eps)) ...
                                             .* (U' .* thisMap(parcelMask)');
                
            else
                warning([mfilename ':EmptySpatialComponentMask'],          ...
                        ['%s: When calculating ROI time-courses, ',        ...
                         'an empty spatial component mask was found for ', ...
                         'component %d. \n',                               ...
                         'The ROI will have a flat zero time-course. \n',  ...
                         'Check this does not cause further problems ',    ...
                         'with the analysis. \n'],                         ...
                         mfilename, iParcel);
                     
                nodeTS = zeros(1, ROInets.cols(weightedTS));
                voxelWeightings(~thisMask, iParcel) = zeros(length(thisMask), 1);
            end%if
            
            nodeData(iParcel, :) = nodeTS;
            
        end%loop over parcels
        
        clear voxelData 
        
        
    otherwise
        error([mfilename ':UnrecognisedTimeCourseMethod'],            ...
              ['Unrecognised method for finding ROI time-course. \n', ...
               'Expected ''PCA'', ''spatialBasis'', or ''mean''. \n']);
end%switch

% ft_progress('close');

end%get_node_tcs
%--------------------------------------------------------------------------
end%run_individual_correlation_analysis
%--------------------------------------------------------------------------






%--------------------------------------------------------------------------
function balanced = balance_correlations(x)
%BALANCE_CORRELATIONS makes matrices symmetric

if ismatrix(x) && ~isvector(x),
    balanced = (x + x') ./ 2.0;
else
    balanced = x;
end%if
end%balance_correlations






%--------------------------------------------------------------------------
function [t, tI, tR] = time_range(time, timeRange, iSession)
%TIME_RANGE selects time range for analysis of each session
% TIME is a vector of times
% TIMERANGE is either a cell array of two-component vectors, a single
% two-component vector, or a null vector

if isempty(timeRange),
    % use the whole time range
    t = time;
    tI = true(size(time));
    tR = [];
else
    % subselect time range
    if iscell(timeRange),
        tR = timeRange{iSession};
    else
        tR = timeRange;
    end%if
    validateattributes(tR, {'numeric'}, ...
                       {'vector', 'numel', 2, 'nonnegative', 'nondecreasing'}, ...
                       'time_range', 'timeRange', 2);
                   
    tI = (time <= tR(2)) & (time >= tR(1));
    t  = time(tI);
end%if
end%time_range





%--------------------------------------------------------------------------
function [] = save_corrected_timecourse_results(nodeData, ...
    allROImask, voxelWeightings, Settings, sessionName, protocol, bandName)
% Save the weightings over voxels used to calculate the time-course
% for each ROI
if Settings.SaveCorrected.ROIweightings,
    allVoxelWeightings                = zeros(ROInets.rows(allROImask), ...
                                              ROInets.cols(voxelWeightings));
    allVoxelWeightings(allROImask, :) = voxelWeightings;

    saveDir             = fullfile(Settings.outputDirectory, ...
                                   'spatialBasis-ROI-weightings', filesep);
    ROItcWeightSaveFile = fullfile(saveDir,                      ...
                                   [sessionName '_' protocol '_' ...
                                    bandName '_ROI_timecourse_weightings']);
    ROInets.make_directory(saveDir);
    try
        nii.quicksave(allVoxelWeightings, ROItcWeightSaveFile, Settings.gridStep);
    catch % perhaps we have a weird number of voxels
        save(ROItcWeightSaveFile, 'allVoxelWeightings');
    end%try
end%if

% save node data
if Settings.SaveCorrected.timeCourses,
    saveDir = fullfile(Settings.outputDirectory, 'corrected-ROI-timecourses', filesep);
    ROInets.make_directory(saveDir);
    saveFile = fullfile(saveDir, [sessionName '_correction-' protocol '_' bandName '_ROI_timecourses.mat']);
    save(saveFile, 'nodeData');
end%if

% save variance in each ROI
if Settings.SaveCorrected.variances,
    saveDir = fullfile(Settings.outputDirectory, 'corrected-ROI-timecourses', filesep);
    ROInets.make_directory(saveDir);
    varSaveFile = fullfile(saveDir, [sessionName '_correction-' protocol '_' bandName '_ROI_variances.mat']);
    ROIvariances = var(nodeData, [], 2);                                   %#ok<NASGU>
    save(varSaveFile, 'ROIvariances');
end%if
end%save_corrected_timecourse_results
%--------------------------------------------------------------------------
% [EOF]
                                                                    martshare_history: historycount = %d, peercount = %d
 smartmem: host->memavail = %llu
 smartmem: MemTotal       = %llu (%f GB)
 smartmem: MemFree        = %llu (%f GB)
 smartmem: Buffers        = %llu (%f GB)
 smartmem: Cached         = %llu (%f GB)
 smartmem: NumPeers       = %u
 smartmem: MemReserved    = %llu (%f GB)
 smartmem: MemSuggested   = %llu (%f GB)
 smartcpu_update: switching to zombie
 smartcpu_update: ProcessorCount = %d
 smartcpu_update: NumPeers       = %d
 smartcpu_update: BogoMips       = %.2f
 smartcpu_update: AvgLoad        = %.2f
 smartcpu_update: CpuLoad        = %.2f %%
 smartcpu_update: host->status   = %u
 smartcpu_update: switching to idle
 open_uds_connection socket error: open_uds_connection socket
 open_uds_connection connect error: open_uds_connection connect
 open_uds_connection: connected to %s on socket %d
 open_tcp_connection: using direct memory copy
 open_tcp_connection: server = %s, port = %u
 open_tcp_connection: nslookup1 failed on '%s'
 open_tcp_connection: nslookup2 failed on '%s'
 open_tcp_connection: socket = %d
 open_tcp_connection error: open_tcp_connection
 open_tcp_connection: connectioncount = %d
 open_tcp_connection: connected to %s:%u on socket %d
 close_connection: socket = %d
 close_connection error: close_connection
 close_connection: connectioncount = %d
                          �  4   4   _�      4                                   zR x�  $      Z������=        A�C       $   D   0Z�������        A�C       $   l   �Z�������       A�C       $   �   �\�������        A�C       $   �   �\������2	       A�C              zR x�  $      �e������i        A�C       $   D   @f������K       A�C       $   l   hg�������       A�C       $   �   0j�������       A�C              zR x�  $      �m������       A�C       $   D   �o������(
       A�C              zR x�  $      xy�������        A�C       $   D   0z������R       A�C       $   l   h��������       A�C              zR x�  $      ؂������V        A�C       $   D   �������U       A�C       $   l   H�������2       A�C              zR x�  $      H�������m        A�C       $   D   ��������:       A�C       $   l   ��������g       A�C       $   �   ��������"       A�C       $   �   ��������p        A�C       $   �   @�������p        A�C       $     ��������_        A�C       $   4  ���������        A�C       $   \  h�������F       A�C       $   �  ���������        A�C       $   �  ��������        A�C       $   �  ���������        A�C       $   �  8��������        A�C       $   $  Д�������        A�C       $   L  h��������        A�C       $   t   ��������        A�C       $   �  ��������        A�C       $   �  ���������        A�C              zR x�  $      ��������>       A�C       $   D   ؗ������v       A�C              zR x�  $      �������"       A�C       $   D    ��������       A�C       $   l   ��������#        A�C              zR x�  $      ��������       A�C       $   D   ���������       A�C              zR x�  $      0��������        A�C       $   D   ض������P       A�C       $   l    �������P       A�C       $   �   (�������P       A�C              zR x�  $      8��������       A�C              zR x�  $      ��������b        A�C       $   D   м������J       A�C       $   l   ��������R       A�C              zR x�  $      ��������        A�C       $   D   ���������
       A�C              zR x�  $      0��������        A�C       $   D   ��������       A�C              zR x�  $      ���������       A�C       $   D   P�������)       A�C       $   l   X�������#        A�C       $   �   `�������       A�C        �_�  �_�          ��      �      �      �      (�      4�      @�      L�      X�      d�      p�      |�      ��      ��      ��      ��      ��      Ĝ      М      ܜ      �      ��       �      �      �      $�      0�      <�      H�      T�      `�      l�      x�      ��      ��      ��      ��      ��      ��      ̝      ؝      �      �      ��      �      �       �      ,�      8�      D�      P�      \�      h�      t�      ��      ��      ��      ��      ��      ��      Ȟ      Ԟ      ��      �      ��      �      �      �              ��      ��              ���2                                                            ���2                                                                   ���<                                            ���2                                                            ���2                                                            ���2                                                            ���2                                                            ���2                                                            ���2                                                            ���2                                                            ���2                                                            ���2                                                            ���2                                                            ���2                                                            ���2                                                            ���2                                                            ���2                                                            ���2                                                            ���2                                                            ���2                                                            ���2                                                            ���2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       (     0     8     @     H     P     X     `     h     p     x     �     �     �     �     �     �     �     �     �     �     �     �     �     �     �     �                           (    0    8    @    H    P    X    `    h    p    x    �    �    �    �    �    �    �    �    �    �    �    �    �    �    �    �                          (    0    H    @    �@����p��������`��p���pp`��������������0�������p��������0        d           +   d           8   f ��U       .  �      x   $  �         $   @          N  @          .         �   $            $   �          N  �          .  �      �   $  �         $   �         N  �         .  �      �   $  �         $   �          N  �          .  @      �   $  @         $   2	         N  2	      �               �               �               �               �                                d             d             d           !  f ��U       .  �      T  $  �         $   p          N  p          .  �      [  $  �         $   P         N  P         .  @      m  $  @         $   �         N  �         .  0!      w  $  0!         $   �         N  �         d             d           �  d           �  f ��U       .  �$      �  $  �$         $             N            .  �&      �  $  �&         $   (
         N  (
         d             d           �  d           �  f ��U       .   1        $   1         $   �          N  �          .   2      *  $   2         $   `         N  `         .  `8      2  $  `8         $   �         N  �         d             d           B  d           K  f ��U    |              �              �              �              �              �              �              �              �              	                            +              @              U              k                            �              �              �              �              �              �              �                                          "              0              =              J              W              h              n              x              �              �              �              �              �              �              �              �              �                                                           d             d           #  d           .  f ��U       .  ;      a  $  ;         $   `          N  `          .  p;      g  $  p;         $   `         N  `         .  �A      q  $  �A         $   2         N  2         d             d           {  d           �  f ��U       .  D      �  $  D         $   p          N  p          .  �D      �  $  �D         $   @         N  @         .  �E      �  $  �E         $   p         N  p         .  0G      �  $  0G         $   0         N  0         .  `I      �  $  `I         $   p          N  p          .  �I      �  $  �I         $   p          N  p          .  @J      �  $  @J         $   `          N  `          .  �J      �  $  �J         $   �          N  �          .  pK      
  $  pK         $   P         N  P         .  �L        $  �L         $   �          N  �          .  `M      /  $  `M         $   �          N  �          .   N      D  $   N         $   �          N  �          .  �N      Z  $  �N         $   �          N  �          .  �O      o  $  �O         $   �          N  �          .  `P      �  $  `P         $   �          N  �          .   Q      �  $   Q         $   �          N  �          .  �Q      �  $  �Q         $             N            .  �Q      �  $  �Q         $   �          N  �          d             d           �  d           �  f ��U       .  pR        $  pR         $   @         N  @         .  �S        $  �S         $   v         N  v         d             d           )  d           5  f ��U       .  0Z      i  $  0Z         $   0         N  0         .  `[      |  $  `[         $   �         N  �         .   b      �  $   b      �  �              $   #          N  #          d             d              d           ,  f ��U       .  Pb      `  $  Pb         $             N            .  pc      s  $  pc         $   �         N  �         d             d           ~  d           �  f ��U       .  @s      �  $  @s         $   �          N  �          .  t      �  $  t         $   P         N  P         .  `u      �  $  `u         $   P         N  P         .  �v      �  $  �v         $   P         N  P         d             d           	  d           	  f ��U       .   x      F	  $   x         $   �         N  �         d             d           W	  d           d	  f ��U       .  �z      �	  $  �z         $   p          N  p          .   {      �	  $   {         $   P         N  P         .  P      �	  $  P         $   R         N  R         d             d           �	  d           �	  f ��U       .  ��      
  $  ��         $   �          N  �          .  P�      
  $  P�         $   �
         N  �
         d             d           /
  d           :
  f ��U       .  0�      m
  $  0�         $   �          N  �          .  Ѝ      |
  $  Ѝ         $            N           d             d           �
  d           �
  f ��U       .  �      �
  $  �         $   �         N  �         .  ��      �
  $  ��         $   0         N  0         .  �      �
  $  �         $   0          N  0          .  @�        $  @�         $            N           d              �      ,    �      ?    �      S           f    �      r    �      {    �      �    �      �    @      �    0!      �    �$      �    �&      �     1      �     2      �    `8      �    ;      �    p;          �A          D          �D      !    �E      +    0G      3    `I      =    �I      H    @J      T    �J      d    pK      s    �L      �    `M      �     N      �    �N      �    �O      �    `P      �     Q          �Q          �Q      %    pR      8    �S      C    0Z      V    `[      a     b      o    Pb      �    pc      �    @s      �    t      �    `u      �    �v      �     x      �    �z      �     {          P           ��      /    P�      @    0�      O    Ѝ      `    �      u    ��      �    �      �    @�      �    X�      �    ��      �    ��      �    ��      �    �      �    P�          ��          ��      &    �      <    P�      G    ��      V    ��      d    �      x    P�      �    ��      �    ��      �    �      �    P�      �    ��      �    ��      �    �          P�          ��      -    ��      6    ��      ?    ��      Q    ��      c    ��      t    ��      �    ��      �    ��      �     �      �    �      �    �      �    �      �    �      �    �      �     �      �    (�          0�          8�      #    @�      3    H�      D    P�      S    X�      c    `�      s    h�      }    ��      �    ��      �    ��      �    ��      �    @      �            �            �            �            �            �                                    #            )            0            7            @            G            Q            W            ^            d            q            z            �            �            �            �            �            �            �            �            �            �            �            �            
                        #            4            J            d            n            �            �            �            �            �            �            �            �            �                                    ,            @            F            O            U            _            e            m            y            �            �            �            �            �            �            �            �            �            �            �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �                     	  
            �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �  �                     	  
              /Users/roboos/matlab/fieldtrip/peer/src/ memprofile.c /Users/roboos/matlab/fieldtrip/peer/src/../private/memprofile.o _memprofile_cleanup _memprofile_sample _memprofile _exitFun _mexFunction _mutexmemlist _mutexmemprofile _reftime _memlist _memprofileStatus _memprofileThread announce.c /Users/roboos/matlab/fieldtrip/peer/src/announce.o _frand _cleanup_announce _announce _announce_once discover.c /Users/roboos/matlab/fieldtrip/peer/src/discover.o _cleanup_discover _discover expire.c /Users/roboos/matlab/fieldtrip/peer/src/expire.o _cleanup_expire _expire _check_watchdog extern.c /Users/roboos/matlab/fieldtrip/peer/src/extern.o _syslog_level _condstatus _mutexstatus _mutexappendcount _mutexsocketcount _mutexthreadcount _mutexconnectioncount _mutexhost _mutexpeerlist _mutexjoblist _mutexallowuserlist _mutexrefuseuserlist _mutexallowgrouplist _mutexrefusegrouplist _mutexallowhostlist _mutexrefusehostlist _mutexwatchdog _mutexsmartmem _mutexsmartcpu _mutexprevcpu _mutexsmartshare _udsserverStatus _tcpserverStatus _announceStatus _discoverStatus _expireStatus _appendcount _socketcount _threadcount _connectioncount _host _peerlist _joblist _allowuserlist _refuseuserlist _allowgrouplist _refusegrouplist _allowhostlist _refusehostlist _smartsharelist _watchdog _smartmem _smartcpu _prevcpu _smartshare peerinit.c /Users/roboos/matlab/fieldtrip/peer/src/peerinit.o _hash _peerinit _peerexit util.c /Users/roboos/matlab/fieldtrip/peer/src/util.o _threadsleep _bufread _bufwrite _append _jobcount _peercount _hoststatus _clear_peerlist _clear_joblist _clear_smartsharelist _clear_allowuserlist _clear_allowgrouplist _clear_allowhostlist _clear_refuseuserlist _clear_refusegrouplist _clear_refusehostlist _check_datatypes _getmem udsserver.c /Users/roboos/matlab/fieldtrip/peer/src/udsserver.o _cleanup_udsserver _udsserver tcpserver.c /Users/roboos/matlab/fieldtrip/peer/src/tcpserver.o _cleanup_tcpserver _tcpserver __OSSwapInt16 /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/usr/include/libkern/i386/_OSByteOrder.h tcpsocket.c /Users/roboos/matlab/fieldtrip/peer/src/tcpsocket.o _cleanup_tcpsocket _tcpsocket security.c /Users/roboos/matlab/fieldtrip/peer/src/security.o _security_check _ismember_userlist _ismember_grouplist _ismember_hostlist localhost.c /Users/roboos/matlab/fieldtrip/peer/src/localhost.o _check_localhost smartshare.c /Users/roboos/matlab/fieldtrip/peer/src/smartshare.o _smartshare_reset _smartshare_check _smartshare_history smartmem.c /Users/roboos/matlab/fieldtrip/peer/src/smartmem.o _smartmem_info _smartmem_update smartcpu.c /U