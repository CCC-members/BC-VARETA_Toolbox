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
end                                                                                                                                                                                                                                              ¼ñhˆp‡ Â'H¬?äÌ'ÓI¡M|Mûœôñn`*m ÃÑ4Çq¦ãş0ç=M8äähK9˜L8Íi¿Ÿ+ydYŞ(ô:BÑù=¨‹÷VIMFÎc‚ä,`wC° o43‚¹‚ 6G(¢0ÄT×ò¥Øƒ1ŠmÁ3lysóXBİ‘9æ™Ü `c<ÌŒkõ‡xµx›É¥xl¾œ€± Ñ+c£°çDÄ…)l6Ï&ƒ)¢›‹G„ 2¯c£9@1ÇÂ9ca]SÈ7šA
lÖ‹hß
Zcn1(Âu àgœÏÍÖrRÜD#0w·ÍX˜t!wœEãXzDº[Øa,†2„¦<ÊÆÂ¡CğÌ1çÎƒh`‹@šED”IÃŠÆà"İ?áq‰†#s’sˆü)y<bba£0øX„L_¨2N”J[‚ÉXY‚“¶„m'"R…-Uáfioà1d¸Êƒcñ8ë_¤·úÇ"£
"¦2ıEö…)T£èD0g
RI0Ó±H9–è`xêL…Êc‘„ÊR©äT¨<:H-eãb‡ˆØá\–©HÌ…œ¦D°„3¦"v¦"ZL¶İºXµÃ»XË4	ş§äh8İÚêòâô‰<†$Èe˜PÀ5UêPUv8‡¡5ü¡oà .­"aÈm^6´J²q—“"Ùxê¨$¡bÚ O²Q,VÎôƒéå>2µëü¸iƒUˆí/ÛñTöÈ–†³ÕùVë´T_ŠMå@Õ¶)£í`‹æb¯Õ
>JÑ¼YñÚŞ;E3;¶­²
a5«ÛûZÌh+Õ°¡U§ìZN~<íái—ÑÈ„á-­ÒtmóUtÛN¿•$Í˜7ÖB`tl4nuOge ÕáÜ/ƒh5XÛI¾ÛÖ#[›¬B&9A
 ­¢o|Ó«˜66,‰Ü&VƒMÄvğø¦—Î­înä©,&f¯YÚ,Õñ‚X±j9[ùPk5“ÛÁ^›ßW·Z¬…¯²Lòç­V«­1ª'İ:ª˜™pø€¨';Dşl‹¸™•å«Ê°UìÍ¾*ÜgˆöÙ·öºyšäîj F*_áİFjæºs|«cÅ5úÙyã‡ÛŠD&6l(Ük{¡Âñå’4ŠAc©¾Ü}4®o9l5~>Xnd7ìĞSDNZo«ou—$mİØ£u~YŠVß÷wš¤uñ®Ê­vŠ|WóT	<ÀHmë‡¦ªñn0X¦ÓlK:U®ö×˜Ëş:S5­&ušÉî0‘İe¦‚ŠIñG¶5Fg4ìDlU½ÈAiŒmV=ÕoÔÈe?^ğğ®[‘ÈÄÁ†öé_ñ'Šd
‡êR/n0TéŠ/—úSÒ·È¨?$ƒ©UxµT[hiX>‘?zP\=Q‹.È°¾X‚?cs’oei±ƒ\½kM‹/xmµîò.ñ'
b¾€ŠÌo½dÅŒ^Ræ-—b	 "ü:L\ew&‰snòôwıD'g”¥Óáéùqñ
]oA;^†)DĞ	1kkxOí¸dï¢ÿJ>TNªjŸÓóGñcGPg«š¿äòU³ÛŞö˜D+¯3æ¤6Ú·ª~ÃÿOŸø//%ü3¿ÂÄ?Å7æ?
õü¿PK    Ú¬¤Nˆ_q<  R     metadata/appScreenshot.pngíXkTSWNi‹Uğ	‚h!¾THâûQõF%TÄAD‰‘©TÄÇØ[T@‡g’©"HD¼B$b$F!B&7ÁÖÎZ3³Æµúcòëæìï|ûì}¾½Ï¹7u7°uú´ùÓÄôíÛ6}@Ä‚oŒu#8SËw„…ÕöMë}b¸ı»ˆf#–@µSkê¿¸´7ùİ#Íw­Û¶u>X¾ŞaÙè‘ßÏIR|ùàAÔyñİ}I;LŠ6ç%U<—'‰Õ—ƒû\D2ü¶oíÎÚHœ×¾µsmÃ²Ÿ'´£»ÃÏœ­—u>°À'Ó6«_­tqqÁ3W“ faNNN*N÷Ännn®/DÍƒ‚‚luã(W(BÕşúcÒò°ŸW5î1ìô¿æsMG£¥IlÜWUª‹ë£ˆ.Ô&.`„wáÇ[F‰RRì@Ÿä¯µ	®Ô<…bâàtâ±àÇÇHÈ`6°³h[ª<ÄšK´ê¾–çÈYğJZPáÆ­õ
{ ³“·¶<j¿Í&*ÁI83Š‘ğÎ<®&ŠJ½½=&Ì””eBÜ©ïÆ"üÖ îuìñéò½U»MTMkÊä'´§¶®Ù¿À”2ÜéÈ6±9x_Úé€™Ø<H?SÑÄãzB±««Å¨t~äH‚-†A¢fM¹±0!­ïg¦è²bˆªµr»ºn9š2Óîµ‚°y¿UnQÀí”RlëúJ¬ª€¬ÊÃös¸¤B\Ö\¢ó’5ã	Xçµ½ª"ÿö@Z€¯Àß8N%c£Ø~˜­ÔÜœœÙ5=QCGW§2ÏÈÌ¬ùš]`	²Á;¬›í5;¤ÔG#vè…úÊ¡VÁö÷ïº#åôF³ö˜Kc	kW¬ÀúU>3EŸÔ,Ş'¯q¶ºnqºgPèf »±‘ääaV@8Ê:İ†ö½£QÆ::®Q!yZù05iš½˜Â"æ¾ş\HÀbW¿¼¡ÿ|é64(1¹1Ÿä)«Í¶»¿ÊG$fL.·*Ê"”U0®‡:ÛuĞ[ED~ŒßDO1¬R³9 Q½áçŒ‡ÿXÒ÷á†ÏÁß°e!Çã¨W Ÿ¼›úñÁÄk{ü¶¼ïpi,ÿEcâÆñ»0ù¦R¨n#YŠ\ˆ[~+9šRƒm”ô_åP¯ª|¨›Ü
.ÒtƒªÇà'YœÜ+ÅŞéç—xwrY¾½îÅÈ]'3µ°¬'§“×ª’ö€q{í—gËÉı§ÑÚ‘¬zÇŞà?Ù[e{ú?ã‘sI;¶-+ohènÙ¡ß†G¬BÊ¢÷UZ“=iT³^Šf Ğ¶åèêÀî/€Îôq5˜š9àD\íPQ]¢H5,>7ŒY¨'=Š¡Î0’ñÏæj8ıÙüÒà–	­"ÏúÍ]ZénWx*nT%¦YjúJ[.uV‘–d½ãÕ0HS¡¼~)£(€¢=m~wq[só‰Y ¥ÓÁoÜ¦¥3áºµ²Ø¿Ú=6X7(c!-½(ã5˜Ä  ha63sÁ‚ì ²=°#)³è?çÕ××w•L¶Ã˜g¬§X°›û–b©¶qjy!ªüZÓ5GFÊp6iäi˜‰±HÙ<yp+\ûéÂ5„)Ë¨'æVÇ«NUÅ|7ª°”å n©8Ì…VœZÒ¦Ü0WpôŞ« ªMte‚3øf}Œüİ ‚İ¼“ª^ôR´†~ÕÕá-yqgíÎ¿ ¥kŠcB'uÉÂ—°3äyTy]3Í²@±¦k<a·>õºfš\E«–Æ´ùÂàT›‡¨ş3ínØ·BjİÅ¤eôµò?£ô™µ[
Ä©{İä!EÙ3… éÚkÌmoŞ¤}5 ¹ÉkÜĞ•7¢âtÊ¦ƒ<Ú©¡ê©†wıİ‹¤>Ò¼97J…e¤S\iô–Êkâü¶‚ÁK°´tmÑ—Ú1^v–ëWÏ_	lHk$<)%½n¦™uÄ¾©hè¾Eâ6HD¿˜4\ÓR ‚^Ì»É”¦!v¡1µµš†ÊÊ¸‚{ßØVrÆÆ
+ƒ'›@í¢©KX!ÃÌo%A:y*µu<'º£â¯:º3Ê[¦+,Éw7Š&íÓÌIËH¨V6r«äòJ<¬é‡ìÏúQ¢—® ¼à*_3JÍQ)íÇ¸ÛÃú#3Ÿ¼–¢Šjı03ıí&Q6ªÂM$Sh¡”;×R?„§úm‘ÖqŸ£K\±p!| °ëebŸºƒÃÕnÎpù F‚_L%¿Œşi’(Ç/áËDA|ıNdğ{YTÆŸ¿Üä@yMeCkÛÌNÀb¨«Lºdf­wW4ƒg• ­]k’’âMlş¸Ë³k&c«E¥2Oo¨Äñ¦7.^ñ¡Òü…YÅJxOÎ°À{å¨»ùlÁ¬Í»ßÂr¸7%b8ì›œ¿osÅe ;Ü]öt¨bÁCÿ×µJÂšAâ)#‘öGz
‡Ñ§=ëx;)'ÛÄw0¢v^Xµ´-Á‰¢å÷$Xè‘ïgJ.ƒ¿¥ûîH„ySìs³ÃäìB¹gµ"ôfûœhIÈÛAES‰‹€—bÕÚ—Á¢^ÎI/2[Àj“îØ±Ô¢²èÈäó•b*¿ :Î¤¾0Ó™#¢2«ÄYéF…)ÅİÃ/¬¦Äu,¨&@­è¢u÷Z"iYN³b†EQŠ©r¥[o+Á+UG0™hˆ)*S•Éëq‚"‰¦äñö§
•|´á`] »|<ÿèkúş¥7‚Ñù]5Ø:T8Yû>OÄê”Cu"³%¨”«üÏ½¬f×ØXÏ1ã¢QĞ
t9“ƒò·)ÑZ•»°,zÈÑ ß^³Ù–A#±%ƒ›"º\o—nxHb‡ˆz_{„Wcp%Eì8Æâoä%Üg„+\A3çØ»µ¯
ûÃÖ”x„¼m ’¤¦A¿RZB¬g»ÔMÃj±­ b°1Dk?¨ğ-ßôµ©o·äèŠ¤á­‰äÎWåŒE‡šF¯!u)š¨îO§ Ç¾­îÇöú"uHıV)$‹´µ·îÕ\lgÌúË¤dı;²<ö˜ƒ†¾…òş;)ôß­
Éı“LÂª‡Å-i8|•<øÛA_G©àõ¯‘‚Ûh›ÄŞ±àúæf`àÊiêI"±¤SØ¸—Á·+^ ´Â4çÇºön5ğúW³>#`­±(d=qùzNGcçÃè´\=Ò;8ôv)Ô1B’•‹ÖÆÊ—­S,á0W¢-¾8”Ÿ¢ë§BsÄd9~±N_`ºıå/‡˜¥l'Õh!ÀqˆùSD{ìœú´V·Ç¨¯¤Ø•v —à½äÄ¬[‰ùôùÑd>)ÒuÀ¡8ç¬dÂìõ‚_2{ö‰gçs"½VióBZ"KxG‰Õ”v9ÃT!Úü{ùÂ•Ï^û%°œ€äÃ\‰b+±÷ŒWz„×wjçøO-y£ã­-_Åäqé°ï÷İQaÓK?è©%Q£ÓPOmŠ¢s÷_ô›~–8Tñ³ñQL”WâºSFÃ``2ñŞÉ8¹j…İõRbFX&ó~Ù…ö
eÚ#R`¡–“¾8J{éû£åÙKÊbVØşÔP	"]RÕŒÍ+–Ø‡÷£/Ÿ¯ól“¦Í[±äÁì«åĞÓÄK°ù#Hšvn··­ÁH`I Ëø‹ñÖÙÙÔe3º²FÎh§-Çñ_eÙ·‡øÌàóç·”´Ô„÷”ù£ã=ınhúÙn9¢à=‘™{[ïšİÈjÈl™^»¶¸üÔÁ,Ì“uÄÿĞQñä	ÉÑ§<—äêƒ =b™Aİyê»ÜÙ „’QqÇ~Êú‘®ñƒ„Qß-wÛHdB‰U=>T£sºİ2ªm´ü–­»qÅ²CN…ØMe	¹PWÚ×º³¢>Jò#³¿Ì€¬JÚGËĞ°Øçºršs¸TÕto_²$_åû}yşºö€"İÆ5Aó†¬¸ÊóWK=v}ß »®0Äziéš?<İØepô)ÁG$iîìÖ¯œ5óÔZ:êkXºŸ±#Û)iÊLıİEÜBKÃÁ¼ù—¡^N‰‰¾ÑÍİ½N¡ˆêŞ>¤Lh+õ3¡‡¾%„¤î'/„qŠ“ŠŠŠ"Ş6]‹P½]¼Xhï?ß`w¼‚Ãdí_ÎdùÄŸğW>PF>¬®Æoª_êm”VL£<¦QÆPz'¨ËR’cÍ¥‘Ó9S,\ï?-[eÜ![5‰—­²@†ÉVY"ãám½*DSÛ‘”±=›6Ö´ÈB81ş‡â—:òœĞYÂ@t…0p9º]ø§†Óã9"_Ó°a-û¶²á\JZÉ/nqĞŞPùêRÆCåşÈĞiğzpXlàÒzÜ}‡ïE&{……l…ÎåŞÍ!ñ÷™,*Î•…{Éí½ËÊ^¦ÏÂm'SàÌ!ty¶GÀDØå4ÃƒIÿß`ú$ôœ®êàí@ êíÿä?`ûºø–9ñû”?ÛGĞ8òB}Œ>ÛKÊa¬œÉ &“ƒ&‚M““`Ó'¡ñ™&6ÄóĞÿ-Ûgéÿ%ÿ¹%ÿ¹$ú_³}õ	üYØ>‚ÖB¶? ˆOÉàË¦v6ıŒ’‡ÊÓàk¼‹—ŸW…şUş¯»éÃ_„rôw|ıËiş]~€ú7Y—w
òü¦ÉIÓä$şÔÄ”,Àı
*4Aè~Û7›Š7:ÿOPK      qVzN~€î   g                  metadata/appMetadata.xmlPK      qVzN
ñ0D  \              $  metadata/coreProperties.xmlPK      qVzNÉ¨wÇÃ   ?              ¡  metadata/mwcoreProperties.xmlPK      qVzNjûÃ†˜   Ò   &            Ÿ  metadata/mwcorePropertiesExtension.xmlPK      qVzN(#GqH  å              {  [Content_Types].xmlPK      qVzN>Ú¾JV  &              ô  _rels/.relsPK      ö“$PÑ2]­O9  Š9               s  appdesigner/appModel.matPK      qVzN“z›Œî  XJ              ø@  matlab/document.xmlPK      Ú¬¤Nˆ_q<  R               Q  metadata/appScreenshot.pngPK    	 	 w  ‹_                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            p] = chol(XtX);
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
                             4   4   _š      4                                   zR x  $      Zÿÿÿÿÿÿ=        A†C       $   D   0Zÿÿÿÿÿÿ¸        A†C       $   l   ÈZÿÿÿÿÿÿî       A†C       $   ”   \ÿÿÿÿÿÿ‚        A†C       $   ¼   ø\ÿÿÿÿÿÿ2	       A†C              zR x  $      øeÿÿÿÿÿÿi        A†C       $   D   @fÿÿÿÿÿÿK       A†C       $   l   hgÿÿÿÿÿÿë       A†C       $   ”   0jÿÿÿÿÿÿ–       A†C              zR x  $      mÿÿÿÿÿÿ       A†C       $   D   ˆoÿÿÿÿÿÿ(
       A†C              zR x  $      xyÿÿÿÿÿÿÖ        A†C       $   D   0zÿÿÿÿÿÿR       A†C       $   l   h€ÿÿÿÿÿÿ¬       A†C              zR x  $      Ø‚ÿÿÿÿÿÿV        A†C       $   D   ƒÿÿÿÿÿÿU       A†C       $   l   H‰ÿÿÿÿÿÿ2       A†C              zR x  $      H‹ÿÿÿÿÿÿm        A†C       $   D   ‹ÿÿÿÿÿÿ:       A†C       $   l   ¨Œÿÿÿÿÿÿg       A†C       $   ”   ğÿÿÿÿÿÿ"       A†C       $   ¼   øÿÿÿÿÿÿp        A†C       $   ä   @ÿÿÿÿÿÿp        A†C       $     ˆÿÿÿÿÿÿ_        A†C       $   4  ÀÿÿÿÿÿÿÂ        A†C       $   \  h‘ÿÿÿÿÿÿF       A†C       $   „  ’ÿÿÿÿÿÿ—        A†C       $   ¬  “ÿÿÿÿÿÿ¿        A†C       $   Ô   “ÿÿÿÿÿÿ¿        A†C       $   ü  8”ÿÿÿÿÿÿ¿        A†C       $   $  Ğ”ÿÿÿÿÿÿ¿        A†C       $   L  h•ÿÿÿÿÿÿ¿        A†C       $   t   –ÿÿÿÿÿÿ¿        A†C       $   œ  ˜–ÿÿÿÿÿÿ        A†C       $   Ä  €–ÿÿÿÿÿÿ€        A†C              zR x  $      À–ÿÿÿÿÿÿ>       A†C       $   D   Ø—ÿÿÿÿÿÿv       A†C              zR x  $      ÿÿÿÿÿÿ"       A†C       $   D    Ÿÿÿÿÿÿÿ¿       A†C       $   l   ¸¥ÿÿÿÿÿÿ#        A†C              zR x  $      ¨¥ÿÿÿÿÿÿ       A†C       $   D    ¦ÿÿÿÿÿÿÏ       A†C              zR x  $      0¶ÿÿÿÿÿÿÊ        A†C       $   D   Ø¶ÿÿÿÿÿÿP       A†C       $   l    ¸ÿÿÿÿÿÿP       A†C       $   ”   (¹ÿÿÿÿÿÿP       A†C              zR x  $      8ºÿÿÿÿÿÿ…       A†C              zR x  $      ˆ¼ÿÿÿÿÿÿb        A†C       $   D   Ğ¼ÿÿÿÿÿÿJ       A†C       $   l   øÀÿÿÿÿÿÿR       A†C              zR x  $      Ãÿÿÿÿÿÿ“        A†C       $   D   ÃÿÿÿÿÿÿÆ
       A†C              zR x  $      0Îÿÿÿÿÿÿ        A†C       $   D   ¨Îÿÿÿÿÿÿ       A†C              zR x  $      ˆÔÿÿÿÿÿÿç       A†C       $   D   PÖÿÿÿÿÿÿ)       A†C       $   l   XÙÿÿÿÿÿÿ#        A†C       $   ”   `Ùÿÿÿÿÿÿ       A†C        À_ÿ  À_ÿ          ø›      œ      œ      œ      (œ      4œ      @œ      Lœ      Xœ      dœ      pœ      |œ      ˆœ      ”œ       œ      ¬œ      ¸œ      Äœ      Ğœ      Üœ      èœ      ôœ                         $      0      <      H      T      `      l      x      „            œ      ¨      ´      À      Ì      Ø      ä      ğ      ü                         ,      8      D      P      \      h      t      €      Œ      ˜      ¤      °      ¼      È      Ô      à      ì      ø      Ÿ      Ÿ      Ÿ              ¨Ÿ      ­Ÿ              §«ª2                                                            §«ª2                                                                   »±°<                                            §«ª2                                                            §«ª2                                                            §«ª2                                                            §«ª2                                                            §«ª2                                                            §«ª2                                                            §«ª2                                                            §«ª2                                                            §«ª2                                                            §«ª2                                                            §«ª2                                                            §«ª2                                                            §«ª2                                                            §«ª2                                                            §«ª2                                                            §«ª2                                                            §«ª2                                                            §«ª2                                                            §«ª2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       (     0     8     @     H     P     X     `     h     p     x     €     ˆ          ˜           ¨     °     ¸     À     È     Ğ     Ø     à     è     ğ     ø                           (    0    8    @    H    P    X    `    h    p    x    €    ˆ        ˜         ¨    °    ¸    À    È    Ğ    Ø    à    è    ğ    ø                          (    0    H    @     @ÀğÀpĞğ  °àà°`àÀpÀğ°pp`ĞĞ ÀÀÀÀÀÀ€À€°À0 ĞĞĞĞĞpĞà à  ğ°0        d           +   d           8   f Á·U       .  À      x   $  À         $   @          N  @          .         Œ   $            $   À          N  À          .  À      Ÿ   $  À         $   ğ         N  ğ         .  °      «   $  °         $             N            .  @      ´   $  @         $   2	         N  2	      Á               Ï               à               é               ò                                d             d             d           !  f €·U       .  €      T  $  €         $   p          N  p          .  ğ      [  $  ğ         $   P         N  P         .  @      m  $  @         $   ğ         N  ğ         .  0!      w  $  0!         $   –         N  –         d             d           †  d           ‘  f €·U       .  Ğ$      Ä  $  Ğ$         $             N            .  ğ&      Ö  $  ğ&         $   (
         N  (
         d             d           à  d           é  f €·U       .   1        $   1         $   à          N  à          .   2      *  $   2         $   `         N  `         .  `8      2  $  `8         $   ¬         N  ¬         d             d           B  d           K  f ·U    |              Š              –              £              µ              Ç              Ù              ï              ú              	                            +              @              U              k                            ”              £              ²              Á              Ï              à              ñ                                          "              0              =              J              W              h              n              x                                                         °              Á              Ğ              à              ğ              ú                                                           d             d           #  d           .  f ·U       .  ;      a  $  ;         $   `          N  `          .  p;      g  $  p;         $   `         N  `         .  ĞA      q  $  ĞA         $   2         N  2         d             d           {  d           ‚  f ·U       .  D      ±  $  D         $   p          N  p          .  €D      ¾  $  €D         $   @         N  @         .  ÀE      Ç  $  ÀE         $   p         N  p         .  0G      Ñ  $  0G         $   0         N  0         .  `I      Ù  $  `I         $   p          N  p          .  ĞI      ã  $  ĞI         $   p          N  p          .  @J      î  $  @J         $   `          N  `          .   J      ú  $   J         $   Ğ          N  Ğ          .  pK      
  $  pK         $   P         N  P         .  ÀL        $  ÀL         $              N             .  `M      /  $  `M         $   À          N  À          .   N      D  $   N         $   À          N  À          .  àN      Z  $  àN         $   À          N  À          .   O      o  $   O         $   À          N  À          .  `P      …  $  `P         $   À          N  À          .   Q      œ  $   Q         $   À          N  À          .  àQ      ²  $  àQ         $             N            .  ğQ      Ã  $  ğQ         $   €          N  €          d             d           Ë  d           ×  f ·U       .  pR        $  pR         $   @         N  @         .  °S        $  °S         $   v         N  v         d             d           )  d           5  f ·U       .  0Z      i  $  0Z         $   0         N  0         .  `[      |  $  `[         $   À         N  À         .   b      ‡  $   b      •  „              $   #          N  #          d             d              d           ,  f ·U       .  Pb      `  $  Pb         $             N            .  pc      s  $  pc         $   Ï         N  Ï         d             d           ~  d           ‰  f ‚·U       .  @s      ¼  $  @s         $   Ğ          N  Ğ          .  t      Ì  $  t         $   P         N  P         .  `u      ß  $  `u         $   P         N  P         .  °v      ó  $  °v         $   P         N  P         d             d           	  d           	  f ‚·U       .   x      F	  $   x         $   …         N  …         d             d           W	  d           d	  f ‚·U       .  z      ™	  $  z         $   p          N  p          .   {      «	  $   {         $   P         N  P         .  P      ½	  $  P         $   R         N  R         d             d           Ñ	  d           Ü	  f ‚·U       .  °      
  $  °         $              N             .  P‚      
  $  P‚         $   Ø
         N  Ø
         d             d           /
  d           :
  f ‚·U       .  0      m
  $  0         $              N             .  Ğ      |
  $  Ğ         $            N           d             d           
  d           —
  f ‚·U       .  ğ“      É
  $  ğ“         $   ğ         N  ğ         .  à•      Ş
  $  à•         $   0         N  0         .  ™      ó
  $  ™         $   0          N  0          .  @™        $  @™         $            N           d                     ,    ´      ?    À      S           f    À      r    °      {    €      ‚    ğ      ”    @          0!      ­    Ğ$      ¿    ğ&      É     1      Ù     2      á    `8      ñ    ;      ÷    p;          ĞA          D          €D      !    ÀE      +    0G      3    `I      =    ĞI      H    @J      T     J      d    pK      s    ÀL      ‰    `M           N      ´    àN      É     O      ß    `P      ö     Q          àQ          ğQ      %    pR      8    °S      C    0Z      V    `[      a     b      o    Pb      ‚    pc          @s          t      °    `u      Ä    °v      ×     x      è    z      ú     {          P           °      /    P‚      @    0      O    Ğ      `    ğ“      u    à•      Š    ™      ˜    @™      ª    XÂ      ¸    ˜Â      É    ØÂ      ×    àÂ      ã    Ã      ğ    PÃ          Ã          ĞÃ      &    Ä      <    PÄ      G    Ä      V    ĞÄ      d    Å      x    PÅ          Å      ¢    ĞÅ      ¸    Æ      Ì    PÆ      á    Æ      ğ    ĞÆ      ÿ    Ç          PÇ          Ç      -    ĞÇ      6    ØÇ      ?    àÇ      Q    èÇ      c    ğÇ      t    ôÇ      …    øÇ      •    üÇ      ¥     È      ³    È      À    È      Í    È      Ú    È      ë    È      ñ     È      û    (È          0È          8È      #    @È      3    HÈ      D    PÈ      S    XÈ      c    `È      s    hÈ      }    ˆÈ      ‡    ˜È      ‘    ¨È      š    ÈÈ      ¦    @      ³            ¼            Ê            Ù            ë            ş                                    #            )            0            7            @            G            Q            W            ^            d            q            z            ‰            –            ¢            ¯            ·            Á            Ì            ×            ß            ï            ÷            ÿ            
                        #            4            J            d            n            †            ‘            ™            ¡            ±            Æ            Ö            æ            ô                                    ,            @            F            O            U            _            e            m            y                        ˆ            ”            œ            ¤            ±            º            Å            Ë            Ó            í  î  ï  ğ  ñ  ê  ë  ì  Ë  Ì  Í  Î  Ğ  Ñ  Ò  Ó  Ô  Õ  Ö  ×  Ø  Ù  Ú  Û  Ü  İ  Ş  ß  à  á  â  ã  ä  å  æ  ç  è  é  ò  ó  ô  õ  ö  ÷  ø  ù  ú  û  ü  ı  ş  ÿ                     	  
            Ï  í  î  ï  ğ  ñ  ê  ë  ì  Ë  Ì  Í  Î  Ğ  Ñ  Ò  Ó  Ô  Õ  Ö  ×  Ø  Ù  Ú  Û  Ü  İ  Ş  ß  à  á  â  ã  ä  å  æ  ç  è  é  ò  ó  ô  õ  ö  ÷  ø  ù  ú  û  ü  ı  ş  ÿ                     	  
              /Users/roboos/matlab/fieldtrip/peer/src/ memprofile.c /Users/roboos/matlab/fieldtrip/peer/src/../private/memprofile.o _memprofile_cleanup _memprofile_sample _memprofile _exitFun _mexFunction _mutexmemlist _mutexmemprofile _reftime _memlist _memprofileStatus _memprofileThread announce.c /Users/roboos/matlab/fieldtrip/peer/src/announce.o _frand _cleanup_announce _announce _announce_once discover.c /Users/roboos/matlab/fieldtrip/peer/src/discover.o _cleanup_discover _discover expire.c /Users/roboos/matlab/fieldtrip/peer/src/expire.o _cleanup_expire _expire _check_watchdog extern.c /Users/roboos/matlab/fieldtrip/peer/src/extern.o _syslog_level _condstatus _mutexstatus _mutexappendcount _mutexsocketcount _mutexthreadcount _mutexconnectioncount _mutexhost _mutexpeerlist _mutexjoblist _mutexallowuserlist _mutexrefuseuserlist _mutexallowgrouplist _mutexrefusegrouplist _mutexallowhostlist _mutexrefusehostlist _mutexwatchdog _mutexsmartmem _mutexsmartcpu _mutexprevcpu _mutexsmartshare _udsserverStatus _tcpserverStatus _announceStatus _discoverStatus _expireStatus _appendcount _socketcount _threadcount _connectioncount _host _peerlist _joblist _allowuserlist _refuseuserlist _allowgrouplist _refusegrouplist _allowhostlist _refusehostlist _smartsharelist _watchdog _smartmem _smartcpu _prevcpu _smartshare peerinit.c /Users/roboos/matlab/fieldtrip/peer/src/peerinit.o _hash _peerinit _peerexit util.c /Users/roboos/matlab/fieldtrip/peer/src/util.o _threadsleep _bufread _bufwrite _append _jobcount _peercount _hoststatus _clear_peerlist _clear_joblist _clear_smartsharelist _clear_allowuserlist _clear_allowgrouplist _clear_allowhostlist _clear_refuseuserlist _clear_refusegrouplist _clear_refusehostlist _check_datatypes _getmem udsserver.c /Users/roboos/matlab/fieldtrip/peer/src/udsserver.o _cleanup_udsserver _udsserver tcpserver.c /Users/roboos/matlab/fieldtrip/peer/src/tcpserver.o _cleanup_tcpserver _tcpserver __OSSwapInt16 /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/usr/include/libkern/i386/_OSByteOrder.h tcpsocket.c /Users/roboos/matlab/fieldtrip/peer/src/tcpsocket.o _cleanup_tcpsocket _tcpsocket security.c /Users/roboos/matlab/fieldtrip/peer/src/security.o _security_check _ismember_userlist _ismember_grouplist _ismember_hostlist localhost.c /Users/roboos/matlab/fieldtrip/peer/src/localhost.o _check_localhost smartshare.c /Users/roboos/matlab/fieldtrip/peer/src/smartshare.o _smartshare_reset _smartshare_check _smartshare_history smartmem.c /Users/roboos/matlab/fieldtrip/peer/src/smartmem.o _smartmem_info _smartmem_update smartcpu.c /U