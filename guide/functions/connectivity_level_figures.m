
%% HiGGS


cmap                = load(properties.general_params.colormap_path);
cmap_a              = cmap.cmap_a;
cmap_c              = cmap.cmap_c;
cmap                = cmap.cmap;


Sc                  = subject.Scortex;
Atlas               = Sc.Atlas(Sc.iAtlas).Scouts;

%%
%% Plotting results
%%


J                   = s2j;
sources_iv          = zeros(length(J),1);
sources_iv(indms)   = sqrt(abs(J(indms)));
sources_iv          = sources_iv/max(sources_iv(:));

figure_name = strcat('BC-VARETA-activation - ',str_band);

figure_BC_VARETA1 = figure('Color','w','Name',figure_name,'NumberTitle','off'); hold on;

define_ico(figure_BC_VARETA1);
patch('Faces',Sc.Faces,'Vertices',Sc.Vertices,'FaceVertexCData',sources_iv,'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
set(gca,'Color','w');
az = 0; el = 0;
view(az,el);
rotate3d on;
colormap(gca,cmap_a);
title('BC-VARETA-activation','Color','k','FontSize',16);

disp('-->> Saving figure');
file_name = strcat('BC_VARETA_activation','_',str_band,'.fig');
saveas(figure_BC_VARETA1,fullfile(pathname,file_name));

close(figure_BC_VARETA1);

pause(1e-12);

%%
%% Plotting results

temp_iv    = abs(s2j(indms));
connect_iv = abs(Thetajj);
temp       = abs(connect_iv);
temp_diag  = diag(diag(temp));
temp_ndiag = temp - temp_diag;
temp_ndiag = temp_ndiag/max(temp_ndiag(:));
temp_diag  = diag(temp_iv);
temp_diag  = temp_diag/max(temp_diag(:));
temp_diag  = diag(diag(temp_diag) + 1);
temp_comp  = temp_diag + temp_ndiag;
label_gen = [];
for ii = 1:length(indms)
    label_gen{ii} = num2str(ii);
end

figure_name = strcat('BC-VARETA-node-wise-conn - ',str_band);

figure_BC_VARETA2 = figure('Color','w','Name',figure_name,'NumberTitle','off');

define_ico(figure_BC_VARETA2);
imagesc(temp_comp);
set(gca,'Color','w','XColor','k','YColor','k','ZColor','k',...
    'XTick',1:length(indms),'YTick',1:length(indms),...
    'XTickLabel',label_gen,'XTickLabelRotation',90,...
    'YTickLabel',label_gen,'YTickLabelRotation',0);
xlabel('sources','Color','k');
ylabel('sources','Color','k');
colorbar;
colormap(gca,cmap_c);
axis square;
title('BC-VARETA-node-wise-conn','Color','k','FontSize',16);

disp('-->> Saving figure');
file_name = strcat('BC_VARETA_node_wise_conn','_',str_band,'.fig');
saveas(figure_BC_VARETA2,fullfile(pathname,file_name));

close(figure_BC_VARETA2);


%% Roi analysis
Thetajj_full              = zeros(length(Ke));
Sjj_full                  = zeros(length(Ke));
Thetajj_full(indms,indms) = Thetajj;
Sjj_full(indms,indms)     = diag(temp_iv);
atlas_label               = cell(1,length(Atlas));
conn_roi                  = zeros(length(Atlas));
act_roi                   = zeros(length(Atlas),1);
for roi1 = 1:length(Atlas)
    for roi2 = 1:length(Atlas)
        conn_tmp             = Thetajj_full(Atlas(roi1).Vertices,Atlas(roi2).Vertices);
        conn_tmp             = mean(abs(conn_tmp(:)));
        conn_roi(roi1,roi2)  = conn_tmp;
    end
    atlas_label{roi1} = Atlas(roi1).Label;
end
for roi1 = 1:length(Atlas)
    act_tmp              = diag(Sjj_full(Atlas(roi1).Vertices,Atlas(roi1).Vertices));
    act_tmp              = mean(abs(act_tmp));
    act_roi(roi1)        = act_tmp;
end
act_roi    = diag(act_roi);
temp_iv    = abs(act_roi);
connect_iv = abs(conn_roi);
temp       = abs(connect_iv);
temp_diag  = diag(diag(temp));
temp_ndiag = temp-temp_diag;
temp_ndiag = temp_ndiag/max(temp_ndiag(:));
temp_diag  = diag(abs(diag(temp_iv)));
temp_diag  = temp_diag/max(temp_diag(:));
temp_diag  = diag(diag(temp_diag)+1);
temp_comp  = temp_diag+temp_ndiag;
figure_name = strcat('BC-VARETA-roi-conn - ',str_band);

figure_BC_VARETA3 = figure('Color','w','Name',figure_name,'NumberTitle','off');

define_ico(figure_BC_VARETA3);
imagesc(temp_comp);
set(gca,'Color','w','XColor','k','YColor','k','ZColor','k',...
    'XTick',1:length(Atlas),'YTick',1:length(Atlas),...
    'XTickLabel',atlas_label,'XTickLabelRotation',90,...
    'YTickLabel',atlas_label,'YTickLabelRotation',0);
xlabel('sources','Color','k');
ylabel('sources','Color','k');
colorbar;
colormap(gca,cmap_c);
axis square;
title('BC-VARETA-roi-conn','Color','k','FontSize',16);

disp('-->> Saving figure');
file_name = strcat('BC_VARETA_roi_conn','_',str_band,'.fig');
saveas(figure_BC_VARETA3,fullfile(pathname,file_name));

close(figure_BC_VARETA3);


%% hg_LASSO
Lvj                 = subject.Ke;
Cdata               = subject.Cdata;
Sh                  = subject.Sh;
cmap                = load(properties.general_params.colormap_path);
cmap_a              = cmap.cmap_a;
cmap_c              = cmap.cmap_c;
Sc                  = subject.Sc;
Atlas               = Sc.Atlas(Sc.iAtlas).Scouts;
%%
%% Plotting results
temp_iv    = abs(Sjj);
connect_iv = abs(Thetajj);
temp       = abs(connect_iv);
temp_diag  = diag(diag(temp));
temp_ndiag = temp-temp_diag;
temp_ndiag = temp_ndiag/max(temp_ndiag(:));
temp_diag  = diag(abs(diag(temp_iv)));
temp_diag  = temp_diag/max(temp_diag(:));
temp_diag  = diag(diag(temp_diag)+1);
temp_comp  = temp_diag+temp_ndiag;
label_gen = [];
for ii = 1:length(indms)
    label_gen{ii} = num2str(ii);
end

figure_name = strcat('BC-VARETA-node-wise-conn - ',str_band);

figure_BC_VARETA2 = figure('Color','k','Name',figure_name,'NumberTitle','off');

define_ico(figure_BC_VARETA2);
imagesc(temp_comp);
set(gca,'Color','k','XColor','w','YColor','w','ZColor','w',...
    'XTick',1:length(indms),'YTick',1:length(indms),...
    'XTickLabel',label_gen,'XTickLabelRotation',90,...
    'YTickLabel',label_gen,'YTickLabelRotation',0);
xlabel('sources','Color','w');
ylabel('sources','Color','w');
colorbar;
colormap(gca,cmap_c);
axis square;
title('BC-VARETA-node-wise-conn','Color','w','FontSize',16);

disp('-->> Saving figure');
file_name = strcat('BC_VARETA_node_wise_conn','_',str_band,'.fig');
saveas(figure_BC_VARETA2,fullfile(pathname,file_name));

close(figure_BC_VARETA2);


%% Roi analysis
Thetajj_full              = zeros(length(Lvj)/3);
Sjj_full                  = zeros(length(Lvj)/3);
Thetajj_full(indms,indms) = Thetajj;
Sjj_full(indms,indms)     = Sjj;
atlas_label               = cell(1,length(Atlas));
conn_roi                  = zeros(length(Atlas));
act_roi                   = zeros(length(Atlas),1);
for roi1 = 1:length(Atlas)
    for roi2 = 1:length(Atlas)
        conn_tmp             = Thetajj_full(Atlas(roi1).Vertices,Atlas(roi2).Vertices);
        conn_tmp             = mean(abs(conn_tmp(:)));
        conn_roi(roi1,roi2)  = conn_tmp;
    end
    atlas_label{roi1} = Atlas(roi1).Label;
end

for roi1 = 1:length(Atlas)
    act_tmp              = diag(Sjj_full(Atlas(roi1).Vertices,Atlas(roi1).Vertices));
    act_tmp              = mean(abs(act_tmp));
    act_roi(roi1)        = act_tmp;
end
act_roi    = diag(act_roi);
temp_iv    = abs(act_roi);
connect_iv = abs(conn_roi);
temp       = abs(connect_iv);
temp_diag  = diag(diag(temp));
temp_ndiag = temp-temp_diag;
temp_ndiag = temp_ndiag/max(temp_ndiag(:));
temp_diag  = diag(abs(diag(temp_iv)));
temp_diag  = temp_diag/max(temp_diag(:));
temp_diag  = diag(diag(temp_diag)+1);
temp_comp  = temp_diag+temp_ndiag;
figure_name = strcat('BC-VARETA-roi-conn - ',str_band);

figure_BC_VARETA3 = figure('Color','k','Name',figure_name,'NumberTitle','off');

define_ico(figure_BC_VARETA3);
imagesc(temp_comp);
set(gca,'Color','k','XColor','w','YColor','w','ZColor','w',...
    'XTick',1:length(Atlas),'YTick',1:length(Atlas),...
    'XTickLabel',atlas_label,'XTickLabelRotation',90,...
    'YTickLabel',atlas_label,'YTickLabelRotation',0);
xlabel('sources','Color','w');
ylabel('sources','Color','w');
colorbar;
colormap(gca,cmap_c);
axis square;
title('BC-VARETA-roi-conn','Color','w','FontSize',16);

disp('-->> Saving figure');
file_name = strcat('BC_VARETA_roi_conn','_',str_band,'.fig');
saveas(figure_BC_VARETA3,fullfile(pathname,file_name));

close(figure_BC_VARETA3);
