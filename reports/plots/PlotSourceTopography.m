function PlotSourceTopography(reportPath,freq,activationData,Cortex)

if(~isfolder(reportPath))
    mkdir(reportPath);
end

template                = load('axes.mat');
currentAxes             = template.axes;
colorMap                = load('tools/mycolormap_brain_basic_conn.mat');
J                       = activationData.J;
sources_iv              = sqrt(abs(J));
sources_iv              = sources_iv/max(sources_iv(:));
smoothValue             = 0.66;
SurfSmoothIterations    = 10;
Vertices                = tess_smooth(Cortex.Vertices, smoothValue, SurfSmoothIterations, Cortex.VertConn, 1);

patch(currentAxes, ...
    'Faces',Cortex.Faces, ...
    'Vertices',Vertices, ...
    'FaceVertexCData',Cortex.SulciMap*0.06+log(1+sources_iv), ...
    'FaceColor','interp', ...
    'EdgeColor','none', ...
    'AlphaDataMapping', 'none', ...
    'EdgeColor',        'none', ...
    'EdgeAlpha',        1, ...
    'BackfaceLighting', 'lit', ...
    'AmbientStrength',  0.5, ...
    'DiffuseStrength',  0.5, ...
    'SpecularStrength', 0.2, ...
    'SpecularExponent', 1, ...
    'SpecularColorReflectance', 0.5, ...
    'FaceLighting',     'gouraud', ...
    'EdgeLighting',     'gouraud', ...
    'FaceAlpha',1);
set(currentAxes,'xcolor','w','ycolor','w','zcolor','w');
colormap(currentAxes,colorMap.cmap_a);

fig = figure("Name","Source Model","Color","w", 'Position', get(0, 'Screensize'),'Visible','off');
set(currentAxes,"Parent",fig);
rotate3d(currentAxes,'on');
axis(currentAxes,'equal');
axis(currentAxes,'off');

% right view
view(currentAxes,0,0);
export_fig(fullfile(reportPath,strcat('sourceTopography_',freq.str_band,'_right')),'-transparent','-png');
% left view
view(currentAxes,180,0);
export_fig(fullfile(reportPath,strcat('sourceTopography_',freq.str_band,'_left')),'-transparent','-png');
% top view
view(currentAxes,0,90);
export_fig(fullfile(reportPath,strcat('sourceTopography_',freq.str_band,'_top')),'-transparent','-png');
% bottom view
view(currentAxes,0,270);
export_fig(fullfile(reportPath,strcat('sourceTopography_',freq.str_band,'_bottom')),'-transparent','-png');
close(fig);
end
