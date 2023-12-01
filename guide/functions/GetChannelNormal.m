%% ===== GET CHANNEL NORMALS =====
% USAGE: GetChannelNormal(sSubject, ChanLoc, Modality)
%        GetChannelNormal(sSubject, ChanLoc, SurfaceType)   % SurfaceType={'scalp','innerskull','cortex','cortexhull','cortexmask'}
function [ChanOrient, ChanLocProj] = GetChannelNormal(Scalp, ChanLoc, SurfaceType, isInteractive)
    % Initialize returned variables
    ChanOrient = [];
    ChanLocProj = ChanLoc;
    % CALL: GetChannelNormal(sSubject, ChanLoc, Modality)
    if ismember(SurfaceType, {'EEG','NIRS'})
        SurfaceType = 'scalp';
    elseif ismember(SurfaceType, {'SEEG','ECOG','ECOG+SEEG'})
        if ~isempty(sSubject.iInnerSkull)
            SurfaceType = 'innerskull';
        elseif ~isempty(sSubject.iCortex)
            SurfaceType = 'cortexmask';
        elseif ~isempty(sSubject.iScalp)
            SurfaceType = 'scalp';
        else
            SurfaceType = 'cortexmask';
        end
    end
    % Get surface
    isConvhull = 0;
    isMask = 0;    
    % Load surface (or get from memory)
    if isMask
        % Compute surface based on MRI mask
        % [sSurf, sOrig] = tess_envelope(SurfaceFile, 'mask_cortex', 5000);
        sSurf = bst_memory('GetSurfaceEnvelope', SurfaceFile, 5000, 0, 1); 
        Vertices    = sSurf.Vertices;
        VertNormals = tess_normals(sSurf.Vertices, sSurf.Faces);
    else
        % Load surface        
        Vertices    = Scalp.Vertices;
        VertNormals = Scalp.VertNormals;
        % Get convex hull
        if isConvhull
            Faces = convhulln(Vertices);
            iVertices = unique(Faces(:));
            Vertices    = Vertices(iVertices, :);
            VertNormals = VertNormals(iVertices, :);
        end
    end
    
    % Project electrodes on the surface 
    ChanLocProj = channel_project_scalp(Vertices, ChanLoc);
    % Get the closest vertex for each channel
    iChanVert = bst_nearest(Vertices, ChanLocProj,1,0);
    % Get the normals at those points
    ChanOrient = VertNormals(iChanVert, :);
% view_surface_matrix(Vertices, sSurf.Faces);

% OTHER OPTIONS WITH SPHERICAL HARMONICS, A BIT FASTER, NOT WORKING WELL
%     % Compute spherical harmonics
%     fvh = hsdig2fv(Vertices, 20, 5/1000, 40*pi/180, 0);
%     VertNormals = tess_normals(fvh.vertices, fvh.faces);
%     % Get the closest vertex for each channel
%     iChanVert = bst_nearest(fvh.vertices, ChanLoc);
%     % Get the normals at those points
%     ChanOrient = VertNormals(iChanVert, :);
%     % Project electrodes on the surface 
%     ChanLocProj = channel_project_scalp(fvh.vertices, ChanLoc);

end
