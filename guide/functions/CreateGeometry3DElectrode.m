function ElectrodeGrid = CreateGeometry3DElectrode(Channel, ChanLoc)
    
    % Initialize returned values   
    ElectrodeGrid  = [];
    % Get subject
    
    % ===== CONTACTS GEOMETRY =====    
    % Define electrode geometry
    nVert = 32;
    [sphereVertex, sphereFaces] = tess_sphere(nVert);    
    FaceLighting = 'flat';
   
    % ===== DISPLAY SEEG/ECOG ELECTRODES =====
    iChanProcessed = [];
    UserData    = [];
    Vertex      = [];
    Faces       = [];
    VertexAlpha = [];
    VertexRGB   = [];   
    
    % ===== ADD SPHERE CONTACTS ======
    % Get the sensors that haven't been displayed yet
    iChanOther = setdiff(1:length(Channel), iChanProcessed);
    isValidLoc = ~any(all(ChanLoc(iChanOther,:)==0,2),1);
    % Display spheres
    if ~isempty(iChanOther) && isValidLoc
        % Get the saved display defaults for this modality

        ctSize    = [1 1 1] .* 0.01 ./ 2;
        tmpVertex = sphereVertex;
        tmpFaces  = sphereFaces;
        ctOrient  = [];
        % Force Gouraud lighting
        FaceLighting = 'gouraud';
        
        % Create contacts geometry
        [ctVertex, ctFaces] = Plot3DContacts(tmpVertex, tmpFaces, ctSize, ChanLoc(iChanOther,:), ctOrient);
        % Display properties
        ctColor   = [.9,.9,0];  % YELLOW
        elecAlpha = 1;
        % Add to global patch
        offsetVert  = size(Vertex,1);
        Vertex      = [Vertex;      ctVertex];
        Faces       = [Faces;       ctFaces + offsetVert];
        VertexAlpha = [VertexAlpha; repmat(elecAlpha, size(ctVertex,1), 1)];
        VertexRGB   = [VertexRGB;   repmat(ctColor,   size(ctVertex,1), 1)];
        % Save the channel index in the UserData
        UserData    = [UserData;    reshape(repmat(iChanOther, size(ctVertex,1)./length(iChanOther), 1), [], 1)];
    end
    % Create patch
    if ~isempty(Vertex)
        ElectrodeGrid.Faces               = Faces;
        ElectrodeGrid.Vertices            = Vertex;
        ElectrodeGrid.FaceVertexCData     = VertexRGB;
        ElectrodeGrid.FaceVertexAlphaData = VertexAlpha;
        ElectrodeGrid.Options = {...
            'EdgeColor',        'none', ...
            'BackfaceLighting', 'unlit', ...
            'AmbientStrength',  0.5, ...
            'DiffuseStrength',  0.6, ...
            'SpecularStrength', 0, ...
            'FaceLighting',     FaceLighting, ...
            'Tag',              'ElectrodeGrid', ...
            'UserData',         UserData};
    end
end

