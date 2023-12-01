function ElectrodeGrid = CreateGeometry3DElectrode(ChanLoc,Modality,Scalp,varargin)

%% Optional params
%   nVert           -number of vertices for each electrode, default: 32
%   ctColor         - sphere color. default: [.9,.9,0];  % YELLOW
%   elecAlpha       - alectrode alpha. default: 1
%   sizeScale       - size scale for electrodes. default: 2
%   Modality        - EEG | MEG | ECoG 
%
%

for k = 4 : nargin
    eval([inputname(k) '=  varargin{k-3};']);
end

% Define electrode geometry
if(~exist('nVert','var'))
    nVert = 32;
end
if(~exist('ctColor','var'))
    ctColor   = [.9,.9,0];  % YELLOW
end

% Initialize returned values
ElectrodeGrid  = [];
% Get subject

% ===== CONTACTS GEOMETRY =====
% SEEG contact cylinder
[seegVertex, seegFaces] = tess_cylinder(nVert, 0.5);
% ECOG contact cylinder: Define electrode geometry (double-layer for Matlab < 2014b)
% if (bst_get('MatlabVersion') < 804)
%     nVert = 66;
%     [ecogVertex, ecogFaces] = tess_cylinder(nVert, 0.8, [], [], 1);
% else
    nVert = 34;
    [ecogVertex, ecogFaces] = tess_cylinder(nVert, 0.8, [], [], 0);
% end
% Define electrode geometry
nVert = 32;
[sphereVertex, sphereFaces] = tess_sphere(nVert);
% Get display configuration from iEEG tab
% Optimal lighting depends on Matlab version
% if (bst_get('MatlabVersion') < 804)
%     FaceLighting = 'gouraud';
% else
    FaceLighting = 'flat';
% end
% Compute contact normals: ECOG and EEG
if (ismember(Modality, {'ECOG','EEG'}) || (~isempty(sElectrodes) && any(strcmpi({sElectrodes.Type}, 'ECOG'))))
    ChanNormal = GetChannelNormal(Scalp, ChanLoc, Modality, 0);
else
    ChanNormal = [];
end

% ===== DISPLAY SEEG/ECOG ELECTRODES =====
iChanProcessed = [];
UserData    = [];
Vertex      = [];
Faces       = [];
VertexAlpha = [];
VertexRGB   = [];

% ===== ADD SPHERE CONTACTS ======
    % Get the sensors that haven't been displayed yet
    iChanOther = setdiff(1:size(ChanLoc,1), iChanProcessed);
    isValidLoc = ~any(all(ChanLoc(iChanOther,:)==0,2),1);
    % Display spheres
    if ~isempty(iChanOther) && isValidLoc
        % Get the saved display defaults for this modality
        ElectrodeConfig.Type            = 'eeg';
        ElectrodeConfig.ContactDiameter = 0.010;
        ElectrodeConfig.ContactLength   = 0.002;
        ElectrodeConfig.ElecDiameter    = [];
        ElectrodeConfig.ElecLength      = [];
        % SEEG: Sphere
        if strcmpi(Modality, 'SEEG')
            ctSize    = [1 1 1] .* ElectrodeConfig.ContactLength ./ 2;
            tmpVertex = sphereVertex;
            tmpFaces  = sphereFaces;
            ctOrient  = [];
            % Force Gouraud lighting
            FaceLighting = 'gouraud';
        % ECOG/EEG: Cylinder (if normals are available)
        elseif ~isempty(ChanNormal)
            ctSize   = [ElectrodeConfig.ContactDiameter ./ 2, ElectrodeConfig.ContactDiameter ./ 2, ElectrodeConfig.ContactLength];
            ctOrient = ChanNormal(iChanOther,:);
            tmpVertex = ecogVertex;
            tmpFaces  = ecogFaces;
        % ECOG/EEG: Sphere (if normals are not available)
        else
            ctSize    = [1 1 1] .* ElectrodeConfig.ContactDiameter ./ 2;
            tmpVertex = sphereVertex;
            tmpFaces  = sphereFaces;
            ctOrient  = [];
            % Force Gouraud lighting
            FaceLighting = 'gouraud';
        end
        % Create contacts geometry
        [ctVertex, ctFaces] = Plot3DContacts(tmpVertex, tmpFaces, ctSize, ChanLoc(iChanOther,:), ctOrient);
        % Display properties
        
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

