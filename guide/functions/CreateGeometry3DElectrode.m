function ElectrodeGrid = CreateGeometry3DElectrode(Channel, ChanLoc, varargin)

%% Optional params
%   nVert           -number of vertices for each electrode, default: 32
%   ctColor         - sphere color. default: [.9,.9,0];  % YELLOW
%   elecAlpha       - alectrode alpha. default: 1
%   sizeScale       - size scale for electrodes. default: 2
%
%

for k = 3 : nargin
    eval([inputname(k) '=  varargin{k-2};']);
end
% Initialize returned values
ElectrodeGrid  = [];
% Get subject

% ===== CONTACTS GEOMETRY =====
% Define electrode geometry
if(~exist('nVert','var'))
    nVert = 32;
end
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
    if(~exist('sizeScale','var'))
        sizeScale = 2;
    end
    ctSize    = [1 1 1] .* 0.01 ./ sizeScale;
    tmpVertex = sphereVertex;
    tmpFaces  = sphereFaces;
    ctOrient  = [];
    % Force Gouraud lighting
    FaceLighting = 'gouraud';

    % Create contacts geometry
    [ctVertex, ctFaces] = Plot3DContacts(tmpVertex, tmpFaces, ctSize, ChanLoc(iChanOther,:), ctOrient);
    % Display properties
    if(~exist('ctColor','var'))
        ctColor   = [.9,.9,0];  % YELLOW
    end
    if(~exist('elecAlpha','var'))
        elecAlpha = 1;
    end
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

