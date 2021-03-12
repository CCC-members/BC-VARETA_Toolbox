function [subject,properties] = get_connectivity_priors(subject,properties)
analysis_method = properties.connectivity_params.methods{1};
fields          = fieldnames(analysis_method);
method_name     = fields{1};

if(analysis_method.(method_name).run)
disp('=================================================================');
disp('BC-V-->> Getting connectivity priors.');

Ke                      = subject.Ke;
Sc                      = subject.Sc;
activation_params       = properties.activation_params;
connectivity_params     = properties.connectivity_params;
aSulc                   = connectivity_params.aSulc.value; % baseline of sulci curvature factor
aGiri                   = connectivity_params.aGiri.value; % baseline of giri curvature factor
bSulc                   = connectivity_params.bSulc.value; % scale of sulci curvature factor
bGiri                   = connectivity_params.bGiri.value; % scale of giri curvature factor
IsCurv                  = connectivity_params.IsCurv.value; % 0 (no compensation) 1 (giri and sulci curvature compensation)
IsNeigh                 = connectivity_params.IsNeigh.value;
IsField                 = connectivity_params.IsField.value; % 1 (projected Lead Field) 3 (3D Lead Field)
IsCurv_act              = activation_params.IsCurv.value; % 0 (no compensation) 1 (giri and sulci curvature compensation)
IsNeigh_act             = activation_params.IsNeigh.value;
IsField_act             = activation_params.IsField.value; % 1 (projected Lead Field) 3 (3D Lead Field)
GridOrient              = subject.GridOrient;
GridAtlas               = subject.GridAtlas;
Faces                   = Sc.Faces; 

%%
%% neigh/field options
%%
if (IsNeigh == IsNeigh_act) && (IsField == IsField_act) && (isfield(subject,'W')) && (isfield(subject,'Winv'))
    disp('-->> Creating Laplacian & Normals');
else
    disp('-->> Creating Laplacian & Normals');
    regLaplacian    = connectivity_params.regLaplacian.value;
    [D,D3D]         = graph_laplacian(Faces,regLaplacian);
    Dinv            = inv(D);
    Dinv            = (Dinv + Dinv)/2;
    D3Dinv          = inv(D3D);
    D3Dinv          = (D3Dinv + D3Dinv)/2;
    if (~IsNeigh)
        if IsField == 1
            W       = speye(length(D));
            Winv    = speye(length(D));
        elseif IsField == 2
            Ninv    = blk_diag(GridOrient',1);
            W       = Ninv;
            Winv    = W';
        elseif IsField == 3
            W      = speye(length(D3D));
            Winv   = speye(length(D3D));
        end
    elseif (IsNeigh)
        if IsField == 1
            W       = Dinv;
            Winv    = D;
        elseif IsField == 2
            Ninv    = blk_diag(GridOrient',1);
            DNinv   = Ninv*Dinv;
            W       = DNinv;
            Winv    = D*Ninv';
        elseif IsField == 3
            W       = D3Dinv;
            Winv    = D3D;
        end
    end
    subject.W       = W;
    subject.Winv    = Winv;
end
%%
%% curv/field options
%%
if (IsCurv == IsCurv_act) && (IsField == IsField_act)
    disp('-->> Creating curvature compensator');
else
    disp('-->> Creating curvature compensator');
    if IsField == 1
        Ke                    = bst_gain_orient(Ke, GridOrient,GridAtlas);
    end   
    if IsCurv == 1
        Curv                  = Sc.Curvature;
        Sulc                  = Sc.SulciMap;
        Curv                  = abs(Curv);
        CurvSulc              = zeros(length(Curv),1);
        CurvGiri              = zeros(length(Curv),1);
        CurvSulc(Sulc == 1)   = aSulc + bSulc.*Curv(Sulc == 1);
        CurvSulc(Sulc == 0)   = 1;
        CurvGiri(Sulc == 0)   = aGiri + bGiri.*Curv(Sulc == 0);
        CurvGiri(Sulc == 1)   = 1;
        if IsField == 1
            Ke_giri               = Ke.*repmat(CurvGiri',size(Ke,1),1);
            Ke_sulc               = Ke.*repmat(CurvSulc',size(Ke,1),1);
        elseif IsField == 2 || IsField == 3
            Sulc3D                = zeros(1,3*length(Sulc));
            CurvSulc3D            = zeros(1,3*length(Curv));
            CurvGiri3D            = zeros(1,3*length(Curv));
            node3 = 1;
            for node = 1:length(Curv)
                CurvSulc3D([node3 node3+1 node3+2]) = repmat(CurvSulc(node),1,3);
                CurvGiri3D([node3 node3+1 node3+2]) = repmat(CurvGiri(node),1,3);
                Sulc3D([node3 node3+1 node3+2])     = repmat(Sulc(node),1,3);
                node3                               = node3 + 3;
            end
            Ke_giri               = Ke.*repmat(CurvGiri3D,size(Ke,1),1);
            Ke_sulc               = Ke.*repmat(CurvSulc3D,size(Ke,1),1);
        end
        subject.Ke_giri = Ke_giri;
        subject.Ke_sulc = Ke_sulc;
    end
    subject.Ke = Ke;
end
end
end

