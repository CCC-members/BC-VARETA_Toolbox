function [Vertex, Faces] = Plot3DContacts(ctVertex, ctFaces, ctSize, ChanLoc, ChanOrient)
    % Apply contact size
    ctVertex = bst_bsxfun(@times, ctVertex, ctSize);
    % Duplicate this contact
    nChan  = size(ChanLoc,1);
    nVert  = size(ctVertex,1);
    nFace  = size(ctFaces,1);
    Vertex = zeros(nChan*nVert, 3);
    Faces  = zeros(nChan*nFace, 3);
    for iChan = 1:nChan
        % Apply orientation
        if ~isempty(ChanOrient) && ~isequal(ChanOrient(iChan,:), [0 0 1])
            v1 = [0;0;1];
            v2 = ChanOrient(iChan,:)';
            % Rotation matrix (Rodrigues formula)
            angle = acos(v1'*v2);
            axis  = cross(v1,v2) / norm(cross(v1,v2));
            axis_skewed = [ 0 -axis(3) axis(2) ; axis(3) 0 -axis(1) ; -axis(2) axis(1) 0];
            R = eye(3) + sin(angle)*axis_skewed + (1-cos(angle))*axis_skewed*axis_skewed;
            % Apply rotation to the vertices of the electrode
            ctVertexOrient = ctVertex * R';
        else
            ctVertexOrient = ctVertex;
        end
        % Set electrode position
        ctVertexOrient = bst_bsxfun(@plus, ChanLoc(iChan,:), ctVertexOrient);
        % Report in final patch
        iVert  = (iChan-1) * nVert + (1:nVert);
        iFace = (iChan-1) * nFace + (1:nFace);
        Vertex(iVert,:) = ctVertexOrient;
        Faces(iFace,:)  = ctFaces + nVert*(iChan-1);
    end
end

