function [ThetaRR] = source_to_region_map(Thetajj,Scouts,Regions,Nsources)

% Initialize R with zeros
R = zeros(length(Regions), Nsources);
for i = 1:length(Regions)
    region = Regions{i};
    ScoutsRegion_ind = ismember({Scouts.Region},{region});
    ScoutsRegion = Scouts(ScoutsRegion_ind);
    for r=1:length(ScoutsRegion)
        R(i, ScoutsRegion(r).Vertices) = 1;
    end
end
for i = 1:size(R, 1)
    norm_i = norm(R(i, :));  % Calculate the norm of the row
    if norm_i > 0
        R(i, :) = R(i, :) / norm_i;  % Normalize the row
    end
end

ThetaRR = R*Thetajj*R';
end



