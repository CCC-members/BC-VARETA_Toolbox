function [rth] = check_likelihood_trend(llh0,llh_grid,rth_grid,ntry)
llh0       = repmat(llh0,length(rth_grid),1);
diff_llh   = llh_grid - [llh0 llh_grid(:,1:(ntry-1))];
id_grow    = zeros(length(rth_grid),ntry);
id_grow(diff_llh > 0) = 1;
for count_try = flip(1:ntry)
    llh_test     = llh_grid(:,1);
    id_grow_tmp  = prod(id_grow(:,1:count_try),2);
    id_grow_tmp  = find(id_grow_tmp > 0);
    if sum(id_grow_tmp) > 0
        llh_test_grow  = llh_test(id_grow_tmp);
        id_rth         = intersect(find(llh_test == max(llh_test_grow)),id_grow_tmp);
        id_rth         = min(id_rth);
        rth            = rth_grid(id_rth);
        break
    elseif (sum(id_grow_tmp) == 0) && (count_try == 1)
        rth          = 0;
    end
end
end