function [Svv_complex_r,K_r] = applying_reference(Svv_complex,K)
H  = eye(size(Svv_complex,1))-ones(size(Svv_complex,1))/size(Svv_complex,1);
K_r = H*K;
Svv_complex_r = Svv_complex;
for ii = 1:size(Svv_complex,3)
    Svv_complex_r(:,:,ii) = H*squeeze(Svv_complex(:,:,ii))*H;
    Svv_complex_r(:,:,ii) = (Svv_complex_r(:,:,ii) + Svv_complex_r(:,:,ii)')/2;
end
end