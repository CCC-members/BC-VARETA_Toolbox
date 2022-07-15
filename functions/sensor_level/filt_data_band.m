function dataFilt = filt_data_band(W,Fs,Fmin,Fmax,var,varargin)
%   FILT_DATA Summary of this function goes here
% 'FourierCoefficients',W,'Fs',Fs,'Freq',F(freq),'varf',varf,'use_gpu',use_gpu
%
for i=1:2:length(varargin)
    eval([varargin{i} '=  varargin{(i+1)};'])
end

if(~exist('use_gpu','var'))
    use_gpu = false;
end
Nf                     = size(W,2);
deltaf                 = Fs/(Nf/2);
F                      = 0:deltaf:(Nf-1)*deltaf;
band                   = find(F == Fmin):find(F == Fmax);
band_window            = exp(-(F-Fmin).^2/(2*var^2)) + ...
                         exp(-(F-(F(end)-Fmin)).^2/(2*var^2)) + ...
                         exp(-(F-Fmax).^2/(2*var^2)) + ...
                         exp(-(F-(F(end)-Fmax)).^2/(2*var^2));
band_window(band)      = 1;
band_window(Nf - band) = 1;
band_window            = band_window/sum(band_window);
if (use_gpu)
    W               = gpuArray(W).*repmat(gpuArray(band_window),size(W,1),1,size(W,3));
    dataFilt        = ifft(W,[],2);
    dataFilt        = gather(real(dataFilt(:,1:(Nf/2),:)));
else
    W               = W.*repmat(band_window,size(W,1),1,size(W,3));
    dataFilt        = ifft(W,[],2);
    dataFilt        = real(dataFilt(:,1:(Nf/2),:));
end
end
