function dataFilt = filt_data(W,Fs,f,var,varargin)
%   FILT_DATA Summary of this function goes here
% 'FourierCoefficients',W,'Fs',Fs,'Freq',F(freq),'varf',varf,'use_gpu',use_gpu
%
for i=1:2:length(varargin)
    eval([varargin{i} '=  varargin{(i+1)};'])
end

if(~exist('use_gpu','var'))
    use_gpu = false;
end
Nf              = size(W,2);
deltaf          = Fs/(Nf/2);
F               = 0:deltaf:(Nf-1)*deltaf;
gaussian_window = exp(-(F-f).^2/(2*var^2))/(sqrt(2*pi)*var) + exp(-(F-(F(end)-f)).^2/(2*var^2))/(sqrt(2*pi)*var);
if (use_gpu)
    W               = gpuArray(W).*repmat(gpuArray(gaussian_window),size(W,1),1,size(W,3));
    dataFilt        = ifft(W,[],2);
    dataFilt        = gather(real(dataFilt(:,1:(Nf/2),:)));
else
    W               = W.*repmat(gaussian_window,size(W,1),1,size(W,3));
    dataFilt        = ifft(W,[],2);
    dataFilt        = real(dataFilt(:,1:(Nf/2),:));
end
end