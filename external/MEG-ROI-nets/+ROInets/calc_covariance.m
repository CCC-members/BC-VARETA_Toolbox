function C= calc_covariance(y)
[m,n,s]=size(y);%m:variabel,n:length,s:segments or trials
C=reshape(cell2mat(arrayfun(@(i) (squeeze(y(:,:,i))*squeeze(y(:,:,i).'))/n,1:s,'UniformOutput',false)),m,m,s);
C=mean(C,3);
end






