function [RC,sv]=SsaYou(signal,m,nrec)
%% Singular spectrum analysis
%% input
% signal: signal
% m:      embedding dimension
% nrec:   number of reconstruction
%% output
% RC:     reconstructed components
%% main body
t=length(signal);
n=t-m+1;
WindowMatrix=hankel(signal(1:n),signal(n:end)); % embed
[u,s,v]=svd(WindowMatrix); % SVD
sv=diag(s); % singular values
% reconstruct
RC=zeros(nrec,t);
for i=1:nrec
    oneMatrix=u(:,i)*sv(i)*v(:,i)';
    buf=flipud(oneMatrix); % anti-diagonal mean
    for p=1:t
        RC(i,p)=mean(diag(buf,p-n));
    end
end
end