function [v,solution,A,fitvalue]=continuous_piecewise_polyfit(t,h,bp,order,pall)
%% continuous piecewise polynomial fitting with/without cycles
%% input
% t        : time ( N )
% h        : signal ( N )
% bp       : location of break point ( l )
% order    : polynomial order of each segment ( l+1 )
% pall     : fitting periods
%% output
% v        : fit error
% solution : solution of parameters of piecewise linear functions
% A        : coefficient matrix
% fitvalue : fit value
%% main body
t=t(:);h=h(:);bp=bp(:);order=order(:);
t=t-t(1); % to avoid too large number
N=length(t);l=length(bp);
if(ismember([1;N],bp))
    error('break point can not be on boundary');
end
bp=[1;bp;N]; % break point
b=zeros(l+2,1);b(2:end)=t(bp(2:end)); % time of break point
% record location of deletion
maxorder=max(order);
loc_order=nan(1,maxorder*(l+1));
for k=1:(l+1) % each part from 1,2,3,...,l,l+1
    loc_order((k-1)*maxorder+(1:order(k)))=1;
end
loc_del=isnan(loc_order);
% matrix A
A=zeros(N,maxorder*(l+1)+1);
for k=1:l+1 % l+1 segment
    for i=bp(k):bp(k+1) % belonging rows
        for j=1:k-1 % the first 1,2,3,...,k-1 column block
            loc=(j-1)*maxorder;
            for n=1:maxorder
                A(i,loc+n)=b(j+1)^n-b(j)^n;
            end
        end
        loc=(k-1)*maxorder; % the last k column block
        for n=1:maxorder
            A(i,loc+n)=t(i)^n-b(k)^n;
        end
        A(i,end)=1; % the last is always 1, corresponding to intercept
    end
end
A(:,loc_del)=[]; % delete those unused parameters
if(nargin>4)
    tspan=t(end)-t(1);
    pall(pall>tspan)=[]; % delete impossible period
    npall=length(pall);
    for i=1:npall
        wt=2*pi*t/pall(i);
        A=[A cos(wt) sin(wt)];
    end
end
solution=(A'*A)\(A'*h);
fitvalue=A*solution;
v=norm(h-fitvalue)/sqrt(N-size(A,2)); % sigma=V'*P*V/r
end