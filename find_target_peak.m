function [loc_target,f0,amp0]=find_target_peak(f,p,target_period,varargin)
%% find the nearest peak to target period
% usually for low-frequency signal
% if input is peak itselt, then simply choose the nearest
f=f(:);p=p(:);target_period=target_period(:);
n_target=length(target_period);
loc_target=nan(n_target,1);
if(nargin<3)
    error('not enough input');
end
if(nargin==3) % findpeaks
    [peak_amp,loc_peak]=findpeaks(p);
    peak_f=f(loc_peak);peak_period=1./peak_f;
    loc=nan(n_target,1);
    for i=1:n_target
        [~,loc(i)]=min(abs(peak_period-target_period(i)));
    end
    f0=peak_f(loc);amp0=peak_amp(loc);
    loc_target=loc_peak(loc);
    return;
end
% nargin>3
if(isa(varargin{1},'double')) % the fourth input is pth
    pth=varargin{1};    
    [peak_amp,loc_peak]=findpeaks(p);peak_f=f(loc_peak);
    locsig=find(peak_amp>pth);
    if(isempty(locsig))
        loc_target=[];f0=[];amp0=[];
        return;
    end
    fsig=peak_f(locsig);psig=peak_amp(locsig);
    period_sig=1./fsig;
    loc=nan(n_target,1);
    for i=1:n_target
        [~,loc(i)]=min(abs(period_sig-target_period(i)));
    end
    f0=fsig(loc);amp0=psig(loc);
    loc_target=loc_peak(locsig(loc));
else % the fourth input is not pth, which means just p is peaks
    periods=1./f;
    for i=1:n_target
        [~,loc_target(i)]=min(abs(periods-target_period(i)));
    end
    f0=f(loc_target);
    amp0=p(loc_target);  
end
