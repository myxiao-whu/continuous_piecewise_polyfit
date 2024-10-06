clear;clc;
%% import data
data=xlsread('global_basin_timeseries.xlsx',1);
tall=data(:,1);
ts=1;pd=0.9;ofac=20;
t=tall;h0=data(:,12);loc=find(isnan(h0));t(loc)=[];h0(loc)=[];n=length(t);
bpall=20:n-20;nbp=length(bpall);tbp=t(bpall);
%% time series: detrend1, detrend2
h=h0;
detrend1=detrend(h,1);detrend2=detrend(h,2);trend1=h-detrend1;trend2=h-detrend2;
% figure 1: plot initial time series
figure;plot(t,h,'k','linewidth',2.5);box on;grid on;axis tight;
xlabel('time/year','fontsize',14);ylabel('sea level/mm','fontsize',14);
% compare detrend1 and detrend2
figure;subplot(2,1,1);hold on;plot(t,h,'k','linewidth',1.5);box on;grid on;axis tight;
plot(t,trend1,'b',t,trend2,'r','linewidth',1.5);
xlabel('time/year','fontsize',14);ylabel('sea level/mm','fontsize',14);
legend('initial','linear trend','parabolic trend','fontsize',14,'location','best');
subplot(2,1,2);plot(t,detrend1,'b',t,detrend2,'r','linewidth',2.5);
box on;grid on;axis tight;
xlabel('time/year','fontsize',14);ylabel('sea level/mm','fontsize',14);
legend('linear detrend','parabolic detrend','fontsize',14,'location','best');
%% EMD
h=h0;
imf=emd(h,'t',t);figure;plot(t,h,t,imf(end,:));
%% SSA
h=h0;
[RC,sv]=SsaYou(h,40,3);
figure;subplot(3,1,1);plot(t,RC(1,:));
subplot(3,1,2);plot(t,RC(2,:));
subplot(3,1,3);plot(t,RC(3,:));
figure;subplot(2,1,1);plot(t,h,t,sum(RC(1:1,:),1));axis tight;
subplot(2,1,2);plot(t,h,t,sum(RC(1:2,:),1));axis tight;
%% piecewise linear fit -> break point at 1953? no no no
h=h0;
vall=nan(nbp,1);
for i=1:nbp
    vall(i)=continuous_piecewise_polyfit(t,h,bpall(i),[1 1],[]);
end
figure;plot(tbp,vall);vplot1=vall;
[~,loc]=min(vall);disp([bpall(loc) t(bpall(loc))]) % 54 -> 1953
[~,~,~,hpiece_linear]=continuous_piecewise_polyfit(t,h,54,[1 1],[]);
figure;plot(t,h,'k',t,hpiece_linear,'b','linewidth',1.5);box on;grid on;axis tight;
xlabel('time/year','fontsize',14);ylabel('sea level/mm','fontsize',14);
legend('initial','piecewise linear trend','fontsize',14,'location','best');
%% break point seem to be 1962, but simple fit bad
h=h0;
hbp=[detrend(h(t<1962));detrend(h(t>=1962),2)];
figure;subplot(2,1,1);plot(t,hbp);axis tight;
subplot(2,1,2);plot(t,h,t,h-hbp);axis tight;
%% piece wise linear fit at 1962: still bad
h=h0;
[~,~,~,fitvalue]=continuous_piecewise_polyfit(t,h,63,[1 1],[]);
figure;plot(t,h,'k',t,fitvalue,'b','linewidth',1.5);box on;grid on;axis tight;
xlabel('time/year','fontsize',14);ylabel('sea level/mm','fontsize',14);
legend('initial','piecewise linear trend','fontsize',14,'location','best');
%% piece wise linear+parabolic fit: 1959? but I still prefer 1962
h=h0;
vall=nan(nbp,1);
for i=1:nbp
    vall(i)=continuous_piecewise_polyfit(t-t(1),h,bpall(i),[1 2],[]);
end
vplot2=vall;
figure;plot(t(bpall),vall,'k');box on;grid on;axis tight;set(gca,'fontsize',12)
xlabel('break point','fontsize',16);ylabel('fit error','fontsize',16);
[~,loc]=min(vall);disp([bpall(loc) t(bpall(loc))]) % 60 -> 1959
%% piece wise linear+parabolic fit at 1962:
h=h0;
[~,~,~,hnocycle]=continuous_piecewise_polyfit(t,h,63,[1 2],[]);hdetrend=h-hnocycle;
tiledlayout(2,1,'padding','tight','tilespacing','tight');
nexttile;hold on;plot(t,h,'color',[1 1 1]*0.8,'linewidth',7.6);
plot(t,hnocycle,'k','linewidth',1.6);set(gca,'fontsize',12);axis tight;box on;grid on;
xlabel('time/year','fontsize',16);ylabel('sea level from glacier/mm','fontsize',16);
legend('initial data','piecewise fit','fontsize',18,'location','best');
nexttile;plot(t,hdetrend,'k',t,hbp,'b','linewidth',1.5);axis tight;box on;grid on;
xlabel('time/year','fontsize',14);ylabel('sea level/mm','fontsize',14);
legend('piecewise continous','piecewise discontinous','fontsize',14,'location','best');
%% compare with simple linear or parabolic trend
tiledlayout(2,2,'padding','tight','tilespacing','tight');
nexttile;plot(t,h,'k',t,trend1,'b','linewidth',1.5);axis tight;box on;grid on;
xlabel('time/year','fontsize',14);ylabel('sea level/mm','fontsize',14);
legend('initial','linear trend','fontsize',14,'location','best');
nexttile;plot(t,h,'k',t,trend2,'r','linewidth',1.5);axis tight;box on;grid on;
xlabel('time/year','fontsize',14);ylabel('sea level/mm','fontsize',14);
legend('initial','parabolic trend','fontsize',14,'location','best');
nexttile;plot(t,h,'k',t,hnocycle,'g','linewidth',1.5);axis tight;box on;grid on;
xlabel('time/year','fontsize',14);ylabel('sea level/mm','fontsize',14);
legend('initial','piecewise trend','fontsize',14,'location','best');
nexttile;plot(t,detrend1,'b',t,detrend2,'r',t,hdetrend,'g','linewidth',2.5);
xlabel('time/year','fontsize',14);ylabel('sea level/mm','fontsize',14);
legend('linear','parabolic','piecewise','fontsize',14,'location','best');
axis tight;box on;grid on;
set(gcf,'position',[0 0 700 350]);
%% compare with simple linear, parabolic, SSA trend
% time series
tiledlayout(5,2,'padding','tight','tilespacing','tight');
h=h0;colork=[1 1 1]*0.8;colorjet=jet(12);i=9;linew=3.6;
nexttile(1,[2 1]);plot(t,h,'color',colork,'linewidth',linew);hold on;
plot(t,trend1,'b',t,trend2,'r',t,RC(1,:),'c','linewidth',1.5);
text(1910,-8,'(a)','fontsize',16);axis tight;box on;grid on;
set(gca,'xticklabel',[]);ylabel('sea level/mm','fontsize',14);
legend('raw data','linear trend','parabolic trend','SSA trend','fontsize',12,'location','best','box','off');
nexttile(2,[2 1]);plot(t,h,'color',colork,'linewidth',linew);hold on;
plot(t,hnocycle,'k','linewidth',1.5);axis tight;box on;grid on;
plot(t,hpiece_linear,'color',colorjet(i,:),'linewidth',1.5);plot(t,h-hbp,'g','linewidth',1.5)
set(gca,'xticklabel',[]);
legend('raw data','trend of CPPFH','trend of CPLF','trend of disc.t. p.w. fit','fontsize',12,'location','best','box','off');
text(1910,-10,'(b)','fontsize',16);
nexttile(5,[3 1]);plot(t,detrend1,'b',t,detrend2,'r',t,h'-RC(1,:),'c','linewidth',2.5);
xlabel('time/year','fontsize',14);ylabel('sea level/mm','fontsize',14);
legend('detrend: linear','detrend: parabolic','detrend: SSA','fontsize',12,'location','south','box','off');
axis tight;box on;grid on;ylim([-8 max(h'-RC(1,:))]);
text(1910,8,'(c)','fontsize',16);
nexttile(6,[3 1]);plot(t,hdetrend,'k','linewidth',1.5);hold on;
plot(t,h-hpiece_linear,'color',colorjet(i,:),'linewidth',1.5);plot(t,hbp,'g','linewidth',1.5);
xlabel('time/year','fontsize',14);axis tight;box on;grid on;ylim([-3 max(h-hpiece_linear)])
legend('detrend: CPPFH','detrend: CPLF','detrend: discontinuous piecewise fit','fontsize',12,'location','south','box','off');
text(1910,1.67,'(d)','fontsize',16);
%% detrended spectrum
plomb(hdetrend,t,[],20);xlim([0 50]);
cwt(hdetrend,'frequencylimits',[0.01 0.1]);
%% piecewise fit with simple last cycles
h=h0;
vall=nan(nbp,1);fitperiod=[52.3 27.3];
for i=1:nbp
    vall(i)=continuous_piecewise_polyfit(t-t(1),h,bpall(i),[1 2],fitperiod);
end
[val,loc]=min(vall);vall2=vall;disp([val bpall(loc) t(bpall(loc))])
figure;plot(t(bpall),vall,'k');box on;grid on;axis tight;set(gca,'fontsize',12)
xlabel('break point','fontsize',16);ylabel('fit error','fontsize',16);title(num2str(val));
% fit effect
[~,~,~,fitvalue]=continuous_piecewise_polyfit(t,h,bpall(loc),[1 2],fitperiod);herr=h-fitvalue;
figure;hold on;plot(t,h,'color',[1 1 1]*0.8,'linewidth',7.6);
plot(t,fitvalue,'k','linewidth',1.6);set(gca,'fontsize',12);axis tight;box on;grid on;
xlabel('time/year','fontsize',16);ylabel('sea level from glacier/mm','fontsize',16);
legend('initial data','piecewise fit','fontsize',18,'location','best');
%% fix at 1962, test exact period
h=h0;
p1all=15:0.1:35;n1=length(p1all);
p2all=45:0.1:60;n2=length(p2all);
vall=nan(n1,n2);
for i=1:n1
    for j=1:n2
        fitp=unique([p1all(i) p2all(j)]);
        vall(i,j)=continuous_piecewise_polyfit(t,h,63,[1 2],fitp);
    end
end
[val,loc]=min(vall(:));disp(val)
[i,j]=ind2sub(size(vall),loc);disp([i j p1all(i) p2all(j)])
% local min
bwmin=imregionalmin(vall);[i,j]=ind2sub(size(vall),find(bwmin==true));disp([i j p1all(i)' p2all(j)'])
figure;imagesc(p2all,p1all,vall);colorbar;hold on;
scatter(p2all(j),p1all(i),[],'r');
text(p2all(j)+0.25,p1all(i)+0.5,num2str([p2all(j)' p1all(i)'],'%.1f+%.1f'));
%% find the best combination of change point and periods
h=h0;loc1all=55:70;nloc=length(loc1all);vall=nan(n1,n2,nloc);
for i=1:n1
    for j=1:n2
        for k=1:nloc
            fitp=unique([p1all(i) p2all(j)]);
            vall(i,j,k)=continuous_piecewise_polyfit(t,h,loc1all(k),[1 2],fitp);
        end
    end
end
% min
[x,y,z]=meshgrid(p1all,p2all,loc1all);
n123=n1*n2*nloc;
xre=reshape(x,n123,1);
yre=reshape(y,n123,1);
zre=reshape(z,n123,1);
vre=reshape(vall,n123,1);
figure;scatter3(xre,yre,zre,[],vre);axis tight;c=colorbar;c.Label.String='拟合误差';
xlabel('周期1');ylabel('周期2');zlabel('突变点');
% find minimal
[val,loc]=min(vall(:));disp(val)
[i,j,k]=ind2sub(size(vall),loc);
disp([i j k]);
disp([p1all(i) p2all(j) loc1all(k)]); % 68 30.4 63, 
% for each combination of periods, determined change point keeps the same?
locdet=nan(n1,n2);vdet=nan(n1,n2);
for i=1:n1
    for j=1:n2
        vone=vall(i,j,:);
        [val,loc]=min(vone);locdet(i,j)=loc1all(loc);vdet(i,j)=val;
    end
end
figure;subplot(2,1,1);imagesc(p1all,p2all,locdet');colorbar;
subplot(2,1,2);imagesc(p1all,p2all,vdet');colorbar;
%% plot for paper about piecewise fit
h=h0;
[~,solution,A,fitvalue]=continuous_piecewise_polyfit(t,h,63,[1 2],[52.3 27.3]);
htrend=A(:,1:4)*solution(1:4);seasonal=A(:,5:end)*solution(5:end);hdetrend=h-htrend;
tiledlayout(3,1,'tilespacing','tight','padding','tight');
nexttile;plot(t,h,t,htrend);
nexttile;plot(t,h,t,fitvalue);
nexttile;plot(t,seasonal);
%% detrended spectrum
plomb(hdetrend,t,[],20);xlim([0 50]);
cwt(hdetrend,'frequencylimits',[0.01 0.1]);
%% plot for paper : show process of real data
tiledlayout(2,3,'TileSpacing','tight','Padding','tight');h=hdetrend;
nexttile(1);plot(t,h0,'k',t,h0-hdetrend-4,'b','linewidth',1.3);axis tight;
xlabel('time/year','fontsize',14);ylabel('sea level/mm','fontsize',14);
text(1905,-1,'(a)','fontsize',16);
legend('initial data','piecewise fitting','location','southeast','fontsize',14)
nexttile(2);plot(t,hdetrend,'k','linewidth',1.3);axis tight;
xlabel('time/year','fontsize',14);ylabel('sea level/mm','fontsize',14);text(1905,0.75,'(b)','fontsize',16);
nexttile(4);[p,f,pth]=plomb(h,t,[],ofac,'pd',0.9);
plot(f,p,'k','linewidth',1.6);box on;hold on;
yline(pth,'k--','linewidth',1.6);xlim([0 0.15]);
xlabel('frequency/cpy','fontsize',14);ylabel('amplitude/mm','fontsize',14);
[~,f0,amp0]=find_target_peak(f,p,[50;25]);scatter(f0,amp0,'r','filled');
text(f0+0.003,amp0,num2str(1./f0,'%.0f'),'color','r');
text(0.004,0.65,'(d)','fontsize',16);legend('spectrum','90% confidence level','fontsize',13,'location','best');
nexttile(3);[wt, freqs] = cwt(h);imagesc(t, freqs, abs(wt));axis xy;
nexttile(6);scatter3(xre,yre,t(zre),[],vre);axis tight;c=colorbar;c.Label.String='fitting error/mm';c.Label.FontSize=14;
xlabel('period 1/year','fontsize',14);ylabel('period 2/year','fontsize',14);zlabel('breakpoint/year','fontsize',14);
text(16,60,1972,'(f)','fontsize',16);
nexttile(5);plot(tbp,vall2,'k',t(63),vall2(44),'ro','linewidth',1.3);axis tight;text(1922,3.05,'(e)','fontsize',16);grid on;
xlabel('time of breakpoint/year','fontsize',14);ylabel('fitting error/mm','fontsize',14);text(1970,0.5,'1962','color','r');
set(gcf,'position',[100 100 900 450]);
%% plot for paper: time series trend and detrend
h=h0;fig=tiledlayout(5,5);colork=[1 1 1]*0.8;colorjet=jet(12);i=9;linew=3.6;
nexttile(1,[2 2]);plot(t,h,'color',colork,'linewidth',linew);hold on;
plot(t,trend1,'b',t,trend2,'r',t,RC(1,:),'c','linewidth',1.5);
text(1910,-9,'(a)','fontsize',16);axis tight;box on;grid on;
set(gca,'xticklabel',[]);ylabel('sea level/mm','fontsize',14);
legend('raw data','linear trend','parabolic trend','SSA trend','fontsize',12,'location','best','box','off');
nexttile(3,[2 2]);plot(t,h,'color',colork,'linewidth',linew);hold on;
plot(t,hpiece_linear,'color',colorjet(i,:),'linewidth',1.5);plot(t,hnocycle,'k','linewidth',1.5);axis tight;box on;grid on;
set(gca,'xticklabel',[]);
legend('raw data','trend: cont. p.w. l. fit','trend: our algorithm','fontsize',12,'location','best','box','off');
text(1910,-10,'(b)','fontsize',16);
nexttile(11,[3 2]);plot(t,detrend1,'b',t,detrend2,'r',t,h'-RC(1,:),'c','linewidth',2.5);
xlabel('time/year','fontsize',14);ylabel('residue and fitting error/mm','fontsize',14);
legend('detrend: linear','detrend: parabolic','detrend: SSA','fontsize',12,'location','south','box','off');
axis tight;box on;grid on;ylim([-8 max(h'-RC(1,:))]);
text(1910,8,'(d)','fontsize',16);
nexttile(13,[3 2]);plot(t,h-hpiece_linear,'color',colorjet(i,:),'linewidth',1.5);hold on;plot(t,herr,'k','linewidth',1.5);
xlabel('time/year','fontsize',14);axis tight;box on;grid on;ylim([-3 max(h-hpiece_linear)])
legend('detrend: continuous piecewise linear fit','detrend: our proposed algorithm','fontsize',12,'location','south','box','off');
text(1910,1.5,'(e)','fontsize',16);
nexttile(5,[2 1]);plot(tbp,h0(bpall),'color',colork,'linewidth',linew);axis tight;set(gca,'xticklabel',[]);
xline(1953,'--','color',colorjet(i,:),'linewidth',1.5);xline(1962,'--','color','k','linewidth',1.5);hold on;
plot(t(54),h0(54),'o','color',colorjet(i,:),'linewidth',1.5);plot(t(63),h0(63),'ko','linewidth',1.5);
text(1930,-18,'(c)','fontsize',16);
nexttile(15,[3 1]);plot(tbp,vplot1,'color',colorjet(i,:),'linewidth',2.5);hold on;
plot(tbp,vall2,'color','k','linewidth',2.5);xlabel('time/year','fontsize',14);axis tight;
plot(t(54),vplot1(35),'o','color',colorjet(i,:),'linewidth',1.5);plot(t(63),vall2(44),'ko','linewidth',1.5);
xline(1953,'--','color',colorjet(i,:),'linewidth',1.5);xline(1962,'--','color','k','linewidth',1.5);
text(1930,3.05,'(f)','fontsize',16);
fig.TileSpacing = 'tight';
fig.Padding = 'compact';
set(gcf,'position',[0 0 800 350]);
