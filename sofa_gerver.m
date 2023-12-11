% --- Basic data

xpg=-sqrt(2)/2;
xkg=sqrt(2)/2;
dxg=(xkg-xpg)/100;

xpd=2*xpg;
xkd=2*xkg;
dxd=(xkd-xpd)/100;
 
xg=xpg:dxg:xkg;
xd=xpd:dxd:xkd;

yg=abs(xg)+sqrt(2);
yd=abs(xd);

% --- Gerver's sofa path and angle

[xtp,ytp,fitp]=gerver_tor(2500,6000); 

% --- Sofa generation

dxtdttp=diff(xtp);
dytdttp=diff(ytp);
dfidttp=diff(fitp);

xt=xtp(1:end-1);
yt=ytp(1:end-1);
fit=fitp(1:end-1);

dytdtt=dytdttp;
dxtdtt=dxtdttp;
dfidtt=dfidttp;

% --- Up-center
xNN=0;
yNN=sqrt(2);
xNsGC=cos(fit).*(xt + xNN) + sin(fit).*(yt - yNN);
yNsGC=sin(fit).*(xt + xNN) - cos(fit).*(yt - yNN);

% --- Up-left
xNN=-(dxtdtt - dytdtt - 2^(1/2)*dfidtt + dfidtt.*xt + dfidtt.*yt)./(2*dfidtt);
yNN=-xNN+sqrt(2);
xNsGL=cos(fit).*(xt + xNN) + sin(fit).*(yt - yNN);
yNsGL=sin(fit).*(xt + xNN) - cos(fit).*(yt - yNN);

% --- Down-left
xNN=-(dxtdtt - dytdtt + dfidtt.*xt + dfidtt.*yt)./(2*dfidtt);
yNN=-xNN;
xNsDL=cos(fit).*(xt + xNN) + sin(fit).*(yt - yNN);
yNsDL=sin(fit).*(xt + xNN) - cos(fit).*(yt - yNN);

% --- Postprocessing
ind=isnan(xNsGL);
xNsGL(ind)=[];
yNsGL(ind)=[];
ind=isnan(yNsGL);
xNsGL(ind)=[];
yNsGL(ind)=[];

ind=isnan(xNsDL);
xNsDL(ind)=[];
yNsDL(ind)=[];
ind=isnan(yNsDL);
xNsDL(ind)=[];
yNsDL(ind)=[];

ind=isnan(xNsGC);
xNsGC(ind)=[];
yNsGC(ind)=[];
ind=isnan(yNsGC);
xNsGC(ind)=[];
yNsGC(ind)=[];

XDL=xNsDL;
YDL=yNsDL;
XGL=xNsGL;
YGL=yNsGL;

GLimx=[-10,0];
GLimy=[0,0];
P = InterX([xNsDL;yNsDL],[GLimx;GLimy]);
xC1=P(1);%xNsDL(1);%
yC1=P(2);%yNsDL(1);
ind=find(yNsDL>0);
xNsDL(ind)=[];
yNsDL(ind)=[];
xNsDL=([xC1,xNsDL]);
yNsDL=([yC1,yNsDL]);

[x0,y0,segments]=selfintersect(xNsDL,yNsDL);
P=segments;
if isempty(P)~=1
xNsDL(P(1,1):P(1,2))=[];
yNsDL(P(1,1):P(1,2))=[];
xNsDL(P(1,1)-1)=x0(1);
yNsDL(P(1,1)-1)=y0(1);
else
end

xNsDL=[xNsDL,0];
yNsDL=[yNsDL,-1];

SLimx=[0,0];
SLimy=[-2,5];
P = InterX([xNsGC;yNsGC],[SLimx;SLimy]);
xC4=P(1);
yC4=P(2);
ind=find(xNsGC>=0);
xNsGC(ind)=[];
yNsGC(ind)=[];
xNsGC=fliplr([xNsGC,xC4]);
yNsGC=fliplr([yNsGC,yC4]);

XGC=xNsGC;
YGC=yNsGC;

P = InterX([xNsGL;yNsGL],[GLimx;GLimy]);
xC6=P(1);
yC6=P(2);
ind=find(yNsGL>0);
xNsGL(ind)=[];
yNsGL(ind)=[];
xNsGL=([xNsGL,xC6]);
yNsGL=([yNsGL,yC6]);

[P,i,j] = InterX2([xNsGC;yNsGC],[xNsGL;yNsGL]);
xNsGC=xNsGC(1:i);
yNsGC=yNsGC(1:i);
xC5=P(1,1);
yC5=P(2,1);
xNsGC=[xNsGC,xC5];
yNsGC=[yNsGC,yC5];

xNsGL=xNsGL(j+1:end);
yNsGL=yNsGL(j+1:end);

xNs=[xNsDL,xNsGC,xNsGL];
yNs=[yNsDL,yNsGC,yNsGL];

% --- Area

sofaGerv=polyshape(xNs,yNs);
A=2*area(sofaGerv);

GERV=2.2195316688;
AG=num2str(A, 12)

figure
plot(sofaGerv)
axis equal; grid on

save sofaGerv.mat sofaGerv

sofaGerv2=polyshape([xNs,-xNs],[yNs,yNs]);
save sofaGerv2.mat sofaGerv2

% --- Visualisation

xpg=-2*sqrt(2);
xkg=2*sqrt(2);
dxg=(xkg-xpg)/100;

xpd=-2*sqrt(2);
xkd=2*sqrt(2);
dxd=(xkd-xpd)/100;
 
xg=xpg:dxg:xkg;
xd=xpd:dxd:xkd;

yg=abs(xg)+sqrt(2);
yd=abs(xd);

figure
for i=1:50:length(xtp)-1
xsd=cos(fitp(i))*(xtp(i) + xd) + sin(fitp(i))*(ytp(i) - yd);
ysd=sin(fitp(i))*(xtp(i) + xd) - cos(fitp(i))*(ytp(i) - yd);
xsg=cos(fitp(i))*(xtp(i) + xg) + sin(fitp(i))*(ytp(i) - yg);
ysg=sin(fitp(i))*(xtp(i) + xg) - cos(fitp(i))*(ytp(i) - yg);

plot(xsd,ysd,'-b')
hold on
plot(xsg,ysg,'-r')
hold on
end
plot(sofaGerv)
axis equal; grid on 

figure
subplot(2,1,1)
plot(XDL,YDL,'-.k','Linewidth',1.5);
hold on
plot(XGC,YGC,':k','Linewidth',1.5);
hold on
plot(XGL,YGL,'-k','Linewidth',1.5);
axis equal; grid on
axis tight
yticks(-1:0.25:0.1)
ylim([-1 0.1]);
xticks(-1.8:1.8/4:0)
xlim([-1.8 0])
xlabel('$x_s$','Interpreter','latex')
ylabel('$y_s$','Interpreter','latex')
legend('$\bar{r}_{sl}^{(s)}$','$\bar{r}_{scu}^{(s)}$','$\bar{r}_{su}^{(s)}$','Interpreter','latex')

subplot(2,1,2)
plot(sofaGerv,'Linewidth',1.5)
axis equal; grid on 
axis tight
yticks(-1:0.25:0.1)
ylim([-1 0.1]);
xticks(-1.8:1.8/4:0)
xlim([-1.8 0])
xlabel('$x_s$','Interpreter','latex')
ylabel('$y_s$','Interpreter','latex')
box on

figure
subplot(2,2,1)
plot(xtp,ytp,'-k','Linewidth',1.5);
grid on
yticks(1.8:0.1:2.1)
xticks(-0.5:0.25:0.5)
xlim([-0.5 0.5])
xlabel('$x_{t}^{(h)}$','Interpreter','latex')
ylabel('$y_{t}^{(h)}$','Interpreter','latex')
title('(a)','FontWeight','normal')

subplot(2,2,3)
plot(xtp,fitp,'-k','Linewidth',1.5);
grid on
yticks(-1:0.5:1)
xticks(-0.5:0.25:0.5)
xlim([-0.5 0.5])
ylabel('$\varphi , rad$','Interpreter','latex')
xlabel('$x_{t}^{(h)}$','Interpreter','latex')
title('(b)','FontWeight','normal')

subplot(2,2,2)
plot(XDL,YDL,'-.k','Linewidth',1.5);
hold on
plot(XGC,YGC,':k','Linewidth',1.5);
hold on
plot(XGL,YGL,'-k','Linewidth',1.5);
axis equal; grid on
axis tight
yticks(-1.2:0.2:0.2)
ylim([-1.2 0.2]);
xticks(-1.8:1.8/4:0)
xlim([-1.8 0])
xlabel('$x_s$','Interpreter','latex')
ylabel('$y_s$','Interpreter','latex')
legend('$\bar{r}_{sl}^{(s)}$','$\bar{r}_{scu}^{(s)}$','$\bar{r}_{su}^{(s)}$','Interpreter','latex','Location','northwest')
title('(c)','FontWeight','normal')

subplot(2,2,4)
plot(sofaGerv,'Linewidth',1.5)
axis equal; grid on 
axis tight
yticks(-1.2:0.2:0.2)
ylim([-1.2 0.2]);
xticks(-1.8:1.8/4:0)
xlim([-1.8 0])
xlabel('$x_s$','Interpreter','latex')
ylabel('$y_s$','Interpreter','latex')
title('(d)','FontWeight','normal')
box on
