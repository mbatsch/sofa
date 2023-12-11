[xt,yt,fi]=gerver_path(20,50);

n_droga=10;
xtp=[linspace(-1.5,xt(1),n_droga),xt,linspace(xt(end),1.5,n_droga)];
ytp=[abs(linspace(-1.5,xt(1),n_droga))+sqrt(2),yt,linspace(xt(end),1.5,n_droga)+sqrt(2)];
fitp=[-ones(1,n_droga)*pi/4, fi, ones(1,n_droga)*pi/4];

xtp = [xtp, fliplr(xtp)];
ytp = [ytp, fliplr(ytp)];
fitp = [fitp, fliplr(fitp)];

load('sofaGerv2.mat')
sofaX=sofaGerv2.Vertices(:,1);
sofaY=sofaGerv2.Vertices(:,2);
sofa=polyshape(sofaX,sofaY);

dr=figure;
dr.Color='white';

plot(-3:0.01:3,abs(-3:0.01:3)+sqrt(2),'-black','LineWidth',1)
hold on
plot(-3:0.01:3,abs(-3:0.01:3),'-black','LineWidth',1)
hold on
ax = gca;
ax.XLim = [-3 3];
ax.YLim = [0 5];
xlabel 'x'
ylabel 'y'
axis equal
ew = animatedline('Color','b','LineWidth',1);
pr = animatedline;
axis tight

loops=length(xtp);

F(loops) = struct('cdata',[],'colormap',[]);
filename = 'sofa_move_Gerver.gif';
hold on
    
for i=1:2:length(xtp)

xsN=sofaX*cos(fitp(i)) - xtp(i) + sofaY*sin(fitp(i));
ysN=ytp(i) + sofaY*cos(fitp(i)) - sofaX*sin(fitp(i));
if i==1
    sofa=fill(xsN,ysN,[0 0.4470 0.7410],'LineWidth',1);
    drawnow limitrate
else
    sofa.XData=xsN;
    sofa.YData=ysN;
    drawnow limitrate
end

F(i) = getframe;
frame = getframe(rys);
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
if i == 1 
 imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
else 
 imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
end 
end
