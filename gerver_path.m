function [xi,yi,fii] = gerver_path(n,m)
  fi=0.03918;
  teta=0.68130;
    
  dalf1=fi/n;
  alf1=0:dalf1:fi;
  
  dalf2=(teta-fi)/n;
  alf2=fi+dalf2:dalf2:teta;
  
  dalf3=(pi/2-teta-teta)/n;
  alf3=teta+dalf3:dalf3:pi/2-teta;
  
  dalf4=(pi/2-fi-(pi/2-teta))/n;
  alf4=pi/2-teta+dalf4:dalf4:pi/2-fi;
  
  dalf5=fi/n;
  alf5=pi/2-fi+dalf5:dalf5:pi/2;
  
  a1=1.210322422072688751;
  a2=-1/4;
  b1=-0.527624598026784624;
  b2=0.920258385160637622;
  c1=0.626045522848465867;
  c2=-0.944750803946430751;
  d1=1.313022761424232933;
  d2=-0.525382670414554437;
  e1=1.210322422072688751;
  e2=1/4;
  
  k11=-0.210322422072688751;
  k12=1/4;
  k21=-0.919179292771593322;
  k22=0.472406619750805465;
  k31=-0.613763229430251668;
  k32=0.889626479003221860;
  k41=-0.308347166088910014;
  k42=0.472406619750805465;
  k51=-1.017204036787814585;
  k52=1/4;
  
  k1=[k11;k12];
  k2=[k21;k22];
  k3=[k31;k32];
  k4=[k41;k42];
  k5=[k51;k52];
  
  i=1;
  for t1=0:dalf1:fi
    w1=Rt(t1)*[a1*cos(t1)+a2*sin(t1)-1;-a2*cos(t1)+a1*sin(t1)-0.5]+k1;
    x1(i)=w1(1);
    y1(i)=w1(2);
    i=i+1;
  end
  
  i=1;
  for t2=fi+dalf2:dalf2:teta
    w2=Rt(t2)*[-0.25*t2.^2+b1*t2+b2;0.5*t2-b1-1]+k2;
    x2(i)=w2(1);
    y2(i)=w2(2);
    i=i+1;
  end
  
  i=1;
  for t3=teta+dalf3:dalf3:pi/2-teta
    w3=Rt(t3)*[c1-t3;c2+t3]+k3;
    x3(i)=w3(1);
    y3(i)=w3(2);
    i=i+1;
  end
  
  i=1;
  for t4=pi/2-teta+dalf4:dalf4:pi/2-fi
    w4=Rt(t4)*[-0.5*t4+d1-1;-0.25*t4.^2+d1*t4+d2]+k4;
    x4(i)=w4(1);
    y4(i)=w4(2);
    i=i+1;
  end
  
  i=1;
  for t5=pi/2-fi+dalf5:dalf5:pi/2
    w5=Rt(t5)*[e1*cos(t5)+e2*sin(t5)-0.5;-e2*cos(t5)+e1*sin(t5)-1]+k5;
    x5(i)=w5(1);
    y5(i)=w5(2);
    i=i+1;
  end
  
  x=[x1 x2 x3 x4 x5];
  y=[y1 y2 y3 y4 y5];
  alf=[alf1 alf2 alf3 alf4 alf5];
  fitp=alf;
  
  xo=x(end)/2;
  yo=0;
  for i=1:length(x)
    xt(i)=- (xo - x(i))*((2^(1/2)*sin(fitp(i)))/2 + (2^(1/2)*cos(fitp(i)))/2) - (yo - y(i))*((2^(1/2)*sin(fitp(i)))/2 - (2^(1/2)*cos(fitp(i)))/2);
    yt(i)=2^(1/2) + (xo - x(i))*((2^(1/2)*sin(fitp(i)))/2 - (2^(1/2)*cos(fitp(i)))/2) - (yo - y(i))*((2^(1/2)*sin(fitp(i)))/2 + (2^(1/2)*cos(fitp(i)))/2);
  end
  fit=-(fitp-pi/4);
  
  xt=fliplr(xt);
  yt=fliplr(yt);
  fit=fliplr(fit);
  
  dxt=(xt(end)-xt(1))/m;
  xi=xt(1):dxt:xt(end);
  fii=spline(xt,fit,xi);
  yi=spline(xt,yt,xi);
end

