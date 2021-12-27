function [q1pp, q2pp]=modelomanipulador(tao1,tao2,q1,q1p,q2,q2p)
//Parámetros del robot
l1=0.45; lc1=0.091; l2=0.45; lc2=0.048;
m1=23.902; m2=3.88; I1=1.266; I2=0.093;
b1=2.288; b2=0.175; fc1=7.17; fc2=1.734;
g=9.81;d1=3.23;d2=1.2;
//Matriz de inercia
M11= m1*lc1*lc1 + m2*l1*l1+m2*lc2*lc2 +2*m2*l1*lc2*cos(q2)+I1+I2;
M12=m2*lc2*lc2+m2*l1*lc2*cos(q2)+I2;
M21=m2*lc2*lc2+m2*l1*lc2*cos(q2)+I2;
M22=m2*lc2*lc2+I2;
//Matriz de Coriolis
C11=-m2*l1*lc2*sin(q2)*q2p;
C12=-m2*l1*lc2*sin(q2)*(q1p+q2p);
C21=m2*l1*lc2*sin(q2)*q1p;
C22=0;
// Matriz de fuerzas hidrodinámicas 
D11=d1*m1*sin(q2)*q2p;
D12=d1*m1
D21=d2*m2**cos(q2)*q1p;
D22=d2*m2
//Par gravitacional
g1=(m1*lc1+m2*l1)*g*sin(q1) + m2*g*lc2*sin(q1+q2);
g2=m2*g*lc2*sin(q1+q2);
//Par de fricción
//fric1=fc1*sign(q1p) + b1*q1p;
//fric2=fc2*sign(q2p) + b2*q2p;
fric1=b1*q1p;
fric2=b2*q2p;

//Ecuaciones diferenciales lazo cerrado
detM=M11*M22-M21*M12;
va1=tao1-C11*q1p-C12*q2p-g1-fric1;
vb1=tao2-C21*q1p-C22*q2p-g2-fric2;
vc1=(M22*va1 - M12*vb1)/detM;
vd1=(-M21*va1 + M11*vb1)/detM;
q1pp=vc1;
q2pp=vd1;
endfunction
