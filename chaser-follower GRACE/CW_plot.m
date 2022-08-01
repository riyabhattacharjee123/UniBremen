close all
clear 
mu =  
n=0.00113;

phi_rr = @(t) [4-3*cos(n*t) 0 0;
            6*(sin(n*t)-n*t) 1 0; 
            0 0 cos(n*t)];
phi_rv = @(t) [1/n*sin(n*t) 2/n*(1-cos(n*t)) 0; 
               2/n*(cos(n*t)-1) 1/n*(4*sin(n*t)-3*n*t) 0;
               0 0 1/n*sin(n*t)];
           
phi_vr = @(t) [3*n*sin(n*t) 0 2;
                6*n*(cos(n*t)-1) 0 0;
                0 0 n*(-1)*sin(n*t)];

phi_vv = @(t) [cos(n*t) 2*sin(n*t) 0;
               (-1)*2*sin(n*t) (4*cos(n*t)-3) 0;
               0 0 cos(n*t)];

r_0 = [6741000;0;0];
v_0 = [0;7689;0];

T=0:10:100;

% plot red circle at origin
plot3(0,0,0,'ro')
hold on

% start loop
for t=T
    % display 50% progres
    if t==T(floor(end/2))
        disp('50 % done');
    end
    r=phi_rr(t)*r_0 + phi_rv(t)*v_0;
    disp(r);
    
    v=phi_vr(t)*r_0 + phi_vv(t)*v_0;
    disp(v);
    
    plot3(r(1),r(2),r(3),'k.');
    
end
rotate3d
xlabel('x'), ylabel('y'), zlabel('z')
axis equal
grid on
hold off