t = 6000 % seconds
n = 0.00115697 % rad/sec
nt = n*t;


phi_rr = [ (4 - 3*cos(nt)) (0) (0); ...
    (6*(sin(nt)-nt)) (1) (0); ...
    (0) (0) (cos(nt))];

disp(phi_rr);

phi_rv = [ ((1/n)*sin(nt)) (0) (0);
    ((2/n)*(cos(nt)-1)) ((1/n)*(4*sin(nt)-3*nt)) (0);
    (0) (0) ((1/n)*sin(nt)) ];

disp(phi_rv);

phi_vr = [ (3*n*sin(nt)) (0) (0);
    (6*n*(cos(nt)-1)) (0) (0);
    (0) (0) (-n*sin(nt)) ];

disp(phi_vr);

phi_vv =[(cos(nt)) (2*sin(nt)) (0);
    (-2*sin(nt)) (4*cos(nt)-3) (0);
    (0) (0) (cos(nt))];
disp(phi_vv);

