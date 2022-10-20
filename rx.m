% 01/17/2022 10:20:36 UTC (MM/DD/YYYY HH:MM:SS)
utc = [2022 1 17 10 20 36];

% Receiver Antenna locations in Lat., Lon., Alt.

rcvr_1=[53.10657108482692, 8.850392209543823 18]; % Bremen Uni 53°N , 8.8°E, 18m 
rcvr_2=[53.106565449092145, 8.850392209543823 18]; 
rcvr_3=[53.10656021590728, 8.850392209543823 18]; 
rcvr_4=[53.10655176231102, 8.850392209543823 18];

rcvr_5=[53.10657108482692, 8.850401597301765 18];
rcvr_6=[53.106565449092145, 8.850401597301765 18];
rcvr_7=[53.10656021590728, 8.850401597301765 18];
rcvr_8=[53.10655176231102, 8.850401597301765 18];

rcvr_9=[53.10657108482692, 8.850410314495912 18];
rcvr_10=[53.106565449092145, 8.850410314495912 18];
rcvr_11=[53.10656021590728, 8.850410314495912 18];
rcvr_12=[53.10655176231102, 8.850410314495912 18];

rcvr_13=[53.10657108482692, 8.850417690583926 18];
rcvr_14=[53.106565449092145, 8.850417690583926 18];
rcvr_15=[53.10656021590728, 8.850417690583926 18];
rcvr_16=[53.10655176231102, 8.850417690583926 18];

% Receiver Antenna locations in ECI
rcvr_pos_eci = [lla2eci(rcvr_1,utc);
    lla2eci(rcvr_2,utc);
    lla2eci(rcvr_3,utc);
    lla2eci(rcvr_4,utc);
    lla2eci(rcvr_5,utc);
    lla2eci(rcvr_6,utc);
    lla2eci(rcvr_7,utc);
    lla2eci(rcvr_8,utc);
    lla2eci(rcvr_9,utc);
    lla2eci(rcvr_10,utc);
    lla2eci(rcvr_11,utc);
    lla2eci(rcvr_12,utc);
    lla2eci(rcvr_13,utc);
    lla2eci(rcvr_14,utc);
    lla2eci(rcvr_15,utc);
    lla2eci(rcvr_16,utc)];


loc1 = lla2eci(rcvr_1,utc);
loc2 = [loc1(1,1)+0,loc1(1,2)+5,loc1(1,3)];
loc3 = [loc1(1,1)+0,loc1(1,2)+10,loc1(1,3)];
loc4 = [loc1(1,1)+0,loc1(1,2)+15,loc1(1,3)];
loc5 = [loc1(1,1)+5,loc1(1,2)+0,loc1(1,3)];
loc6 = [loc1(1,1)+5,loc1(1,2)+5,loc1(1,3)];
loc7 = [loc1(1,1)+5,loc1(1,2)+10,loc1(1,3)];
loc8 = [loc1(1,1)+5,loc1(1,2)+15,loc1(1,3)];
loc9 = [loc1(1,1)+10,loc1(1,2)+00,loc1(1,3)];
loc10 = [loc1(1,1)+10,loc1(1,2)+5,loc1(1,3)];
loc11 = [loc1(1,1)+10,loc1(1,2)+10,loc1(1,3)];
loc12 = [loc1(1,1)+10,loc1(1,2)+15,loc1(1,3)];
loc13 = [loc1(1,1)+15,loc1(1,2)+00,loc1(1,3)];
loc14 = [loc1(1,1)+15,loc1(1,2)+5,loc1(1,3)];
loc15 = [loc1(1,1)+15,loc1(1,2)+10,loc1(1,3)];
loc16 = [loc1(1,1)+15,loc1(1,2)+15,loc1(1,3)];

rcvr_locs = [loc1;loc2;loc3;loc4;
    loc5;loc6;loc7;loc8;
    loc9;loc10;loc11;loc12;
    loc13;loc14;loc15;loc16];

% relative movement with respect to the leader satellite
% ideal scenario, without external perturbations

figure();
hold on;
plot3(rcvr_locs(:,1),rcvr_locs(:,2),rcvr_locs(:,3),'r^','LineWidth',2,'MarkerSize',10);

title('ECI coordinates of receiver antennas on Earth')
rotate3d
xlabel('LVLH-x (metres)'), ylabel('LVLH-y (metres)'), zlabel('LVLH-z (metres)')
legend('4X4 MIMO ECI coordinates of Rx')
grid on;
axis equal;







