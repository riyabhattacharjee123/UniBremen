
loc1 = lla2eci(rcvr_1,utc);
ak=1;
bk=1;
rcvr_pos_eci = zeros(Nrx*Nrx,3)

for b = 1:Nrx
    for a = 1:Nrx
        rcvr_pos_eci(ak,:,:) = [(loc1(1,1)+((a-1)*gap)),... % x cood ECI
                (loc1(1,2)+((b-1)*gap)),... % y cood ECI
                (loc1(1,3)) ]; % z cood ECI
        ak=ak+1;        
    end    
end
figure();
hold on;
plot3(rcvr_pos_eci(:,1),rcvr_pos_eci(:,2),rcvr_pos_eci(:,3),'r^','LineWidth',2,'MarkerSize',10);

title('ECI coordinates of receiver antennas on Earth')
rotate3d
xlabel('ECI-x (metres)'), ylabel('ECI-y (metres)'), zlabel('ECI-z (metres)')
legend('MIMO ECI coordinates')
grid on;
axis equal;

