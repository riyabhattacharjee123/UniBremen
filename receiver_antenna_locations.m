%%%% receiver_antenna_locations.m %%%%
%%% gives the lcations of receiving antennas arranged in Unified Planar
%%% Array (UPA) .
%%% All the coordinates are in ECI

ant_location = lla2eci(rcvr_1,utc);
ak=1;
bk=1;
rcvr_pos_eci = zeros(Nrx*Nrx,3)
Nr = Nrx*Nrx; % number of receiving antennas at GS
Ngs = Nr ; % number of antennas at GS
Ngsx = Nrx; % number of antennas at x_axis in GS

for b = 1:Nrx
    for a = 1:Nrx
        rcvr_pos_eci(ak,:,:) = [(ant_location(1,1)+((a-1)*D_Ant)),... % x cood ECI
                (ant_location(1,2)+((b-1)*D_Ant)),... % y cood ECI
                (ant_location(1,3)) ]; % z cood ECI
        ak=ak+1;        
    end    
end


%figure();
%hold on;
%plot3(rcvr_pos_eci(:,1),rcvr_pos_eci(:,2),rcvr_pos_eci(:,3),'r^'...
 %   ,'LineWidth',2,'MarkerSize',10);

%title('ECI coordinates of receiver antennas on Earth')
%rotate3d
%xlabel('ECI-x (metres)'), ylabel('ECI-y (metres)'),...
 %   zlabel('ECI-z (metres)')
%legend('MIMO ECI coordinates')
%grid on;
%axis equal;

