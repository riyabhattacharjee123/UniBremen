
% Moon acceleration graph
for s=1:(length(r_start)/3)
    
    figure();
    
    % ECI-x Acceleration
    nexttile
    acc_x_1 = moon_acc(:,(3*s-2));
    ts_acc_x_1 = timeseries(acc_x_1,1:sim_points);
    ts_acc_x_1.Time = ts_acc_x_1.Time - ts_acc_x_1.Time(1); 
    plot(ts_acc_x_1)   
    
    title('Moon gravitation Perturbed Orbit Acceleration: ECI-x: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-x (m/s^2)')
    legend('Perturbed path'); 
    
    % ECI-y Acceleration
    nexttile
    acc_y_1 = moon_acc(:,(3*s-1));
    ts_acc_y_1 = timeseries(acc_y_1,1:sim_points);
    ts_acc_y_1.Time = ts_acc_y_1.Time - ts_acc_y_1.Time(1); 
    plot(ts_acc_y_1)   
    
    title('Moon gravitation Perturbed Orbit Acceleration: ECI-y: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-y (m/s^2)')
    legend('Perturbed path');
    
     % ECI-z Acceleration
    nexttile
    acc_z_1 = moon_acc(:,(3*s-0));
    ts_acc_z_1 = timeseries(acc_z_1,1:sim_points);
    ts_acc_z_1.Time = ts_acc_z_1.Time - ts_acc_z_1.Time(1); 
    plot(ts_acc_z_1)
       
    title('Moon gravitation Perturbed Orbit Acceleration: ECI-z: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-z (m/s^2)')
    legend('Perturbed path');   
    
end

% Sun acceleration graph
for s=1:(length(r_start)/3)
    
    figure();
    
    % ECI-x Acceleration
    nexttile
    acc_x_1 = sun_acc(:,(3*s-2));
    ts_acc_x_1 = timeseries(acc_x_1,1:sim_points);
    ts_acc_x_1.Time = ts_acc_x_1.Time - ts_acc_x_1.Time(1); 
    plot(ts_acc_x_1)   
    
    title('Sun gravitation Perturbed Orbit Acceleration: ECI-x: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-x (m/s^2)')
    legend('Perturbed path'); 
    
    % ECI-y Acceleration
    nexttile
    acc_y_1 = sun_acc(:,(3*s-1));
    ts_acc_y_1 = timeseries(acc_y_1,1:sim_points);
    ts_acc_y_1.Time = ts_acc_y_1.Time - ts_acc_y_1.Time(1); 
    plot(ts_acc_y_1)   
    
    title('Sun gravitation Perturbed Orbit Acceleration: ECI-y: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-y (m/s^2)')
    legend('Perturbed path');
    
     % ECI-z Acceleration
    nexttile
    acc_z_1 = sun_acc(:,(3*s-0));
    ts_acc_z_1 = timeseries(acc_z_1,1:sim_points);
    ts_acc_z_1.Time = ts_acc_z_1.Time - ts_acc_z_1.Time(1); 
    plot(ts_acc_z_1)
       
    title('Sun gravitation Perturbed Orbit Acceleration: ECI-z: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-z (m/s^2)')
    legend('Perturbed path');   
    
end

% SRP acceleration graph
for s=1:(length(r_start)/3)
    
    figure();
    
    % ECI-x Acceleration
    nexttile
    acc_x_1 = srp_acc(:,(3*s-2));
    ts_acc_x_1 = timeseries(acc_x_1,1:sim_points);
    ts_acc_x_1.Time = ts_acc_x_1.Time - ts_acc_x_1.Time(1); 
    plot(ts_acc_x_1)   
    
    title('SRP Perturbed Orbit Acceleration: ECI-x: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-x (m/s^2)')
    legend('Perturbed path'); 
    
    % ECI-y Acceleration
    nexttile
    acc_y_1 = srp_acc(:,(3*s-1));
    ts_acc_y_1 = timeseries(acc_y_1,1:sim_points);
    ts_acc_y_1.Time = ts_acc_y_1.Time - ts_acc_y_1.Time(1); 
    plot(ts_acc_y_1)   
    
    title('SRP Perturbed Orbit Acceleration: ECI-y: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-y (m/s^2)')
    legend('Perturbed path');
    
     % ECI-z Acceleration
    nexttile
    acc_z_1 = srp_acc(:,(3*s-0));
    ts_acc_z_1 = timeseries(acc_z_1,1:sim_points);
    ts_acc_z_1.Time = ts_acc_z_1.Time - ts_acc_z_1.Time(1); 
    plot(ts_acc_z_1)
       
    title('SRP Perturbed Orbit Acceleration: ECI-z: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-z (m/s^2)')
    legend('Perturbed path');   
    
end

% J2 Effects acceleration graph
for s=1:(length(r_start)/3)
    
    figure();
    
    % ECI-x Acceleration
    nexttile
    acc_x_1 = J2_acc(:,(3*s-2));
    ts_acc_x_1 = timeseries(acc_x_1,1:sim_points);
    ts_acc_x_1.Time = ts_acc_x_1.Time - ts_acc_x_1.Time(1); 
    plot(ts_acc_x_1)   
    
    title('J2 Effects Perturbed Orbit Acceleration: ECI-x: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-x (m/s^2)')
    legend('Perturbed path'); 
    
    % ECI-y Acceleration
    nexttile
    acc_y_1 = J2_acc(:,(3*s-1));
    ts_acc_y_1 = timeseries(acc_y_1,1:sim_points);
    ts_acc_y_1.Time = ts_acc_y_1.Time - ts_acc_y_1.Time(1); 
    plot(ts_acc_y_1)   
    
    title('J2 Effects Perturbed Orbit Acceleration: ECI-y: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-y (m/s^2)')
    legend('Perturbed path');
    
     % ECI-z Acceleration
    nexttile
    acc_z_1 = J2_acc(:,(3*s-0));
    ts_acc_z_1 = timeseries(acc_z_1,1:sim_points);
    ts_acc_z_1.Time = ts_acc_z_1.Time - ts_acc_z_1.Time(1); 
    plot(ts_acc_z_1)
       
    title('J2 Effects Perturbed Orbit Acceleration: ECI-z: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-z (m/s^2)')
    legend('Perturbed path');   
    
end

% Air Drag Effects acceleration graph
for s=1:(length(r_start)/3)
    
    figure();
    
    % ECI-x Acceleration
    nexttile
    acc_x_1 = airdrag_acc(:,(3*s-2));
    ts_acc_x_1 = timeseries(acc_x_1,1:sim_points);
    ts_acc_x_1.Time = ts_acc_x_1.Time - ts_acc_x_1.Time(1); 
    plot(ts_acc_x_1)   
    
    title('Air Drag Perturbed Orbit Acceleration: ECI-x: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-x (m/s^2)')
    legend('Perturbed path'); 
    
    % ECI-y Acceleration
    nexttile
    acc_y_1 = airdrag_acc(:,(3*s-1));
    ts_acc_y_1 = timeseries(acc_y_1,1:sim_points);
    ts_acc_y_1.Time = ts_acc_y_1.Time - ts_acc_y_1.Time(1); 
    plot(ts_acc_y_1)   
    
    title('Air Drag Perturbed Orbit Acceleration: ECI-y: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-y (m/s^2)')
    legend('Perturbed path');
    
     % ECI-z Acceleration
    nexttile
    acc_z_1 = airdrag_acc(:,(3*s-0));
    ts_acc_z_1 = timeseries(acc_z_1,1:sim_points);
    ts_acc_z_1.Time = ts_acc_z_1.Time - ts_acc_z_1.Time(1); 
    plot(ts_acc_z_1)
       
    title('Air Drag Perturbed Orbit Acceleration: ECI-z: satellite- ',num2str(s))
    xlabel('Time(seconds)'), ylabel('ECI-z (m/s^2)')
    legend('Perturbed path');   
    
end
