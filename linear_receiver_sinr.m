%%%%% linear_receiver_sinr.m %%%%
% This program calculates the SINR and Sum_Rate of the GS receiver antennas %

%% Calculate the SINR of each receiver antenna at GS %%
Sum_H_dl_i_1=zeros();
%for ut_k = 1:Ngs
 %   mult3 = W_lin(:,ut_k) * Channel_matrix_H(ut_k,:) ;
  %  Numerator1(ut_k)= det(mult3)^2 ;
    
   % for ut_i = 1:Ngs
    %    if ut_i ~= ut_k
     %       G_dl_i_1 = G_dl(ut_i,:);
      %      mult4=W_k(:,ut_k) * H_dl_i_1 ;
       %     element = det(mult4)^2;
        %    Sum_H_dl_i_1 = Sum_H_dl_i_1 + element;             
        %end       
    %end
    %W_k_norm = norm(transpose(conj(W_k(:,ut_k))));
    %Denominator1(ut_k)= Sum_H_dl_i_1+rho_ul_inv*W_k_norm^2;    
    %SINR(ut_k)=Numerator1(ut_k)/Denominator1(ut_k) ; %Tau
%end


for ut_k = 1:Ns
    
    H_ut_k = Channel_matrix_H(:,ut_k);
    H_ut_k_H = conj(H_ut_k.');
    
    mult3 = W_lin(ut_k,:) * H_ut_k;
    Numerator1(ut_k)= (det(mult3))^2;
    
    for ut_i = 1:Ns
        if ut_i ~= ut_k
            mult4 = W_lin(ut_k,:) * Channel_matrix_H(:,ut_i);
            element = (det(mult4))^2 ;
            Sum_H_dl_i_1 = Sum_H_dl_i_1 + element;
        end        
    end
    
    Denominator1 = Sum_H_dl_i_1 + noiseVar ...
                        * (W_lin(ut_k,:))' ...
                        * (W_lin(ut_k,:)) ;
                    
    SINR(ut_k) = Numerator1(ut_k)/Denominator1(ut_k) ; %Tau    
end

R_sum_rate = zeros();
MSE = zeros();
for ut_k = 1:Ns 
   %% Calculate Mean Squared Error %%
   MSE(ut_k) = inv(1+SINR(ut_k));   
   
   %% Calculate the Sum Rate of each antenna at GS %%
   R_sum_rate = R_sum_rate+(log2(1+SINR(ut_k)));   
end



