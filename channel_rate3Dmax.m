%%%% channel_rate3Dmax.m %%%%
%% calculate hermitian of steering matrix A downlink %%
steering_matrix_A_hermitian = transpose(conj(steering_matrix_A));
prod_A = steering_matrix_A_hermitian*steering_matrix_A ;
scaled_I_matrix_check = prod_A/Nr;
I_Nt = eye(Nt,Nt); % identity matrix
SNR_dB = 160; 
SNR = 10^(SNR_dB/10); % converting SNR from dB to linear scale
R_opt = log2(real(det(I_Nt + (SNR/Nt)*prod_A ))) ;
