%%%%%%%%%%%%%% project 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate the realizations
A=4;
N=700; %number of sample in each realization
random_rel=randi([0,1],500,101);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store every relization matrix
n = 500; % Number of matrices
m = 7; % Size of each matrix (number of rows)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
undelayed_data=zeros(500,707);
for i=1:500
    undelayed_data(i,:)=repelem(random_rel(i,:),1,7);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%construct all line codes
random_rel_polar=2*A*undelayed_data-A;
random_rel_polar_rz=zeros(500,707);

for i=1:500
    h1=1;
    for j=1:707
        if h1<=4
            h1=h1+1;
            random_rel_polar_rz(i,j)=random_rel_polar(i,j);
        else    
            random_rel_polar_rz(i,j)=0;   
        if rem(j,7)==0            
            h1=1;
        end    
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%take our bits only  
pol_RZ_delayed=zeros(500,700);
delay = randi([0 6],500,1);
for i=1:500
    
    pol_RZ_delayed(i,:) = random_rel_polar_rz (i, 8-delay(i) : 707-delay(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot first realization
time = 1:100;  % Define time variable from 1 to 100
figure(1)
stairs(time,pol_RZ_delayed(1,1:100)); % Plot the data using a stair-step plot (stairs) for the first 100 points of data_t 
axis([1 100 -8 8]) % Set the axis limits for the plot
title('realization for Polar RZ');
ylabel('magnitude');
xlabel('time');
grid on
%%
% Statistical mean
statistical_mean=zeros(1,700);
for i=1:500
    statistical_mean=statistical_mean+pol_RZ_delayed(i,:);
end
statistical_mean=statistical_mean/500;
figure(2);
plot(time,statistical_mean(1,1:100)); %%plot the stastistical mean
axis([1 100 -8 8]) % Set the axis limits for the plot
title('statistical mean for polar RZ');
ylabel('magnitude');
xlabel('time');
grid on;

% time mean
average_mean=zeros(1,500);
pol_rz_trans=pol_RZ_delayed';
for i=1:700
    average_mean=average_mean+pol_rz_trans(i,:);
end
average_mean=average_mean/700;
figure(3);
plot(time,average_mean(1,1:100)); %plots the time average mean over the range of time
axis([1 100 -8 8]) % Set the axis limits for the plot
title('Average time mean for polar RZ');
ylabel('magnitude');
xlabel('time');
grid on;
%%
% statistical autocorrelation
statistical_autocorrelation_matrix= zeros(700,700);%define matrix 700*700
for m = 1:700
for n = m:700
%each column contains the mean of multiplying column m by each column n
statistical_autocorrelation_matrix(m,n-m+1) = mean(pol_RZ_delayed(:,m).*pol_RZ_delayed(:,n));
end
end
statistical_autocorrelation_matrix = [fliplr(mean(statistical_autocorrelation_matrix)) mean(statistical_autocorrelation_matrix)];%flip the matrix and concatenate , as ACF is even function
statistical_autocorrelation_matrix = cat(2,statistical_autocorrelation_matrix(1,1:699),statistical_autocorrelation_matrix(1,701:1400));%remove statistical_autocorrelation_matrix(0) as it repeated twice
time = -699:699;
figure(4);
plot(time*10,statistical_autocorrelation_matrix);
axis([-200 200 -3 20])  % Set the axis limits for the plot
title('statistical autocorrelation for polar RZ');
ylabel('magnitude');
xlabel('tau');
grid on;

% time autocorrelation
realization_before = zeros(500,700);
new_vect = pol_RZ_delayed(1,:);%take the first realization(waveform)
for i=1:500
%we add random bits by getting last bits of realiztions according to value of shifting and put them
%in the beginning of realiztion
realization_before(i,:) = horzcat(new_vect(700-i:700),new_vect(1:700-i-1));
end
time_autocorrelation_matrix_T = zeros(1,501);
for m = 1:500
%each column contains the mean of multiplying first realiztion by
%shifted version of the first realiztion(tau from 1 to 500)
%we take the transpose as mean works on vector columns
time_autocorrelation_matrix_T(1,m) = mean(transpose(realization_before(1,:).*realization_before(m,:)));
end
time_autocorrelation_matrix_T = [fliplr(time_autocorrelation_matrix_T) time_autocorrelation_matrix_T];%flip the matrix and concatenate as ACF is even function
time_autocorrelation_matrix_T = cat(2,time_autocorrelation_matrix_T(1,1:501),time_autocorrelation_matrix_T(1,503:1002));%remove time_autocorrelation_matrix(0) as it repeated twice
time = -500:500;
figure(5);
plot(time*10,time_autocorrelation_matrix_T);%we multiply the time with the sample duration
title('time autocorrelationfor polar RZ');
ylabel('magnitude');
xlabel('tau');
axis([-500 500 -20 20])  % Set the axis limits for the plot
grid on;

% calculate bandwidth (use fourier transform)
ft_Rx = fft(statistical_autocorrelation_matrix.');
shift_Rx = fftshift(ft_Rx);
fs = 1399;%lenghth of statistical autocorrelation function
f = fs/2*linspace(-1,1,fs);
figure(6);
plot(f*100/fs,abs(shift_Rx));
title('PSD for polar RZ');
ylabel('magnitude');
xlabel('frequency');
axis([-100 100 0 40]) % Set the axis limits for the plot
grid on;
