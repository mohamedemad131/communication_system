data=randi([0,1],500,101); %create the ensemble data    
A=4; % amplitude to represent one  
data=A*data; % map 1 to 4 and zero to 0
data_1=zeros(500,707);% define the sampling matrix
data_2=zeros(500,707);% define the final data matrix with delay
for i=1:500 
 data_1(i,:)=repelem(data(i,:),1,7); %sampling each element by seven samples
 delay_value=randi([0,6],1,1); % generate random delay
 data_2(i,:)=circshift( data_1(i,:),delay_value); % shift right the matrix to add the delay at the start of each realization
end
data_t=data_2(:, 1:700); %we removed the extra bits where added to generate random delay
data_t=data_t.'; %to make each coloum represent a realization

% plot first realization
time = 1:100;  % Define time variable from 1 to 100
figure(1)
stairs(time,data_t(1:100)); % Plot the data using a stair-step plot (stairs) for the first 100 points of data_t 
axis([1 100 -8 8]) % Set the axis limits for the plot
title('realization');
ylabel('magnitude');
xlabel('time');
grid on

% Statistical mean

data_t = data_t.'; %transpose the matrix to make every row is a realization
statistical_mean=zeros(1,700);
for i=1:500
    statistical_mean=statistical_mean+data_t(i,:);
end
statistical_mean=statistical_mean/500;
figure(2);
plot(time,statistical_mean(1,1:100)); %%plot the stastistical mean
axis([1 100 -8 8]) % Set the axis limits for the plot
title('statistical mean');
ylabel('magnitude');
xlabel('time');
grid on;

% time mean
average_mean=zeros(1,500);
pol_rz_trans=data_t';
for i=1:700
    average_mean=average_mean+pol_rz_trans(i,:);
end
average_mean=average_mean/700;
figure(3);
plot(time,average_mean(1,1:100)); %plots the time average mean over the range of time
axis([1 100 -8 8]) % Set the axis limits for the plot
title('Average time mean');
ylabel('magnitude');
xlabel('time');
grid on;

% statistical autocorrelation
statistical_autocorrelation_matrix= zeros(700,700);%define matrix 700*700
for m = 1:700
for n = m:700
%each column contains the mean of multiplying column m by each column n
statistical_autocorrelation_matrix(m,n-m+1) = mean(data_t(:,m).*data_t(:,n));
end
end
statistical_autocorrelation_matrix = [fliplr(mean(statistical_autocorrelation_matrix)) mean(statistical_autocorrelation_matrix)];%flip the matrix and concatenate , as ACF is even function
statistical_autocorrelation_matrix = cat(2,statistical_autocorrelation_matrix(1,1:699),statistical_autocorrelation_matrix(1,701:1400));%remove statistical_autocorrelation_matrix(0) as it repeated twice
time = -699:699;
figure(4);
plot(time*10,statistical_autocorrelation_matrix);
axis([-200 200 -3 20])  % Set the axis limits for the plot
title('statistical autocorrelation');
ylabel('magnitude');
xlabel('tau');
grid on;

% time autocorrelation
realization_before = zeros(500,700);
new_vect = data_t(1,:);%take the first realization(waveform)
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
title('time autocorrelation');
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
title('PSD');
ylabel('magnitude');
xlabel('frequency');
axis([-100 100 0 40]) % Set the axis limits for the plot
grid on;
