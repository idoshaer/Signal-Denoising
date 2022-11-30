%% -------------- Example 1 --------------
% In this example we will run a MATLAB function that do automatic
% denoising, we will only give the function a number of input parameters:
% 1. 
close all; clear all; clc;

% Loading audio file
[male, Fs1] = audioread('male.wav'); 
[female, Fs2] = audioread('female.wav');

% Transpose audio signals
male = male.';
female = female.';
male = male(1:20e3);
female = female(1:20e3);

SNR = 3; % SNR value in dB for adding awgn noise
[male_noisy,~] = addGaussianNoise(male,SNR);
[female_noisy,~] = addGaussianNoise(female,SNR);

MW = ["db8","db20","coif2","coif5","bior3.1","bior3.3"];

denoised1_vec = zeros(length(MW),length(female));
denoised2_vec = zeros(length(MW),length(female));
denoised3_vec = zeros(length(MW),length(female));
denoised4_vec = zeros(length(MW), length(female));

snr_values = zeros(4,length(MW));

% male, level = 3, threshold = SURE
figure('Name','Male, Level = 3, Threshold method: SURE');
for i = 1:(length(MW))
  
    [denoised1] = wdenoise(male_noisy,3,'Wavelet',MW(i),'NoiseEstimate','LevelDependent','DenoisingMethod','SURE','ThresholdRule','Soft');
    denoised1_vec(i,:) = denoised1;
    [snr_values(1,i),~,~,~] = Eval_SNR_Corr(denoised1,male);
    
    subplot(2,3,i);
    plot(denoised1_vec(i,:));
    grid on; ylabel('Amplitude'); xlabel('t');
    title({'Male, Mother Wavelet = ' + MW(i) ; ['SNR = ', num2str(snr_values(1,i)), 'dB']});
end

% female, level = 3, threshold = SURE
figure('Name','Female, Level = 3, Threshold method: SURE');
for i = 1:(length(MW))
  
    [denoised2] = wdenoise(female_noisy,3,'Wavelet',MW(i),'NoiseEstimate','LevelDependent','DenoisingMethod','SURE','ThresholdRule','Soft');
    denoised2_vec(i,:) = denoised2;
    [snr_values(2,i),~,~,~] = Eval_SNR_Corr(denoised2,female);
    
    subplot(2,3,i);
    plot(denoised2_vec(i,:));
    grid on; ylabel('Amplitude'); xlabel('t');
    title({'Female, Mother Wavelet = ' + MW(i) ; ['SNR = ', num2str(snr_values(2,i)), 'dB']});
end

% male, level = 4, threshold = SURE
figure('Name','Male, Level = 4, Threshold method: SURE');
for i = 1:(length(MW))
        
    [denoised3] = wdenoise(male_noisy,4,'Wavelet',MW(i),'NoiseEstimate','LevelDependent','DenoisingMethod','SURE','ThresholdRule','Soft');
    denoised3_vec(i,:) = denoised3;
    [snr_values(3,i),~,~,~] = Eval_SNR_Corr(denoised3,male);
    
    subplot(2,3,i);
    plot(denoised3_vec(i,:));
    grid on; ylabel('Amplitude'); xlabel('t');
    title({'Male, Mother Wavelet = ' + MW(i) ; ['SNR = ', num2str(snr_values(3,i)), 'dB']});
end

% female, level = 4, threshold = SURE
figure('Name','Female, Level = 4, Threshold method: SURE');
for i = 1:(length(MW))
        
    [denoised4] = wdenoise(female_noisy,4,'Wavelet',MW(i),'NoiseEstimate','LevelDependent','DenoisingMethod','SURE','ThresholdRule','Soft');
    denoised4_vec(i,:) = denoised4;
    [snr_values(4,i),~,~,~] = Eval_SNR_Corr(denoised4,female);
    
    subplot(2,3,i);
    plot(denoised1_vec(i,:));
    grid on; ylabel('Amplitude'); xlabel('t');
    title({'Female, Mother Wavelet = ' + MW(i) ; ['SNR = ', num2str(snr_values(4,i)), 'dB']});
end

figure(5);
subplot(2,1,1)
plot(male);
grid on; axis tight; ylabel('Amplitude'); xlabel('t');
title('Male Original Signal');

subplot(2,1,2)
plot(female);
grid on; axis tight; ylabel('Amplitude'); xlabel('t');
title('Female Original Signal');

figure(6)
subplot(2,1,1)
plot(male_noisy);
grid on; axis tight; ylabel('Amplitude'); xlabel('t');
title('Male Noisy Signal with SNR = 3dB');

subplot(2,1,2)
plot(female_noisy);
grid on; axis tight; ylabel('Amplitude'); xlabel('t');
title('Female Noisy Signal with SNR = 3dB');

figure(7);
for j = 1:4
    
    subplot(2,2,j);
    snr = snr_values(j,:);
    [maxVal,index] = max(snr);
    if j == 1
        plot(denoised1_vec(index,:))
        grid on; ylabel('Amplitude'); xlabel('t');
        title({'Male, Level = 3, Mother Wavelet = ' + MW(index) ; ['SNR = ', num2str(maxVal), 'dB']});
    elseif j == 2
        plot(denoised2_vec(index,:))
        grid on; ylabel('Amplitude'); xlabel('t');
        title({'Female, Level = 3, Mother Wavelet = ' + MW(index) ; ['SNR = ', num2str(maxVal), 'dB']});
    elseif j == 3
        plot(denoised3_vec(index,:))
        grid on; ylabel('Amplitude'); xlabel('t');
        title({'Male, Level = 4, Mother Wavelet = ' + MW(index) ; ['SNR = ', num2str(maxVal), 'dB']});
    elseif j == 4
        plot(denoised4_vec(index,:))
        grid on; ylabel('Amplitude'); xlabel('t');
        title({'Female, Level = 4, Mother Wavelet = ' + MW(index) ; ['SNR = ', num2str(maxVal), 'dB']});
    end
end

%% -------------- Example 2 --------------
% In this example we will run a MATLAB function that do automatic

%% Part 1, Male audio
close all; clear all; clc;

[male, Fs1] = audioread('male.wav');

% Transpose audio signals
male = male.';
male = male(1:20e3);

SNR = 5; % SNR value in dB for adding awgn noise
[male_noisy,~] = addGaussianNoise(male,SNR);

figure('Name','Noisy signal');
plot(male_noisy); axis tight; xlabel('t'); ylabel('Amplitude'); grid on;
title('Noisy signal with SNR = 5dB');

dwtmode('status','nodisplay');
wname = ["db8","db20","coif2","coif5","bior3.1","bior3.3"];
THR = ["sqtwolog","rigrsure"];
denoised_Sqtwolog = zeros(length(wname),length(male));
denoised_Rigrsure = zeros(length(wname),length(male));
denoised_Sqtwolog2 = zeros(length(wname),length(male));
denoised_Rigrsure2 = zeros(length(wname),length(male));
SNRVal = zeros(length(THR),length(wname));
SNRVal2 = zeros(length(THR),length(wname));

% Level = 3
for i = 1:length(THR)
    
    for j = 1:length(wname)
        if i == 1
            denoised_Sqtwolog(j,:) = wden(male_noisy,THR(i),'s','sln',3,wname(j));
            SNRVal(i,j) = -20*log10(norm(abs(male-denoised_Sqtwolog(j,:)))/norm(male));
        else
            denoised_Rigrsure(j,:) = wden(male_noisy,THR(i),'s','sln',3,wname(j));
            SNRVal(i,j) = -20*log10(norm(abs(male-denoised_Rigrsure(j,:)))/norm(male));
        end
    end
    
end

% Level = 4
for i = 1:length(THR)
    
    for j = 1:length(wname)
        if i == 1
            denoised_Sqtwolog2(j,:) = wden(male_noisy,THR(i),'s','sln',4,wname(j));
            SNRVal2(i,j) = -20*log10(norm(abs(male-denoised_Sqtwolog2(j,:)))/norm(male));
        else
            denoised_Rigrsure2(j,:) = wden(male_noisy,THR(i),'s','sln',4,wname(j));
            SNRVal2(i,j) = -20*log10(norm(abs(male-denoised_Rigrsure2(j,:)))/norm(male));
        end
    end
    
end

PlotAll(denoised_Sqtwolog,denoised_Sqtwolog2,denoised_Rigrsure,denoised_Rigrsure2,wname,SNRVal,SNRVal2);

% Find best SNR values and plot
figure('Name','Best SNR values for every Mother wavelet/THR/Level');

[maxSNR1,I1] = max(SNRVal(1,:)); % Sqtwolog, Level = 3
best_den1 = denoised_Sqtwolog(I1,:);
[maxSNR2,I2] = max(SNRVal(2,:)); % Rigrsure, Level = 3
best_den2 = denoised_Rigrsure(I2,:);
[maxSNR3,I3] = max(SNRVal2(1,:)); % Sqtwolog, Level = 4
best_den3 = denoised_Sqtwolog2(I3,:);
[maxSNR4,I4] = max(SNRVal2(2,:)); % Rigrsure, Level = 4
best_den4 = denoised_Rigrsure2(I4,:);

subplot(2,2,1);
plot(best_den1); axis tight; xlabel('t'); ylabel('Amplitude'); grid on;
title({['Level = 3, THR: Sqtwolog, SNR = ', num2str(maxSNR1)] ; 'Mother Wavelet: ' + wname(I1)});
subplot(2,2,2);
plot(best_den2); axis tight; xlabel('t'); ylabel('Amplitude'); grid on;
title({['Level = 3, THR: Rigrsure, SNR = ', num2str(maxSNR2)] ; 'Mother Wavelet: ' + wname(I2)});
subplot(2,2,3);
plot(best_den3); axis tight; xlabel('t'); ylabel('Amplitude'); grid on;
title({['Level = 4, THR: Sqtwolog, SNR = ', num2str(maxSNR3)] ; 'Mother Wavlet: ' + wname(I3)});
subplot(2,2,4);
plot(best_den4); axis tight; xlabel('t'); ylabel('Amplitude'); grid on;
title({['Level = 4, THR: Rigrsure, SNR = ', num2str(maxSNR4)] ; 'Mother Wavelet: ' + wname(I4)});

%% Part 2, Female audio
close all; clear all; clc;

[female, Fs1] = audioread('female.wav');

% Transpose audio signals
female = female.';
female = female(1:20e3);

SNR = 5; % SNR value in dB for adding awgn noise
[male_noisy,~] = addGaussianNoise(female,SNR);

figure('Name','Noisy signal');
plot(male_noisy); axis tight; xlabel('t'); ylabel('Amplitude'); grid on;
title('Noisy signal with SNR = 5dB');

dwtmode('status','nodisplay');
wname = ["db8","db20","coif2","coif5","bior3.1","bior3.3"];
THR = ["sqtwolog","rigrsure"];
denoised_Sqtwolog = zeros(length(wname),length(female));
denoised_Rigrsure = zeros(length(wname),length(female));
denoised_Sqtwolog2 = zeros(length(wname),length(female));
denoised_Rigrsure2 = zeros(length(wname),length(female));
SNRVal = zeros(length(THR),length(wname));
SNRVal2 = zeros(length(THR),length(wname));

% Level = 3
for i = 1:length(THR)
    
    for j = 1:length(wname)
        if i == 1
            denoised_Sqtwolog(j,:) = wden(male_noisy,THR(i),'s','sln',3,wname(j));
            SNRVal(i,j) = -20*log10(norm(abs(female-denoised_Sqtwolog(j,:)))/norm(female));
        else
            denoised_Rigrsure(j,:) = wden(male_noisy,THR(i),'s','sln',3,wname(j));
            SNRVal(i,j) = -20*log10(norm(abs(female-denoised_Rigrsure(j,:)))/norm(female));
        end
    end
    
end

% Level = 4
for i = 1:length(THR)
    
    for j = 1:length(wname)
        if i == 1
            denoised_Sqtwolog2(j,:) = wden(male_noisy,THR(i),'s','sln',4,wname(j));
            SNRVal2(i,j) = -20*log10(norm(abs(female-denoised_Sqtwolog2(j,:)))/norm(female));
        else
            denoised_Rigrsure2(j,:) = wden(male_noisy,THR(i),'s','sln',4,wname(j));
            SNRVal2(i,j) = -20*log10(norm(abs(female-denoised_Rigrsure2(j,:)))/norm(female));
        end
    end
    
end

PlotAll(denoised_Sqtwolog,denoised_Sqtwolog2,denoised_Rigrsure,denoised_Rigrsure2,wname,SNRVal,SNRVal2);

% Find best SNR values and plot
figure('Name','Best SNR values for every Mother wavelet/THR/Level');

[maxSNR1,I1] = max(SNRVal(1,:)); % Sqtwolog, Level = 3
best_den1 = denoised_Sqtwolog(I1,:);
[maxSNR2,I2] = max(SNRVal(2,:)); % Rigrsure, Level = 3
best_den2 = denoised_Rigrsure(I2,:);
[maxSNR3,I3] = max(SNRVal2(1,:)); % Sqtwolog, Level = 4
best_den3 = denoised_Sqtwolog2(I3,:);
[maxSNR4,I4] = max(SNRVal2(2,:)); % Rigrsure, Level = 4
best_den4 = denoised_Rigrsure2(I4,:);

subplot(2,2,1);
plot(best_den1); axis tight; xlabel('t'); ylabel('Amplitude'); grid on;
title({['Level = 3, THR: Sqtwolog, SNR = ', num2str(maxSNR1)] ; 'Mother Wavelet: ' + wname(I1)});
subplot(2,2,2);
plot(best_den2); axis tight; xlabel('t'); ylabel('Amplitude'); grid on;
title({['Level = 3, THR: Rigrsure, SNR = ', num2str(maxSNR2)] ; 'Mother Wavelet: ' + wname(I2)});
subplot(2,2,3);
plot(best_den3); axis tight; xlabel('t'); ylabel('Amplitude'); grid on;
title({['Level = 4, THR: Sqtwolog, SNR = ', num2str(maxSNR3)] ; 'Mother Wavlet: ' + wname(I3)});
subplot(2,2,4);
plot(best_den4); axis tight; xlabel('t'); ylabel('Amplitude'); grid on;
title({['Level = 4, THR: Rigrsure, SNR = ', num2str(maxSNR4)] ; 'Mother Wavelet: ' + wname(I4)});