clear all;
close all;

% Load the signals
load('x_plus_v1.mat'); % Load the signal x_plus_v1
load('v2.mat'); % Load the signal v2

% Extract the signals
x_plus_v1 = x_plus_v1; % Signal containing music and noise
v2 = v2; % Signal containing noise picked up by the external microphone

% Sampling rate
fs = 14700;

% Play the original mixed signal
%sound(x_plus_v1, fs);
%pause(10); % Pause to allow the sound to play

% cross-correlation to find the delay
%[corr, lags] = xcorr(x_plus_v1, v2);

% Find the maximum correlation and the corresponding delay
%[~, I] = max(abs(corr));
%max_delay = abs(lags(I));

% determine the filter length
%filter_length = max_delay;


% RLS Adaptive Filter parameters
filter_length = round(0.01 * fs); % 10 ms filter length
lambda = 1; % Forgetting factor
delta = 0.01; % Initialization parameter

% Initialize filter coefficients, inverse correlation matrix, and error signal
w = zeros(filter_length, 1);
P = (1/delta) * eye(filter_length);
e = zeros(length(x_plus_v1), 1);

% Implement RLS Adaptive Filter
for n = filter_length:length(x_plus_v1)
    x_n = v2(n:-1:n-filter_length+1); % Input vector
    k = (P * x_n) / (lambda + x_n' * P * x_n); % Gain vector
    alpha = x_plus_v1(n) - w' * x_n; % A priori error
    w = w + k * alpha; % Update filter coefficients
    P = (P - k * x_n' * P) / lambda; % Update inverse correlation matrix
    e(n) = alpha; % Error signal
end

% Play the filtered signal
sound(e, fs);

% Plot the results
t = (0:length(x_plus_v1)-1) / fs;
figure;
plot(t, x_plus_v1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Mixed Signal');

figure;
plot(t, v2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Noise Signal (v_2(n))');

figure;
plot(t, e);
xlabel('Time (s)');
ylabel('Amplitude');
ylim([-0.8 0.8]);
title('Filtered Signal with RLS');

% Save the filtered signal
audiowrite('filtered_signal_rls.wav', e, fs);