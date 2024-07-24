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

% LMS Adaptive Filter parameters
filter_length = round(0.01 * fs); % 10 ms filter length
mu = 0.01; % Step size for the LMS algorithm

% Initialize filter coefficients and error signal
w = zeros(filter_length, 1);
e = zeros(length(x_plus_v1), 1);

% Implement LMS Adaptive Filter
for n = filter_length:length(x_plus_v1)
    x_n = v2(n:-1:n-filter_length+1); % Input vector
    y = w' * x_n; % Filter output
    e(n) = x_plus_v1(n) - y; % Error signal
    w = w + mu * e(n) * x_n; % Update filter coefficients
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
title('Noise Signal v_2(n)');

figure;
plot(t, e);
xlabel('Time (s)');
ylabel('Amplitude');
ylim([-0.8 0.8]);
title('Filtered Signal with LMS');

% Save the filtered signal
audiowrite('filtered_signal.wav', e, fs);