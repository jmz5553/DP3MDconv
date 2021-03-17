function [avg, std] = zAngles(direction)
%This program calculates each of the angles between a set of direction vectors and
%the z axis. It returns the average angle and standard deviation of those
%angles

[x, N] = size(direction);

if x ~= 3
    error('Please enter a matrix of N 3-dimensional direction vectors');
end

z = [0 0 1]';

angles = zeros(1, N);

% Calculate the angle of each vector with the vertical
for i = 1:N
    dDotZ = dot(direction(:, i), z);
    angles(i) = acos(dDotZ/norm(direction(:, i)));
end

% Calculate average angle
avg = sum(angles)/N;

% Calculate standard deviation
std = sqrt(sum((angles - avg).^2)/N);

h =  findobj('type','figure');
n = length(h);

figure(3 + n)
histogram(angles, ceil(sqrt(N)), 'BinLimits', [0, pi], 'Normalization', 'probability');
end

