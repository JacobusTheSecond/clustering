function [] = tw_velocityPlot(fmat)
% id 4 name ltibia
% id 9 name rtibia
% id 17 name head
% id 21 name lwrist
% id 28 name rwrist
%[4,9,17,21,28]

frame_count = size(fmat,2);
frames = 1:frame_count;
velocities = zeros(5, frame_count);
for i=1:5
    offset = 15 + (i-1) * 3;
    vx = diff([fmat(1 + offset, :) fmat(1 + offset, end)]);
    vy = diff([fmat(2 + offset, :) fmat(2 + offset, end)]);
    vz = diff([fmat(3 + offset, :) fmat(3 + offset, end)]);
    velocities(i, :) = sqrt( vx.^2 + vy.^2 + vz.^2 );
end

vmax = max(max(velocities));

subplot(3,2,1);
plot(frames, velocities(1, :));
title('left tibia');
ylim([0 vmax]);

subplot(3,2,2);
plot(frames, velocities(2, :));
title('right tibia');
ylim([0 vmax]);

subplot(3,2,3);
plot(frames, velocities(4, :));
title('left wrist');
ylim([0 vmax]);

subplot(3,2,4);
plot(frames, velocities(5, :));
title('right wrist');
ylim([0 vmax]);

subplot(3,2,5);
plot(frames, velocities(3, :));
title('head');
ylim([0 vmax]);

%plot(frames, velocities(1, :), frames, velocities(2, :), frames, velocities(3, :), frames, velocities(4, :), frames, velocities(5, :));
%legend('ltibia', 'rtibia', 'head', 'lwrist', 'rwrist');
end