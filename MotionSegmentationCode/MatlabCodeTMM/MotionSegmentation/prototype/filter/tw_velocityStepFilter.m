function [ nnidx_filtered, nndists_filtered ] = tw_velocityStepFilter( fmat, radius, nnidx, nndists )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    k_nn = size(nnidx, 1);
    frame_count = size(nnidx, 2);
    velocities = zeros(5, frame_count);
    for i=1:5
        offset = 15 + (i-1) * 3;
        vx = diff([fmat(1 + offset, :) fmat(1 + offset, end)]);
        vy = diff([fmat(2 + offset, :) fmat(2 + offset, end)]);
        vz = diff([fmat(3 + offset, :) fmat(3 + offset, end)]);
        velocities(i, :) = sqrt( vx.^2 + vy.^2 + vz.^2 );
    end

    winSize = 20;
    coefs = gausswin(winSize);
    coefs = coefs / sum(coefs);
    for i=1:5
        velocities(i, :) = conv(velocities(i, :), coefs, 'same');
    end

    threshold = 0.5;
    bodyPartsActive = sum(velocities > threshold);
    
    radius_filtered = [0.4, 0.4, 0.6, 0.8, 1.0, 1.0];
    radius_filtered = radius_filtered * radius;

    nnidx_filtered = nan(k_nn, frame_count);
    nndists_filtered = inf(k_nn, frame_count);
    
    for i=1:frame_count
        r = radius_filtered( bodyPartsActive(i) + 1 );  
        idx = nndists(:, i) <= r;
        nnidx_filtered(idx, i) = nnidx(idx, i);
        nndists_filtered(idx, i) = nndists(idx, i);
    end
end