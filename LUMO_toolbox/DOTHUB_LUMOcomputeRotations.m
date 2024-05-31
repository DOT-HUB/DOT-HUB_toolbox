function angle_deg = DOTHUB_LUMOcomputeRotations(xy_original, xy_rotated, plotFlag)
% This function estimates the rotation of source A from the loaded 2D coordinates
% after comparing the 2D coordinates of the original (fixed position) and rotated (loaded) tiles. 
% The function does the following:
%   1. Find the centroid for the original and rotated sources (A,B,C),
%      treating them as triangles.
%   2. Calculate the vector length from the centroid to the reference
%      vertex (source A) for both triangles.
%   3. Compute the cosine of the angle between the two vectors, return the
%      angle in radians and convert it into degrees. 
%
% ######################## INPUTS #########################################
%
% xy_original        :   The 2D coordinates of the original ABC sources,
%                        making a triangle.
% xy_rotaded         :   The 2D coordinates of the rotated (loaded) ABC 
%                        sources, making a triangle.
% plotFlag           :   A boolean value to enable a visualisation of the
%                        original and rotated vectors
%
% ####################### OUTPUTS #########################################
% angle_deg          :   The estimated degree of rotation
%
% ####################### Dependencies ####################################
%
% #########################################################################
% GL, UCL, April 2024

% #########################################################################    

    % extract x, y original and rotated
    x_original = xy_original(:, 1);
    y_original = xy_original(:, 2);
    x_rotated = xy_rotated(:, 1);
    y_rotated = xy_rotated(:, 2);

    % Find centroid for original and rotated coordinates
    xor_centroid = mean(x_original);
    yor_centroid = mean(y_original);
    xrot_centroid = mean(x_rotated);
    yrot_centroid = mean(y_rotated);

    % Calculate vector components from centroid to source 1 (reference source)
    vx_original = x_original(1) - xrot_centroid;
    vy_original = y_original(1) - yrot_centroid;
    vx_rotated = x_rotated(1) - xrot_centroid;
    vy_rotated = y_rotated(1) - yrot_centroid;

    %% Plot
    if plotFlag
        figure(1)
        scatter(x_original, y_original, 'r');
        hold on
        scatter(xor_centroid, yor_centroid, 'black');
        % Plot the vector from centroid to source 1
        quiver(xor_centroid, yor_centroid, vx_original, vy_original, 0, 'red', 'LineWidth', 2);
        hold on
        % Overlay the rotated triangle
        scatter(x_rotated, y_rotated, 'blue');
        hold on
        scatter(xrot_centroid, yrot_centroid, 'blue');
        % Plot the vector from centroid to rotated source
        quiver(xrot_centroid, yrot_centroid, vx_rotated, vy_rotated, 0, 'blue', 'LineWidth', 2);
        hold off
    end

    %%
    % Compute the dot product of vectors v1 and v2
    dot_product = dot([vx_original, vy_original], [vx_rotated, vy_rotated]);
    
    % Establish if rotation is clockwise or anti-clockwise
    if vx_rotated < vx_original
       % counter-clockwise rotation
       rotation_direction = 'left-side';
    elseif vx_rotated > vx_original
        % clockwise rotation
        rotation_direction = 'right-side';
    else
        rotation_direction = '180';
    end 
    
    % Compute the magnitudes of vectors v1 and v2
    magnitude_v1 = norm([vx_original, vy_original]);
    magnitude_v2 = norm([vx_rotated, vy_rotated]);

    % Compute the cosine of the angle between v1 and v2
    cos_angle = dot_product / (magnitude_v1 * magnitude_v2);

    % Compute the angle in radians
    angle_rad = acos(cos_angle);

    % Convert the angle from radians to degrees
    angle_deg = rad2deg(angle_rad);
    
    % Adjust for rotation on the right-side of original vector
    if strcmp(rotation_direction, 'right-side')
        angle_deg = 360 - angle_deg;
    end 