function [com1, vect0] = detectFirstCOM(path_to_video)

%% Specify manually the two points allowing tail tracking
    
    % Plot first frame of video
    fig = figure;
    vid = VideoReader(path_to_video);
    im = readFrame(vid);
    im = mean(im, 3);
    image(im, 'CDataMapping', 'scaled')
    axis equal
    
    % Get the two points from user
    fprintf('Select first point on fish tail \n');
    [x_1, y_1] = getpts(fig);
    fprintf('Select second point on fish tail \n');
    [x_2, y_2] = getpts(fig);
    
    % Compute first center of mass and initial vector
    com1 = [y_1, x_1];
    vect0 = [y_2-y_1, x_2-x_1];
    vect0 = vect0 ./ sqrt(sum(vect0.^2));
    close(fig);
end
