%% Video Analysis Template
% Sharbat 2020
% Thanks to Hipppolyte 2019
clear; close all; clc;

%% Insert video folder
% path2video = '/run/user/1000/gvfs/sftp:host=134.157.132.157/home/ljp/Zebraid2/User/Natalia/Backup_Natalia_Coolcat_/Rolling/2017-12/2017-12-13/run.1/video.avi';
videofolder = '/run/user/1000/gvfs/sftp:host=134.157.132.157/home/ljp/Zebraid2/User/Natalia/Backup_Natalia_Coolcat_/Rolling/2017-12';
datesfolder = dir([videofolder,'/2017*']);

% path2video = '/home/ljp/Science/Projects/Sharbat/Data/2021-03-18/Run 02/Images/tailconcat.avi'


%% Import video and create videoreader, invert greyscale if needed

for k = 1:length(datesfolder)
    
    runs = dir([datesfolder(k).folder,'/',datesfolder(k).name,'/*']);
    runs = runs(~startsWith({runs.name}, '.'));
    for i = 1:length(runs)
        path2video = [runs(i).folder, '/', runs(i).name, '/video.avi'];
        try
            vid = VideoReader(path2video);
            numframes = floor(vid.Duration * vid.FrameRate);
            pause_frame = 0.05;

            frames = read(vid);
            clear corr;
    %         size(frames,3)

            for j=1:numframes
                if size(frames,3) == 3
                    corr(:,:,:,j)=imcomplement(rgb2gray(frames(:,:,:,j)));
                end

                if size(frames,3) == 1
                    corr(:,:,:,j)=imcomplement(frames(:,:,:,j));
                end

            end
            writepath = [runs(i).folder, '/', runs(i).name,'/video_inv.avi']
            vid_new = VideoWriter(writepath,'Grayscale AVI');

            open(vid_new);
            writeVideo(vid_new, corr);
            close(vid_new);
        catch
            warning('no video file');
        end
    end
end
%% Parameters for tail tracking

num_segments = 10;
inertia = 0.5;
body_length = 0; %don't need this, head not recorded
tail_length = 125; %change for different fish
initial_box = 1.0; %1.2 initially
box_increment = 0.001; %0.01 initially

%% Check interactive tracking of tail
trackingPlot_manual(path2video, 'num_segments', num_segments, 'inertia', inertia, ...
                    'body_length', body_length, 'tail_length', tail_length, 'initial_box', ...
                    initial_box, 'box_increment', box_increment);
              
%% Save COMs of the single concatenated video 
% check whether COM already exists
clear; close all; clc;
videofolder = '/run/user/1000/gvfs/sftp:host=134.157.132.157/home/ljp/Zebraid2/User/Natalia/Backup_Natalia_Coolcat_/Rolling/2017-12/2017-12-13';
runs = dir([videofolder,'/run*']);
for i = 1:length(runs)
    path2video = [runs(i).folder, '/', runs(i).name, '/video_inv.avi'];
    [com1, vect0] = detectFirstCOM(path2video);
    [path, savefilename, ext] = fileparts(path2video);
    save([path, '/',savefilename, 'COM1.mat'],'com1','vect0');
end



%% Final tracking structure

dates = [19];
for d = dates
    videofolder = ['/home/ljp/Science/Projects/Sharbat/Data/2021-02-',num2str(d)];
    runs = dir([videofolder,'/Run*']);

    for i = 1:length(runs)
        runpath = [runs(i).folder, '/', runs(i).name, '/Images/'];
        path2video = [runpath,'/tailconcat.avi'];

        [path, loadfilename, ext] = fileparts(path2video);
        
        com1 = load([path, '/', loadfilename, 'COM1.mat'], 'com1'); %careful, saves a structure
        vect0 = load([path, '/', loadfilename, 'COM1.mat'], 'vect0');
        
        tracking = boutsWrapper_manual(path2video, com1.com1, vect0.vect0, 'num_segments', num_segments, 'inertia', inertia, ...
                               'body_length', body_length, 'tail_length', tail_length, 'initial_box', ...
                               initial_box, 'box_increment', box_increment);
        save([path, '/',loadfilename, '.mat'],'-struct','tracking')
  

    end
end

%% Save first COM and vectors without concatenation


% for i = 1:length(runs)
%     runpath = [runs(i).folder, '/', runs(i).name, '/Images/'];
%     vidpath = dir([runpath,'/*.avi']);
%     for k = 1:length(vidpath)
%         path2video = [vidpath(k).folder,'/',vidpath(k).name];
%         [com1, vect0] = detectFirstCOM(path2video);
%         [path, savefilename, ext] = fileparts(path2video);
%         save([path, '/',savefilename, 'COM1.mat'],'com1','vect0');
%     end
% end


%% Save video times (from pixels) and tracking structures

videofolder = '/home/ljp/Science/Projects/Sharbat/Data/2021-03-19';
runs = dir([videofolder,'/Run*']);

% save the behavioural time

vid = VideoReader(path2video);
j=1;

while hasFrame(vid) 
    im = readFrame(vid);
%    imshow(im)
    timestamp = im(1,1:4);
    tbeha(j) = fast_extract_timestamp(timestamp);
    j=j+1;
end

offset = 0;
tbeha_corr(1) = tbeha(1);
for i = 2 : length(tbeha)
    
    if tbeha(i) - tbeha(i-1) < 0
        offset = offset + 128; %128 for the maximum size of the offset
    end

    tbeha_corr(i) = tbeha(i) + offset;

end

tbeha_corr = tbeha_corr - tbeha_corr(1); %time in sec
save([path,'/tbeha.mat'],'tbeha');

for i = 1:length(runs)
    runpath = [runs(i).folder, '/', runs(i).name, '/Images/'];
    vidpath = dir([runpath,'/tailconcat.avi']);
    for k = 1:length(vidpath)
        path2video = [vidpath(k).folder,'/',vidpath(k).name];
        [path, loadfilename, ext] = fileparts(path2video);
        com1 = load([path, '/', loadfilename, 'COM1.mat'], 'com1'); %careful, saves a structure
        vect0 = load([path, '/', loadfilename, 'COM1.mat'], 'vect0');
        tracking = boutsWrapper_manual(path2video, com1.com1, vect0.vect0, 'num_segments', num_segments, 'inertia', inertia, ...
                               'body_length', body_length, 'tail_length', tail_length, 'initial_box', ...
                               initial_box, 'box_increment', box_increment);
        save([path, '/',loadfilename, '.mat'],'-struct','tracking')
    end

end

% 
% [path, savefilename, ext] = fileparts(path2video)
% save([path, '/',savefilename, '.mat'],'-struct','tracking')

% %% Save the tracking in a new video
% saveTracking_manual(path2video, 'num_segments', num_segments, 'inertia', inertia, ...
%                     'body_length', body_length, 'tail_length', tail_length, 'initial_box', ...
%                     initial_box, 'box_increment', box_increment); 
%%

% figure
% hold on
% plot(diff(tracking.total_angle) / max(abs(diff(tracking.total_angle))) * max(abs(tracking.total_angle)))
% plot(sum([tracking.Angle0(:), tracking.Angles(:, 1:9)], 2))
% legend('End of tail cumulative angle', 'mid tail cumulative angle')
% title('Cumulative angle against video frame')
% xlabel('video frame', 'Interpreter', 'latex')
% ylabel('cumulative angle [$^{\circ}$]', 'Interpreter', 'latex')



%% Save behavioural time stamps (tbeha) for already processed videos

% path2video = '/home/ljp/Science/Projects/Sharbat/Data/2020-11-20/Run 02/fc2_save_2020-11-20-131550-0000.avi';


dates = [16 17 19 22 23 24];

for d = dates
    videofolder = ['/home/ljp/Science/Projects/Sharbat/Data/2021-03-',num2str(d)];
    runs = dir([videofolder,'/Run*']);

    for i = 1:length(runs)
        runpath = [runs(i).folder, '/', runs(i).name, '/Images/'];
        path2video = [runpath,'/tailconcat.avi'];

        vid = VideoReader(path2video);
        j=1;
        tbeha = 0;
        tbeha_corr = 0;

        while hasFrame(vid) 
            im = readFrame(vid);
        %    imshow(im)
            timestamp = im(1,1:4);
            tbeha(j) = fast_extract_timestamp(timestamp);
            j=j+1;
        end

        offset = 0;
        tbeha_corr(1) = tbeha(1);
        for i = 2 : length(tbeha)

            if tbeha(i) - tbeha(i-1) < 0
                offset = offset + 128; %128 for the maximum size of the offset
            end

            tbeha_corr(i) = tbeha(i) + offset;

        end

        tbeha_corr = tbeha_corr - tbeha_corr(1); %time in sec

        [path, loadfilename, ext] = fileparts(path2video);
        save([path,'/tbeha.mat'],'tbeha_corr');
    end
end


%% testing the saving of timestamps (Volker & Sharbat 2021)



% filepath=dir([path,'/sinus*.avi'])
% 
% for i=1:length(filepath)
%     
%     path2video = [filepath(i).folder, '/', filepath(i).name]
%     vid = VideoReader(path2video);
% 
%     
%     j = 1;
% 
% 
%     while hasFrame(vid)
%     
%         im = readFrame(vid);
% %    imshow(im)
%         timestamp = im(1,1:4);
%         tbeha(j) = fast_extract_timestamp(timestamp);
%         j=j+1;
%     end
%     
%     [pks,locs] = findpeaks(tbeha)
%     
%     for k=1:length(locs)
%         
%         
%         if k==length(locs)
%             tbeha(locs(k)+1:end) = tbeha(locs(k)+1:end)+tbeha(locs(k))
%         else
%             tbeha(locs(k)+1:locs(k+1)) = tbeha(locs(k)+1:locs(k+1))+tbeha(locs(k))
%         end
%     end
%     
% [path, savefilename, ext] = fileparts(path2video)
% 
% plot(tbeha) 
% save([path,'/',savefilename,'tbeha.mat'],'tbeha')
% 
%     
% end


%%
% figure
% subplot(2,1,1)
% plot(tbeha)
% xlabel('Frame No')
% ylabel('Time (s)')
% legend('tbeha')
% subplot(2,1,2)
% plot(tbeha_corr)
% xlabel('Frame No')
% ylabel('Time (s)')
% legend('tbeha_corr')
% 

%%
%     [pks,locs] = findpeaks(tbeha)
%     
%     for k=1:length(locs)
%         
%         
%         if k==length(locs)
%             tbeha(locs(k)+1:end) = tbeha(locs(k)+1:end)+tbeha(locs(k))
%         else
%             tbeha(locs(k)+1:locs(k+1)) = tbeha(locs(k)+1:locs(k+1))+tbeha(locs(k))
%         end
%     end
%     
% 
% plot(tbeha) 


%%
% 
% offset = 0
% for i = 2 : length(tbeha)
% 
% if tbeha(i) - tbeha(i-1) <0
% 
% offset = offset + 127;
% 
% end
% 
% tbeha(i) = tbeha(i) + offset;
% 
% end


                