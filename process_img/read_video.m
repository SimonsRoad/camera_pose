% v = VideoReader('move.mp4');
% 
% count = 0;
% framenum = 1;
% if 0
%     while hasFrame(v)
%         frame = readFrame(v);
%         if count == 27
%             imwrite(frame, ['move_', num2str(framenum), '.png']);
%             count = 0;
%             framenum = framenum + 1;
%         end
%         count = count + 1;
%     end
% end
% 
% 
% ref0 = imread('ref.jpg') 
% ref1 = zeros(1920, 1080, 3); 


for i = 0:27
    frame_name = ['n', num2str(i), '.jpg'];
    frame = imread(frame_name);
    frame = imresize(frame, [1440, 1920]);
    frame = frame(181:180+1080,:,:);
    new_frame_name = ['n', num2str(i), '.png'];
    imwrite(frame, new_frame_name);
end