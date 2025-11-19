% MakeMP4_ordered.m
clear; clc;

folder = 'ex01_CH_2D_data';
% folder = 'ex16_Epitaxy_no_slpoe_epsilon0.03';
% folder = 'ex16_Epitaxy_no_slpoe_epsilon0.01';
% folder = 'ex16_Epitaxy_no_slpoe_epsilon0.005';
% folder = 'ex16_Epitaxy_no_slpoe_epsilon0.0025';
% folder = 'ex16_Epitaxy_no_slpoe_epsilon0.001';
% folder = 'ex16_Epitaxy_no_slpoe_epsilon0.0005';

pattern = fullfile(folder, '*.png');
files = dir(pattern);

% 提取文件名中的数字（t=数字部分）
timeSteps = zeros(1,numel(files));
for k = 1:numel(files)
    name = files(k).name;
    tok = regexp(name,'t=(\d+)','tokens');
    if ~isempty(tok)
        timeSteps(k) = str2double(tok{1}{1});
    else
        timeSteps(k) = NaN; % 没有匹配到数字
    end
end

% 按数字排序
[~, idx] = sort(timeSteps);
files = files(idx);

% 输出检查
disp({files.name}');

% 设置输出视频
outmp4  = [folder,'/',folder,'.mp4'];
fps = 10;
v = VideoWriter(outmp4, 'MPEG-4'); % H.264
v.FrameRate = fps;
open(v);

% 可选缩放
targetSize = [720, NaN]; % 高度 720，保持宽高比

for k = 1:numel(files)
    img = imread(fullfile(folder, files(k).name));
    if ~isempty(targetSize)
        img = imresize(img, targetSize);
    end
    % 灰度图转 3 通道
    if size(img,3) == 1
        img = repmat(img,1,1,3);
    end
    writeVideo(v, img);
end

close(v);
fprintf('MP4 written to %s\n', outmp4);