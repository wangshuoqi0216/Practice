clear;
clc;
close all;

% add path
addpath('../','-begin');

% 设置数据目录
dirname = 'ex02_CH_2D_data';
datadir = [dirname, '/data'];

figdir  = [dirname,'/',dirname];
data_files = dir([datadir, '/*.txt']);
fprintf('找到 %d 个数据文件\n', length(data_files));

% 计算域
domain.left = -pi;
domain.right = pi;
domain.bottom = -pi;
domain.top = pi;

% 设置网格点
N = 64;
Nx = 2*N; 
Ny = 2*N;
x = domain.left + (domain.right-domain.left)/Nx*(0:Nx-1);
y = domain.bottom + (domain.top-domain.bottom)/Ny*(0:Ny-1);
[xx, yy] = ndgrid(x, y);

% 检测时间点
time_points = [];
for i = 1:length(data_files)
    filename = data_files(i).name;
    [match, tokens] = regexp(filename, 'phi1_t=([0-9.]+)\.txt', 'match', 'tokens');
    if ~isempty(match)
        time_points(end+1) = str2double(tokens{1}{1});
    end
end

time_points = sort(time_points);
fprintf('时间点: %s\n', mat2str(time_points));

% 创建动态图
figure(1);
set(gcf, 'Position', [100, 100, 800, 600]);

for k = 1:length(time_points)
    t = time_points(k);
    filename = [datadir '/phi_t=' num2str(t)];
    % 读取数据
    data_file = [datadir, '/phi1_t=', num2str(t), '.txt'];
    if ~exist(data_file, 'file')
        data_file = [datadir, '/', data_files(k).name];
    end
    
    phi = load(data_file);
    phi_2d = reshape(phi,Nx,Ny);
    
    % 绘制彩色图
    pcolor(xx', yy', phi_2d');
    shading interp;
    colormap(jet);
    colorbar;
    clim([-1, 1]);
     %caxis([-1.5, 1.5]);
    
    axis equal;
    xlim([domain.left, domain.right]);
    ylim([domain.bottom, domain.top]);
    
    xlabel('x');
    ylabel('y');
    title(['CahnHilliard Equation: t = ', num2str(t)]);
    
    drawnow;
    
    % 保存PNG图片
    figname = [figdir '_phi_t=' num2str(t) '.png'];
    saveas(gcf, figname);
   % fprintf('已保存图片: %s\n', figname);
    
    pause(0.1); % 控制动画的速度
end

