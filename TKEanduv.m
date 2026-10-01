clc; clear; close all;

% =========================================================================
% 0. 時間區間設定 (全域控制開關)
% =========================================================================
load('uvt_grid.mat');

is_zoom_in = true; %false true
if is_zoom_in
    t_start = datenum('2016-03-23 07:00'); 
    t_end   = datenum('2016-03-23 08:00'); 
    dt_10min = 10 / (24 * 60); 
    t_ticks = t_start : dt_10min : t_end;
else
    t_start = time(1);
    t_end   = time(end);
    t_ticks = linspace(t_start, t_end, 8);
end

idx = find(time >= t_start & time <= t_end);

% =========================================================================
% 1. 計算 Thorpe Scale 翻轉與 TKE (epsilon_raw)
% =========================================================================
[z_col, z_sort_idx] = sort(z(:), 'descend'); 
temp_sorted = temp(z_sort_idx, :);
num_z = length(z_col);
num_t_raw = size(temp_sorted, 2);
time_axis_raw = time(:)';
p_col = -z_col;
s_col = 34.5 * ones(num_z, 1);

fprintf('正在計算 Thorpe Scale 翻轉與 TKE (epsilon_raw)...\n');
epsilon_raw = zeros(num_z, num_t_raw);
opts = {'lat', 22, 'usetemp', 0, 'minotsize', 1, 'sigma', 5e-4}; 

for t = 1:num_t_raw
    t_profile = temp_sorted(:, t);
    if all(isnan(t_profile))
        epsilon_raw(:, t) = NaN;
        continue;
    end
    [eps_t, ~, ~, ~, ~] = compute_overturns(p_col, t_profile, s_col, opts{:});
    epsilon_raw(:, t) = eps_t;
end
epsilon_raw(isnan(epsilon_raw) | epsilon_raw < 1e-11) = 1e-11;

% =========================================================================
% Figure 1: ADCP 與 TKE 5 合 1 複合剖面圖
% =========================================================================
figure('Name', 'ADCP & TKE 5-Subplot Combined Field');

all_axes = zeros(5, 1);
all_cb   = zeros(5, 1);

% 預先建立 5 個獨立 axes，避免 tiledlayout 綁架位置
for i = 1:5
    all_axes(i) = axes();
end

% --- (1) v-velocity ---
ax1 = all_axes(1); axes(ax1);
contourf(time, z, v, 50, 'LineColor', 'none', 'FaceAlpha', 0.8); hold on;
contour(time, z, temp, 14:1:19, 'k', 'LineWidth', 0.8); hold off;
c1 = colorbar(ax1); all_cb(1) = c1; colormap(ax1, "jet"); caxis([-1, 1]); set(c1, 'Ticks', -1 : 1 : 1);
ylabel(c1, 'v (m/s)', 'FontSize', 9); ylabel('z (m)');

% --- (2) w-velocity ---
ax2 = all_axes(2); axes(ax2);
contourf(time, z, w, 50, 'LineColor', 'none', 'FaceAlpha', 0.8); hold on;
contour(time, z, temp, 14:1:19, 'k', 'LineWidth', 0.8); hold off;
c2 = colorbar(ax2); all_cb(2) = c2; colormap(ax2, "jet"); caxis([-0.1, 0.1]); set(c2, 'Ticks', -0.1 : 0.1 : 0.1);
ylabel(c2, 'w (m/s)', 'FontSize', 9); ylabel('z (m)');

% --- (3) Temperature ---
ax3 = all_axes(3); axes(ax3);
contourf(time, z, temp, 50, 'LineColor', 'none', 'FaceAlpha', 0.8); hold on; 
contour(time, z, temp, 14:0.25:19, 'k', 'LineWidth', 0.8); hold off;
c3 = colorbar(ax3); all_cb(3) = c3; colormap(ax3, "jet"); caxis([14, 19]);set(c3, 'Ticks', 14 : 2 : 19);
ylabel(c3, 'Temperature (\circC)', 'FontSize', 9); ylabel('z (m)');

% --- (4) Echo 1 ---
ax4 = all_axes(4); axes(ax4);
contourf(time, z, e1, 50, 'LineColor', 'none');
c4 = colorbar(ax4); all_cb(4) = c4; colormap(ax4, "parula"); caxis([80, 200]);
ylabel(c4, 'echo', 'FontSize', 9); ylabel('z (m)');

% --- (5) TKE \epsilon ---
ax5 = all_axes(5); axes(ax5);
eps_raw_log = log10(max(epsilon_raw, 1e-11));
img5 = imagesc('XData', time_axis_raw, 'YData', z_col, 'CData', eps_raw_log);
hold on;
contour(time, z, temp, 14:1:19, 'w', 'LineWidth', 0.8); hold off;
set(img5, 'AlphaData', 0.85); 
set(ax5, 'YDir', 'normal'); 
c5 = colorbar(ax5); all_cb(5) = c5; colormap(ax5, "jet"); caxis([-8, -5]); set(c5, 'Ticks', -8 : 1 : -5);
ylabel(c5, 'log_{10}(\epsilon) (W/kg)', 'FontSize', 9); 
ylabel('z (m)');

% =========================================================================
% 2. 座標軸刻度與範圍設定
% =========================================================================
y_ticks = -240 : 40 : -140;

for i = 1:5
    ax = all_axes(i);
    cb = all_cb(i);
    
    set(ax, 'TickDir', 'out', 'TickLength', [0.006, 0.01], 'Box', 'on', 'FontSize', 9);
    set(cb, 'TickDir', 'out', 'TickLength', 0.015);
    
    ylim(ax, [-241, -139]); % 略留 1m 邊界確保 -240 完全顯示
    yticks(ax, y_ticks);
    
    xlim(ax, [t_start, t_end]);             
    set(ax, 'XTick', t_ticks);         
    
    if i < 5
        xticklabels(ax, {}); 
    else
        datetick(ax, 'x', 'HH:MM', 'keeplimits', 'keepticks');
        xlabel(ax, 'Time', 'FontSize', 10, 'FontWeight', 'bold');
    end
end

linkaxes(all_axes, 'x');

% =========================================================================
% 3. ★ 🎛️ 絕對位置與數字控制區 (改數字直接生效，絕對不貼圖) ★
% =========================================================================
v_gap    = 0.020;  % ★ 上下子圖之間的間距 (數字越小間距越緊密，可調小如 0.008)
cb_gap   = 0.015;  % ★ Colorbar 與主圖右邊界的距離 (加大即可徹底遠離子圖)
cb_width = 0.008;  % ★ Colorbar 的條狀寬度 (數字越小越細)

left_m   = 0.080;  % 圖形左邊界 (留給 z (m) 標籤)
right_m  = 0.100;  % 圖形右邊界 (留給 Colorbar 標籤文字)
bottom_m = 0.070;  % 圖形底邊界 (留給 Time 標籤)
top_m    = 0.030;  % 圖形頂邊界

% 計算主圖寬度與每個子圖高度
plot_w = 1 - left_m - right_m - cb_gap - cb_width;
total_h = 1 - bottom_m - top_m;
sub_h = (total_h - (4 * v_gap)) / 5;

% 5 個 Colorbar 統一靠右對齊的 X 座標
cb_x_pos = left_m + plot_w + cb_gap;

for i = 1:5
    % 從第 1 張圖 (頂端) 排到第 5 張圖 (底端)
    y_pos = bottom_m + (5 - i) * (sub_h + v_gap);
    
    % 精確設置主圖 Position: [左, 下, 寬, 高]
    set(all_axes(i), 'Position', [left_m, y_pos, plot_w, sub_h]);
    
    % 精確設置 Colorbar Position (獨立於主圖右側，絕不蓋圖)
    set(all_cb(i), 'Position', [cb_x_pos, y_pos, cb_width, sub_h]);
end

