clc; clear; close all;

% =========================================================================
% 0. 載入資料與網格初始化 (確保 z 軸由淺至深排序)
% =========================================================================
load('uvt_grid.mat');

% 確保深度 z 由淺到深排序 (例如: 0, -2, -4 ... 或 0, 2, 4 ...)
[z_col, z_sort_idx] = sort(z(:), 'descend'); 
temp = temp(z_sort_idx, :);
u    = u(z_sort_idx, :);
v    = v(z_sort_idx, :);

num_z = length(z_col);
num_t_raw = size(temp, 2);
time_axis_raw = time(:)';

p_col = -z_col; % 壓力 (dbar)，取正值傳入 compute_overturns
s_col = 34.5 * ones(num_z, 1); % 假定背景鹽度

% 1. 第一階段：呼叫 compute_overturns 計算高頻 \epsilon_raw 與 Thorpe Scale
fprintf('正在計算 Thorpe Scale 翻轉與 TKE (epsilon_raw)...\n');

epsilon_raw = zeros(num_z, num_t_raw);
Lt_raw      = zeros(num_z, num_t_raw);

% 設定 compute_overturns 的參數
% sigma: 密度噪聲門檻 (CTD 通常設定 5e-4 ~ 1e-3)
% minotsize: 最小採納翻轉尺寸 (單位: m)
opts = {'lat', 22, 'usetemp', 0, 'minotsize', 1, 'sigma', 5e-4}; 

for t = 1:num_t_raw
    t_profile = temp(:, t);
    
    % 若剖面全是 NaN 則跳過
    if all(isnan(t_profile))
        epsilon_raw(:, t) = NaN;
        Lt_raw(:, t)      = NaN;
        continue;
    end
    
    % 呼叫 compute_overturns (內部會透過 nanfilt.m 自動過濾假翻轉)
    [eps_t, ~, ~, ~, Lt_t] = compute_overturns(p_col, t_profile, s_col, opts{:});
    
    epsilon_raw(:, t) = eps_t;
    Lt_raw(:, t)      = Lt_t;
end

% 背景極小值處理 (無翻轉區域設為 1e-11 W/kg)
epsilon_raw(isnan(epsilon_raw) | epsilon_raw < 1e-11) = 1e-11;

% 計算高頻背景的水文欄位 (供剖面圖第 1 子圖繪製)
p_raw_mat = repmat(p_col, 1, num_t_raw);
S_prac_raw = 34.5 * ones(num_z, num_t_raw);
lon_raw = 122 * ones(num_z, num_t_raw);
lat_raw = 22 * ones(num_z, num_t_raw);

SA_raw  = gsw_SA_from_SP(S_prac_raw, p_raw_mat, lon_raw, lat_raw);
CT_raw  = gsw_CT_from_t(SA_raw, temp, p_raw_mat);
rho_raw = gsw_rho(SA_raw, CT_raw, p_raw_mat);

% =========================================================================
% 2. 第二階段：5 分鐘 Bin 時間平均
% =========================================================================
fprintf('正在執行 5 分鐘時間 Bin 平均...\n');

dt_5min = 5 / (24 * 60); 
t_start = time_axis_raw(1);
t_end   = time_axis_raw(end);
t_edges = t_start : dt_5min : t_end;
if t_edges(end) < t_end, t_edges = [t_edges, t_end]; end

num_bins = length(t_edges) - 1;
time_5min = (t_edges(1:end-1) + t_edges(2:end)) / 2; 

u_5min       = zeros(num_z, num_bins);
v_5min       = zeros(num_z, num_bins);
temp_5min    = zeros(num_z, num_bins);
eps_5min_raw = zeros(num_z, num_bins); 

bin_idx = discretize(time_axis_raw, t_edges);

for b = 1:num_bins
    mask = (bin_idx == b);
    if any(mask)
        u_5min(:, b)    = mean(u(:, mask), 2, 'omitnan');
        v_5min(:, b)    = mean(v(:, mask), 2, 'omitnan');
        temp_5min(:, b) = mean(temp(:, mask), 2, 'omitnan');
        
        % TKE 做 Log 空間平均 (幾何平均)，能更精確呈現湍流能量的代表值
        eps_block = epsilon_raw(:, mask);
        eps_5min_raw(:, b) = 10.^(mean(log10(max(eps_block, 1e-11)), 2, 'omitnan'));
    else
        u_5min(:, b)       = NaN;
        v_5min(:, b)       = NaN;
        temp_5min(:, b)    = NaN;
        eps_5min_raw(:, b) = NaN;
    end
end

% =========================================================================
% 3. 第三階段：計算 5 分鐘平均下之 S^2, N^2 與深度網格精確對齊
% =========================================================================
fprintf('正在計算平滑流場剪切 S^2 與浮力頻率 N^2...\n');

p_5min = repmat(p_col, 1, num_bins);
S_prac_5min = 34.5 * ones(num_z, num_bins);
lon_5min = 122 * ones(num_z, num_bins);
lat_5min = 22 * ones(num_z, num_bins);

SA_5min = gsw_SA_from_SP(S_prac_5min, p_5min, lon_5min, lat_5min);
CT_5min = gsw_CT_from_t(SA_5min, temp_5min, p_5min);

% 計算 5 分鐘平均下的 N^2 與中間層網格 (z_mid_axis)
[N2_5min, p_mid_5min] = gsw_Nsquared(SA_5min, CT_5min, p_5min, lat_5min);
z_mid_axis = -p_mid_5min(:, 1); % 中間層深度網格

% 對流場做垂直適度平滑後計算 S^2
u_smooth = smoothdata(u_5min, 1, 'movmean', 3);
v_smooth = smoothdata(v_5min, 1, 'movmean', 3);

du = diff(u_smooth, 1, 1);
dv = diff(v_smooth, 1, 1);
dz_raw = abs(diff(z_col)); % 確保深度差為正值
dz = repmat(dz_raw, 1, num_bins); 

du_dz = du ./ dz;
dv_dz = dv ./ dz;

% S2 的深度網格天然對齊在 z_mid_axis 上
S2_5min = du_dz.^2 + dv_dz.^2;               
shear_instability_5min = S2_5min - 4 * N2_5min; 

% 將 5 分鐘 TKE 插值至 z_mid_axis 中間層網格，達成完全對齊
eps_5min_mid = zeros(size(S2_5min));
for b = 1:num_bins
    eps_5min_mid(:, b) = interp1(z_col, eps_5min_raw(:, b), z_mid_axis, 'linear', 'extrap');
end

% =========================================================================
% 4. 繪圖一：時空剖面圖 (Physics Fields & TKE Profile)
% =========================================================================
% =========================================================================
% 4. 繪圖一：時空剖面圖 (Physics Fields & TKE Profile)
% =========================================================================
figure('Position', [100, 100, 1500, 1200], 'Name', 'Physics Fields & TKE Profiles'); 

% 建立包含精確頭尾的 8 個均勻時間刻度
t_min = time_axis_raw(1);
t_max = time_axis_raw(end);
t_ticks_custom = linspace(t_min, t_max, 8); 

% (1) \rho (原始高頻)
subplot(5, 1, 1);
levels_rho = 1025:0.1:1027; 
contourf(time_axis_raw, z_col, rho_raw, levels_rho, 'LineColor', 'k', 'LineWidth', 0.2); 
set(gca, 'YDir', 'normal');
c1 = colorbar; colormap(gca, "jet"); caxis([1025, 1027]); 
ylabel(c1, '\rho (kg/m^3)', 'FontSize', 8); ylabel('z (m)'); 

% (2) N^2 (5 分鐘)
subplot(5, 1, 2);
pcolor(time_5min, z_mid_axis, log10(max(N2_5min, 1e-10))); shading interp; set(gca, 'YDir', 'normal');
c2 = colorbar; colormap(gca, "jet"); caxis([-5, -2.75]); 
ylabel(c2, 'log_{10} N^2 (s^{-2})', 'FontSize', 8); ylabel('z (m)');

% (3) S^2 (5 分鐘)
subplot(5, 1, 3);
pcolor(time_5min, z_mid_axis, log10(max(S2_5min, 1e-10))); shading interp; set(gca, 'YDir', 'normal');
c3 = colorbar; colormap(gca, "jet"); caxis([-5, -2.75]); 
ylabel(c3, 'log_{10} S^2 (s^{-2})', 'FontSize', 8); ylabel('z (m)'); 

% (4) S^2 - 4N^2 (5 分鐘)
subplot(5, 1, 4);
pcolor(time_5min, z_mid_axis, shear_instability_5min); shading interp; set(gca, 'YDir', 'normal');
c4 = colorbar; colormap(gca, "jet"); caxis([-5e-4, 5e-4]); 
ylabel(c4, 'S^2 - 4N^2 (s^{-2})', 'FontSize', 8); ylabel('z (m)'); 

% (5) TKE \epsilon (原始高頻 compute_overturns 計算結果)
subplot(5, 1, 5);
eps_raw_log = log10(max(epsilon_raw, 1e-11));
img5 = imagesc(time_axis_raw, z_col, eps_raw_log);
set(img5, 'AlphaData', 0.85); set(gca, 'YDir', 'normal');
c5 = colorbar; colormap(gca, "jet"); caxis([-8.5, -5]); 
ylabel(c5, 'log_{10}(\epsilon) (W/kg)', 'FontSize', 8); 
xlabel('Time'); ylabel('z (m)'); 

% ★ 統一強制鎖定 5 張子圖的時間軸頭尾範圍與刻度 ★
for i = 1:5
    subplot(5, 1, i);
    xlim([t_min, t_max]);                       % 1. 強制設定 X 軸邊界為真正的頭尾
    set(gca, 'XTick', t_ticks_custom);          % 2. 指定包含頭尾的時間刻度
    datetick('x', 'mm-dd HH:MM', 'keeplimits', 'keepticks'); % 3. 轉為時間格式並保留邊界與刻度
end
% =========================================================================
% 5. 繪圖二：散點圖 (Richardson Number Scatter Plot)
% =========================================================================
x_data = S2_5min(:);
y_data = 4 * N2_5min(:);
eps_data = eps_5min_mid(:); 

% 過濾無效與極端異常值 (保留有效數據點)
valid_idx = (x_data > 0) & (y_data > 0) & (eps_data > 1e-11) & ~isnan(x_data) & ~isnan(y_data);

x_plot = x_data(valid_idx);
y_plot = y_data(valid_idx);
c_plot = log10(eps_data(valid_idx)); 

% 按 TKE 強弱排序（高 TKE 暖色點最後畫，避免圖層遮擋）
[c_plot, sort_order] = sort(c_plot, 'ascend');
x_plot = x_plot(sort_order);
y_plot = y_plot(sort_order);

figure('Name', 'Richardson Number Scatter Plot', 'Position', [200, 200, 700, 600]);
scatter(x_plot, y_plot, 14, c_plot, 'filled', 'MarkerFaceAlpha', 0.65);
hold on;

% 1:1 臨界線 (Ri = 0.25 Line)
plot([1e-8, 1], [1e-8, 1], '--k', 'LineWidth', 1.5); 

set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([1e-8, 1e-1]); ylim([1e-8, 1e-1]); % 聚焦在主要數據集中區
grid on; box on;

xlabel('S^2 (s^{-2})', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('4N^2 (s^{-2})', 'FontSize', 11, 'FontWeight', 'bold');

c_scatter = colorbar; colormap(gca, 'jet'); caxis([-8, -5]); 
ylabel(c_scatter, 'log_{10}(\epsilon) (W/kg)', 'FontSize', 10, 'Rotation', 90);

text(1e-4, 5e-7, 'Unstable Region (Ri < 0.25)', 'FontSize', 10, 'FontWeight', 'bold');
text(1e-5, 5e-3, 'Stable Region (Ri >= 0.25)', 'FontSize', 10, 'FontWeight', 'bold');
hold off;

fprintf('處理完成！\n');