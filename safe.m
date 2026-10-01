clc; clear; close all;

filename = 'SAFE_GPS.xlsx'; 
sheet1_name = 'Rover-21390';  
sheet2_name = 'Apollo-00110'; 
opts = detectImportOptions(filename);
opts.VariableNamingRule = 'preserve';

% SAFE 投放點與時間
deploy_lon = 117.350000;
deploy_lat = 21.098185;
deploy_time_str = '2026-07-28 11:59'; % 投放時間

% Rover
opts.Sheet = sheet1_name;
data1 = readtable(filename, opts);
lat1_raw  = data1.Latitude;   
lon1_raw  = data1.Longitude;  
time1_raw = datetime(data1.('Time (mm-dd HH:MM)')); 

% Apollo
opts.Sheet = sheet2_name;
data2 = readtable(filename, opts);
lat2_raw  = data2.Latitude;   
lon2_raw  = data2.Longitude;  
time2_raw = datetime(data2.('Time (mm-dd HH:MM)')); 

% 2. 納入投放點並計算區間速度
deploy_time = datetime(deploy_time_str);

% 將投放點插入到第 1 筆
lat1  = [deploy_lat; lat1_raw];   lon1  = [deploy_lon; lon1_raw];   time1 = [deploy_time; time1_raw];
lat2  = [deploy_lat; lat2_raw];   lon2  = [deploy_lon; lon2_raw];   time2 = [deploy_time; time2_raw];

% 呼叫 sw_dist 計算兩點之間的真實距離 (km) 與速度 (m/s)
dist_km1 = sw_dist(lat1, lon1, 'km');
dt_sec1  = seconds(diff(time1));
speed1   = (dist_km1 * 1000) ./ dt_sec1;

dist_km2 = sw_dist(lat2, lon2, 'km');
dt_sec2  = seconds(diff(time2));
speed2   = (dist_km2 * 1000) ./ dt_sec2;

% 時間點使用前後兩點的中點，代表區間平均速度的時間
time_speed1 = time1(1:end-1) + diff(time1)/2;
time_speed2 = time2(1:end-1) + diff(time2)/2;

% 3. 儀表板印出兩台 GPS 的最新狀態 (含度分秒顯示)
fprintf(' SAFE 投放點座標 : %s, %s (投放時間: %s)\n', ...
        deg2dms_str(deploy_lat, 'lat'), deg2dms_str(deploy_lon, 'lon'), datestr(deploy_time, 'yyyy-mm-dd HH:MM'));
fprintf('----------------------------------------------------------\n');
fprintf('【(%s)】\n', sheet1_name);
fprintf('  最新時間 : %s\n', datestr(time1(end), 'yyyy-mm-dd HH:MM'));
fprintf('  最新座標 : %s, %s (%.6f°N, %.6f°E)\n', ...
        deg2dms_str(lat1(end), 'lat'), deg2dms_str(lon1(end), 'lon'), lat1(end), lon1(end));
fprintf('  最新速度 : %.2f m/s\n', speed1(end));
fprintf('----------------------------------------------------------\n');
fprintf('【(%s)】\n', sheet2_name);
fprintf('  最新時間 : %s\n', datestr(time2(end), 'yyyy-mm-dd HH:MM'));
fprintf('  最新座標 : %s, %s (%.6f°N, %.6f°E)\n', ...
        deg2dms_str(lat2(end), 'lat'), deg2dms_str(lon2(end), 'lon'), lat2(end), lon2(end));
fprintf('  最新速度 : %.2f m/s\n', speed2(end));

% 4. 圖一：軌跡對比圖 (乾淨無地形底圖)
% figure('Name', 'SAFE 雙 GPS 漂流軌跡圖', 'Position', [100, 100, 700, 600]);
figure(1);

% 繪制 GPS 1 與 GPS 2 的漂流軌跡
plot(lon1_raw, lat1_raw, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
plot(lon2_raw, lat2_raw, 'm-s', 'LineWidth', 1.2, 'MarkerSize', 4);

% 標示 SAFE 投放點與最新位置
plot(deploy_lon, deploy_lat, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8); 
plot(lon1_raw(end), lat1_raw(end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8); 
plot(lon2_raw(end), lat2_raw(end), 'rx', 'LineWidth', 2, 'MarkerSize', 10); 

% 投放點與第一點連線
plot([deploy_lon, lon1_raw(1)], [deploy_lat, lat1_raw(1)], 'b--', 'LineWidth', 1.2);
plot([deploy_lon, lon2_raw(1)], [deploy_lat, lat2_raw(1)], 'm--', 'LineWidth', 1.2);

grid on; axis equal;
xlabel('Longitude'); 
ylabel('Latitude');
legend('Rover-21390 軌跡', 'Apollo-00110軌跡','SAFE 投放點', 'Rover-21390 now', 'Apollo-00110 now', 'Location', 'best');

% 💡 將軌跡圖的 X/Y 軸刻度自動轉為度分秒 (DMS)
ax = gca;
xticks_val = ax.XTick;
yticks_val = ax.YTick;
xtick_labels = cell(size(xticks_val));
for i = 1:length(xticks_val)
    xtick_labels{i} = deg2dms_str(xticks_val(i), 'lon');
end
ytick_labels = cell(size(yticks_val));
for i = 1:length(yticks_val)
    ytick_labels{i} = deg2dms_str(yticks_val(i), 'lat');
end
ax.XTickLabel = xtick_labels;
ax.YTickLabel = ytick_labels;

% 5. 圖二：速度隨時間變化圖 (獨立視窗)
% figure('Name', 'SAFE 雙 GPS 移動速度變化圖', 'Position', [820, 100, 700, 500]);
figure(2);
plot(time_speed1, speed1, '-bo', 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
plot(time_speed2, speed2, '--ms', 'LineWidth', 1.5, 'MarkerSize', 4);
grid on;
xlabel('Time(UTC)');
ylabel('m/s');
legend('Rover-21390', 'Apollo-00110', 'Location', 'best');
xtickformat('MM-dd HH:mm');

% 6. 轉換成「度分秒 (DMS)」的輔助函式 (修正四捨五入滿60秒自動進位問題)
function dms_str = deg2dms_str(deg, type)
    abs_deg = abs(deg);
    d = floor(abs_deg);
    m = floor((abs_deg - d) * 60);
    s = ((abs_deg - d) * 60 - m) * 60;
    
    % 四捨五入處理，防止出現 60.00"
    if round(s, 2) >= 60
        s = 0;
        m = m + 1;
        if m >= 60
            m = 0;
            d = d + 1;
        end
    end
    
    if strcmp(type, 'lat')
        dir = 'N'; if deg < 0, dir = 'S'; end
    else
        dir = 'E'; if deg < 0, dir = 'W'; end
    end
    
    dms_str = sprintf('%d°%02d''%05.2f"%s', d, m, s, dir);
end