% ПАРАМЕТРЫ МОДЕЛИ
% ЗНАЧЕНИЯ R0 ДЛЯ АНАЛИЗА
R0_list = [4, 8, 12];  

% Фиксированное значение w 
w_value = 10; 

% Параметры модели
I0E = 0.03;
I0H = 9e-6;
time_total = 100.0;

% Параметры латентности
r0 = 0.001;
pL0 = 0.4;
E0L = 0.01;

% Временные параметры
c = 2.0;
dT = 0.3;
dI = 1.0;
dE = 0.15;
dH = 0.15;
rmax = 0.2 * dI;

% Масштабирование
T0 = 1.0;
b = T0 * dT;
H0 = 0.01;
E0 = 0.01;
k = dI / E0;
Em = 1.0;

% Начальные условия
Hin = 1e-5;
Ein = 1e-5;
Iin = 1e-11 * T0;
Lin = Iin;

% Эффективность терапии
effic = 0.97;

% Диапазоны для tARTstart и tARTend
tARTstart_range = 0.5:0.25:10;
tARTend_range = 4:0.5:20;

% Порог для хелперов
H_threshold = 1e-9;

% Массивы для хранения результатов
points_cell = cell(length(R0_list), 1);

for idx = 1:length(R0_list)
    points_cell{idx} = [];
end

fprintf('=== АНАЛИЗ ДЛЯ R0 = [%d, %d, %d], w = %d ===\n\n', R0_list(1), R0_list(2), R0_list(3), w_value);


for R_idx = 1:length(R0_list)
    current_R0 = R0_list(R_idx);
    fprintf('Обработка R0 = %d...\n', current_R0);
    
    % Вычисляем p для текущего R0
    p = current_R0 * dT * dI / b;
    
    % Начальные условия
    y0 = [T0, Iin, Lin, Hin, Ein];
    
    % Структура параметров
    params.R0 = current_R0;
    params.I0H = I0H;
    params.I0E = I0E;
    params.time = time_total;
    params.r0 = r0;
    params.pL0 = pL0;
    params.E0L = E0L;
    params.c = c;
    params.dT = dT;
    params.dI = dI;
    params.dE = dE;
    params.dH = dH;
    params.rmax = rmax;
    params.b = b;
    params.p = p;
    params.H0 = H0;
    params.E0 = E0;
    params.k = k;
    params.Em = Em;
    params.Iin = Iin;
    params.y0 = y0;
    params.effic = effic;
    
    points_count = 0;
    
    % Цикл по всем комбинациям
    for i = 1:length(tARTstart_range)
        for j = 1:length(tARTend_range)
            t_start = tARTstart_range(i);
            t_end = tARTend_range(j);
            
            if t_end <= t_start
                continue;
            end
            
            % Проверка спонтанного контроля
            [spontaneous, ~, ~] = check_spontaneous(t_start, t_end, params, 1.0);
            
            if spontaneous
                continue;
            end
            
            % Проверка посттерапевтического контроля
            I_lower = compute_I_lower(params, w_value);
            I_upper = compute_I_upper(params);
            
            therapy = check_therapy(t_start, t_end, params, w_value, I_lower, I_upper, H_threshold);
            
            if therapy
                points_cell{R_idx} = [points_cell{R_idx}; t_start, t_end];
                points_count = points_count + 1;
            end
        end
    end
    
    fprintf('  Найдено точек: %d\n', points_count);
end

fprintf('\nАнализ завершен успешно!\n');


% ГРАФИКИ
figure('Position', [100, 100, 1200, 400]);

for R_idx = 1:length(R0_list)
    current_R0 = R0_list(R_idx);
    
    subplot(1, 3, R_idx);
    hold on;
    grid on;
    
    title(sprintf('R0 = %d I0E = %.1e IOH = %.1e', current_R0, I0E, I0H), 'FontSize', 14);
    xlabel('tARTend', 'FontSize', 12);
    ylabel('tARTstart', 'FontSize', 12);
    
    if ~isempty(points_cell{R_idx})
        plot(points_cell{R_idx}(:,2), points_cell{R_idx}(:,1), 'r+', 'MarkerSize', 6);
        text(15, 1, sprintf('Точек: %d', size(points_cell{R_idx}, 1)), 'FontSize', 10);
    else
        text(12, 5, 'Нет точек', 'HorizontalAlignment', 'center');
    end
    
    xlim([3, 21]);
    ylim([0, 10.5]);
    hold off;
end

sgtitle(sprintf('w = %d', w_value), 'FontSize', 16);

% ФУНКЦИИ
function dxdt = model_no_therapy(t, x, R0, I0H, I0E, w, r0, pL0, E0L, c, dT, dI, dE, dH, rmax, b, p, H0, E0, k, Em, Iin)
    T = max(x(1), 0); I = max(x(2), 0); L = max(x(3), 0); H = max(x(4), 0); E = max(x(5), 0);
    
    if I > 0 && I0H > 0
        alphaH = 1 - exp(-I / I0H);
    else
        alphaH = 0;
    end
    
    if Em > 0
        alphaH = alphaH * max(0, min(1, (1 - H / Em)))^(1/3);
    else
        alphaH = 0;
    end
    
    if H > 0 && H0 > 0 && I > 0 && I0E > 0
        sigma = 1 - exp(-alphaH * H / H0 - I / I0E);
    else
        sigma = 0;
    end
    
    if Em > 0
        sigma = sigma * max(0, min(1, (1 - E / Em)))^(1/3);
    else
        sigma = 0;
    end
    
    if E > 0 && E0L > 0
        r = r0 + rmax * E / (E0L + E);
    else
        r = r0;
    end
    
    if E > 0 && E0L > 0
        pL = pL0 * E0L / (E0L + E);
    else
        pL = pL0;
    end
    
    if I > Iin / 2
        p1 = p;
    else
        p1 = 0;
    end
    
    if H <= Iin / 2
        alphaH = 0;
    end
    
    if E <= Iin / 2
        sigma = 0;
    end
    
    dTdt = b - p1 * I * T - dT * T;
    dIdt = ((1 - pL) * p1 * I * T - k * I * E - dI * I + r * L);
    dLdt = pL * p1 * I * T - r * L;
    dHdt = (c * alphaH * H - p1 * w * I * H - dH * H);
    dEdt = (c * sigma * E - dE * E);
    
    dxdt = [dTdt; dIdt; dLdt; dHdt; dEdt];
end

function dxdt = model_with_therapy(t, x, R0, I0H, I0E, w, t_start, t_end, effic, r0, pL0, E0L, c, dT, dI, dE, dH, rmax, b, p, H0, E0, k, Em, Iin)
    T = max(x(1), 0); I = max(x(2), 0); L = max(x(3), 0); H = max(x(4), 0); E = max(x(5), 0);
    
    if I > 0 && I0H > 0
        alphaH = 1 - exp(-I / I0H);
    else
        alphaH = 0;
    end
    
    if Em > 0
        alphaH = alphaH * max(0, min(1, (1 - H / Em)))^(1/3);
    else
        alphaH = 0;
    end
    
    if H > 0 && H0 > 0 && I > 0 && I0E > 0
        sigma = 1 - exp(-alphaH * H / H0 - I / I0E);
    else
        sigma = 0;
    end
    
    if Em > 0
        sigma = sigma * max(0, min(1, (1 - E / Em)))^(1/3);
    else
        sigma = 0;
    end
    
    if E > 0 && E0L > 0
        r = r0 + rmax * E / (E0L + E);
    else
        r = r0;
    end
    
    if E > 0 && E0L > 0
        pL = pL0 * E0L / (E0L + E);
    else
        pL = pL0;
    end
    
    if I > Iin / 2
        p1 = p;
    else
        p1 = 0;
    end
    
    if H <= Iin / 2
        alphaH = 0;
    end
    
    if E <= Iin / 2
        sigma = 0;
    end
    
    if t_start < t && t < t_end
        p1 = p1 * (1 - effic);
    end
    
    dTdt = b - p1 * I * T - dT * T;
    dIdt = ((1 - pL) * p1 * I * T - k * I * E - dI * I + r * L);
    dLdt = pL * p1 * I * T - r * L;
    dHdt = (c * alphaH * H - p1 * w * I * H - dH * H);
    dEdt = (c * sigma * E - dE * E);
    
    dxdt = [dTdt; dIdt; dLdt; dHdt; dEdt];
end

function [spontaneous, I_final, H_final] = check_spontaneous(t_start, t_end, params, w)
    R0 = params.R0; I0H = params.I0H; I0E = params.I0E;
    r0 = params.r0; pL0 = params.pL0; E0L = params.E0L;
    c = params.c; dT = params.dT; dI = params.dI; dE = params.dE; dH = params.dH;
    rmax = params.rmax; b = params.b; p = params.p; H0 = params.H0; E0 = params.E0;
    k = params.k; Em = params.Em; Iin = params.Iin;
    
    tspan = [0, t_start];
    y0 = params.y0;
    
    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-6, 'MaxStep', 0.1);
    
    try
        [~, x] = ode15s(@(t,x) model_no_therapy(t, x, R0, I0H, I0E, w, r0, pL0, E0L, c, dT, dI, dE, dH, rmax, b, p, H0, E0, k, Em, Iin), ...
                        tspan, y0, opts);
        
        I_final = max(x(end, 2), 0);
        H_final = max(x(end, 4), 0);
    catch
        I_final = 1e10;
        H_final = 0;
    end
    
    if H_final > 10 * I_final && I_final < 1e-3
        spontaneous = true;
    else
        spontaneous = false;
    end
end

function therapy = check_therapy(t_start, t_end, params, w, I_low, I_up, H_thresh)
    R0 = params.R0; I0H = params.I0H; I0E = params.I0E;
    r0 = params.r0; pL0 = params.pL0; E0L = params.E0L;
    c = params.c; dT = params.dT; dI = params.dI; dE = params.dE; dH = params.dH;
    rmax = params.rmax; b = params.b; p = params.p; H0 = params.H0; E0 = params.E0;
    k = params.k; Em = params.Em; Iin = params.Iin;
    effic = params.effic; time_total = params.time;
    
    tspan = [0, time_total];
    y0 = params.y0;
    
    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-6, 'MaxStep', 0.1);
    
    try
        [~, x] = ode15s(@(t,x) model_with_therapy(t, x, R0, I0H, I0E, w, t_start, t_end, effic, r0, pL0, E0L, c, dT, dI, dE, dH, rmax, b, p, H0, E0, k, Em, Iin), ...
                        tspan, y0, opts);
        
        I_final = max(x(end, 2), 0);
        H_final = max(x(end, 4), 0);
    catch
        therapy = false;
        return;
    end
    
    if H_final > H_thresh && I_final > I_low && I_final < I_up
        therapy = true;
    else
        therapy = false;
    end
end

function I_up = compute_I_upper(params)
    c = params.c;
    dE = params.dE;
    I0E = params.I0E;
    
    I_up = I0E * log(c/(c - dE));
end

function I_low = compute_I_lower(params, w)
    c = params.c;
    dH = params.dH;
    I0H = params.I0H;
    p = params.p;
    
    I_guess = I0H;
    max_iter = 1000;
    tolerance = 1e-8;
    
    for iter = 1:max_iter
        argument = c/(c - dH - p*w*I_guess);
        if argument > 0
            I_new = I0H * log(argument);
        else
            I_new = 0;
        end
        
        if abs(I_new - I_guess) < tolerance
            I_low = I_new;
            return;
        end
        
        I_guess = I_new;
    end
    
    I_low = I_guess;
end
