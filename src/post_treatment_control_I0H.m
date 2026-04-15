% ПАРАМЕТРЫ МОДЕЛИ

R0 = 8.0;
I0E = 0.03;  


I0H_values = [1e-6, 1e-5, 1e-4, 1e-3, 3e-2, 1.0];

time = 100.0;

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
p = R0 * dT * dI / b;
H0 = 0.01;
E0 = 0.01;
k = dI / E0;
Em = 1.0;

% Начальные условия
Hin = 1e-5;
Ein = 1e-5;
Iin = 1e-11 * T0;
Lin = Iin;
y0 = [T0, Iin, Lin, Hin, Ein];

% Эффективность терапии
effic = 0.97;

% АЛГОРИТМ ПО ТОЧКАМ
% Значение w для спонтанного контроля
w_spontaneous = 1.0;

% Диапазоны для tARTstart и tARTend
tARTstart_vec = 0.5:0.25:10;
tARTend_vec = 4:0.5:20;

% Порог для хелперов
H_threshold = 1e-9;

% Массивы для хранения результатов
spontaneous_points = []; 
therapy_points = cell(length(I0H_values), 1); 

for I0H_idx = 1:length(I0H_values)
    therapy_points{I0H_idx} = [];
end


% ОСНОВНОЙ ЦИКЛ ПО ТОЧКАМ
fprintf('Анализ точек для разных I0H (всего %d значений)...\n', length(I0H_values));
fprintf('Значения I0H: 1e-6, 1e-5, 1e-4, 1e-3, 3e-2, 1.0\n');

base_params.R0 = R0;
base_params.I0E = I0E;
base_params.time = time;
base_params.r0 = r0;
base_params.pL0 = pL0;
base_params.E0L = E0L;
base_params.c = c;
base_params.dT = dT;
base_params.dI = dI;
base_params.dE = dE;
base_params.dH = dH;
base_params.rmax = rmax;
base_params.b = b;
base_params.p = p;
base_params.H0 = H0;
base_params.E0 = E0;
base_params.k = k;
base_params.Em = Em;
base_params.Iin = Iin;
base_params.y0 = y0;
base_params.effic = effic;

for i = 1:length(tARTstart_vec)
    for j = 1:length(tARTend_vec)
        tARTstart_val = tARTstart_vec(i);
        tARTend_val = tARTend_vec(j);
        
        if tARTend_val <= tARTstart_val
            continue;
        end
        
        % ШАГ 1: Проверка на спонтанный контроль
        params_spontaneous = base_params;
        params_spontaneous.I0H = 0.03;
        
        [is_spontaneous, I_final, H_final] = check_spontaneous_control_extended(...
            tARTstart_val, tARTend_val, params_spontaneous, w_spontaneous);
        
        if is_spontaneous
            spontaneous_points = [spontaneous_points; tARTstart_val, tARTend_val, I_final, H_final];
            continue;
        end
        
        
        for I0H_idx = 1:length(I0H_values)
            current_I0H = I0H_values(I0H_idx);
            
            params = base_params;
            params.I0H = current_I0H;
            params.I0E = I0E;
            
            I1 = compute_lower_bound(params, 10);
            I2 = compute_upper_bound(params);
            
            is_therapy_control = check_post_therapy_control_with_bounds(...
                tARTstart_val, tARTend_val, params, 10, I1, I2, H_threshold);
            
            if is_therapy_control
                therapy_points{I0H_idx} = [therapy_points{I0H_idx}; tARTstart_val, tARTend_val];
            end
        end
    end
end

fprintf('Анализ завершен.\n');



if ~isempty(spontaneous_points)
    save_spontaneous_therapy_points(spontaneous_points);
    fprintf('Точки спонтанного контроля сохранены в файл\n');
else
    fprintf('Нет точек спонтанного контроля\n');
end



figure('Position', [50, 50, 1600, 900]);


titles = {
    sprintf('I0H = 1e-6'), ...
    sprintf('I0H = 1e-5'), ...
    sprintf('I0H = 1e-4'), ...
    sprintf('I0H = 1e-3'), ...
    sprintf('I0H = 3e-2'), ...
    sprintf('I0H = 1.0')
};

for I0H_idx = 1:length(I0H_values)
    current_I0H = I0H_values(I0H_idx);
    
    subplot(2, 3, I0H_idx);
    hold on;
    grid on;
    
    title(titles{I0H_idx}, 'FontSize', 12);
    xlabel('tARTend (days)', 'FontSize', 10);
    ylabel('tARTstart (days)', 'FontSize', 10);
    
    if ~isempty(therapy_points{I0H_idx})
        
        plot(therapy_points{I0H_idx}(:,2), therapy_points{I0H_idx}(:,1), ...
            'r+', 'MarkerSize', 8, 'LineWidth', 1.5);
        

    else
        text(mean(xlim), mean(ylim), 'Нет точек контроля', ...
            'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', 'black');
    end
    
    xlim([min(tARTend_vec)-0.5, max(tARTend_vec)+0.5]);
    ylim([min(tARTstart_vec)-0.25, max(tARTstart_vec)+0.25]);
    set(gca, 'FontSize', 9);
    
    hold off;
end

sgtitle(sprintf('Посттерапевтический контроль для I0H (I0E=%.2f, w=10, R0=%d)', I0E, R0), ...
    'FontSize', 14, 'FontWeight', 'bold');


% ФУНКЦИИ 

function dxdt = odefun_no_therapy(t, x, R0, I0H, I0E, w, r0, pL0, E0L, c, dT, dI, dE, dH, rmax, b, p, H0, E0, k, Em, Iin)
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

function dxdt = odefun_with_therapy(t, x, R0, I0H, I0E, w, tARTstart, tARTend, effic, r0, pL0, E0L, c, dT, dI, dE, dH, rmax, b, p, H0, E0, k, Em, Iin)
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
    
    if tARTstart < t && t < tARTend
        p1 = p1 * (1 - effic);
    end
    
    dTdt = b - p1 * I * T - dT * T;
    dIdt = ((1 - pL) * p1 * I * T - k * I * E - dI * I + r * L);
    dLdt = pL * p1 * I * T - r * L;
    dHdt = (c * alphaH * H - p1 * w * I * H - dH * H);
    dEdt = (c * sigma * E - dE * E);
    
    dxdt = [dTdt; dIdt; dLdt; dHdt; dEdt];
end

function [is_spontaneous, I_final, H_final] = check_spontaneous_control_extended(tARTstart, tARTend, params, w)
    R0 = params.R0; I0H = params.I0H; I0E = params.I0E;
    r0 = params.r0; pL0 = params.pL0; E0L = params.E0L;
    c = params.c; dT = params.dT; dI = params.dI; dE = params.dE; dH = params.dH;
    rmax = params.rmax; b = params.b; p = params.p; H0 = params.H0; E0 = params.E0;
    k = params.k; Em = params.Em; Iin = params.Iin;
    
    tspan = [0, tARTstart];
    y0 = params.y0;
    
    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-6, 'MaxStep', 0.1);
    
    try
        [t, x] = ode15s(@(t,x) odefun_no_therapy(t, x, R0, I0H, I0E, w, r0, pL0, E0L, c, dT, dI, dE, dH, rmax, b, p, H0, E0, k, Em, Iin), ...
                        tspan, y0, opts);
        
        I_final = max(x(end, 2), 0);
        H_final = max(x(end, 4), 0);
    catch
        I_final = 1e10;
        H_final = 0;
    end
    
    if H_final > 10 * I_final && I_final < 1e-3
        is_spontaneous = true;
    else
        is_spontaneous = false;
    end
end

function is_therapy_control = check_post_therapy_control_with_bounds(tARTstart, tARTend, params, w, I1, I2, H_threshold)
    R0 = params.R0; I0H = params.I0H; I0E = params.I0E;
    r0 = params.r0; pL0 = params.pL0; E0L = params.E0L;
    c = params.c; dT = params.dT; dI = params.dI; dE = params.dE; dH = params.dH;
    rmax = params.rmax; b = params.b; p = params.p; H0 = params.H0; E0 = params.E0;
    k = params.k; Em = params.Em; Iin = params.Iin;
    effic = params.effic; time = params.time;
    
    tspan = [0, time];
    y0 = params.y0;
    
    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-6, 'MaxStep', 0.1);
    
    try
        [t, x] = ode15s(@(t,x) odefun_with_therapy(t, x, R0, I0H, I0E, w, tARTstart, tARTend, effic, r0, pL0, E0L, c, dT, dI, dE, dH, rmax, b, p, H0, E0, k, Em, Iin), ...
                        tspan, y0, opts);
        
        I_final = max(x(end, 2), 0);
        H_final = max(x(end, 4), 0);
    catch
        is_therapy_control = false;
        return;
    end
    
    if H_final > H_threshold && I_final > I1 && I_final < I2
        is_therapy_control = true;
    else
        is_therapy_control = false;
    end
end

function I2 = compute_upper_bound(params)
    c = params.c;
    dE = params.dE;
    I0E = params.I0E;
    
    I2 = I0E * log(c/(c - dE));
end

function I1 = compute_lower_bound(params, w)
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
            I1 = I_new;
            return;
        end
        
        I_guess = I_new;
    end
    
    I1 = I_guess;
end

function save_spontaneous_therapy_points(spontaneous_points)
    filename = 'spontaneous_control_points.txt';
    
    fid = fopen(filename, 'w');
    if fid == -1
        error('Не удалось открыть файл для записи');
    end
    
    fprintf(fid, 'ТОЧКИ СПОНТАННОГО КОНТРОЛЯ (w=1, без терапии)\n\n');
    fprintf(fid, 'Всего точек: %d\n\n', size(spontaneous_points, 1));
    fprintf(fid, 'tARTstart\ttARTend\t\tI_final\t\t\tH_final\n');
    fprintf(fid, '------------------------------------------------------------\n');
    
    for k = 1:size(spontaneous_points, 1)
        fprintf(fid, '%.2f\t\t%.2f\t\t%.2e\t%.2e\n', ...
            spontaneous_points(k,1), spontaneous_points(k,2), ...
            spontaneous_points(k,3), spontaneous_points(k,4));
    end
    
    fclose(fid);
end
