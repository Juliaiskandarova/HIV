% АНАЛИЗ ЧУВСТВИТЕЛЬНОСТИ: c, dH, dE, dT
clear; clc; close all;


R0 = 8.0;
I0H = 9e-6;
I0E = 0.03;
time = 100.0;
r0 = 0.001;
pL0 = 0.4;
E0L = 0.01;

dT = 0.3;        % Target cell death rate
dI = 1.0;
dE = 0.15;       % Killer cell death rate
dH = 0.15;       % Helper cell death rate

c = 2.0;         % Max. proliferation rate (из диапазона 1.5-3)

rmax = 0.2;
T0 = 1.0;
b = T0 * dT;
p = R0 * dT * dI / b;
H0 = 0.01;
E0 = 0.01;
k = dI / E0;
Em = 1.0;
Iin = 1e-11 * T0;
Lin = Iin;
Hin = 1e-5;
Ein = 1e-5;
effic = 0.97;

% Диапазоны для терапии
tARTstart_vec = 0.5:0.5:10;
tARTend_vec = 4:1:20;
H_threshold = 1e-9;
w_values = [10, 100];
w_spontaneous = 1.0;

fprintf('=== АНАЛИЗ ЧУВСТВИТЕЛЬНОСТИ: c, dH, dE, dT ===\n');
fprintf('Базовые значения:\n');
fprintf('  c  = %.2f (диапазон: 1.5-3)\n', c);
fprintf('  dH = %.2f (диапазон: 0.07-0.3)\n', dH);
fprintf('  dE = %.2f (диапазон: 0.07-0.3)\n', dE);
fprintf('  dT = %.2f (диапазон: 0.1-0.5)\n', dT);
fprintf('\n');

c_values = [2.0, 3];

dH_values = [0.10, 0.30];

dE_values = [0.10, 0.30];

dT_values = [0.05, 0.20];


% ГРАФИК 1: Влияние c 
figure('Position', [100, 100, 1200, 500]);

for c_idx = 1:length(c_values)
    c_current = c_values(c_idx);
    
    for w_idx = 1:length(w_values)
        w_current = w_values(w_idx);
        
        % Параметры
        y0 = [T0, Iin, Lin, Hin, Ein];
        
        params.R0 = R0; params.I0H = I0H; params.I0E = I0E;
        params.time = time; params.r0 = r0; params.pL0 = pL0; params.E0L = E0L;
        params.c = c_current; params.dT = dT; params.dI = dI; 
        params.dE = dE; params.dH = dH;
        params.rmax = rmax; params.b = b; params.p = p;
        params.H0 = H0; params.E0 = E0; params.k = k;
        params.Em = Em; params.Iin = Iin; params.y0 = y0; params.effic = effic;
        
        points = find_control_points(params, w_current, tARTstart_vec, ...
                                      tARTend_vec, w_spontaneous, H_threshold);
        
        subplot(2, 2, (c_idx-1)*2 + w_idx);
        if ~isempty(points)
            scatter(points(:,2), points(:,1), 40, 'filled', ...
                'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');
            title(sprintf('c = %.2f, w = %d\n%d точек', c_current, w_current, ...
                size(points,1)), 'FontSize', 11, 'FontWeight', 'bold');
        else
            text(0.5, 0.5, 'Нет контроля', 'HorizontalAlignment', 'center', ...
                'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
            title(sprintf('c = %.2f, w = %d', c_current, w_current), ...
                'FontSize', 11, 'FontWeight', 'bold');
        end
        
        xlabel('tARTend (days)', 'FontSize', 10);
        ylabel('tARTstart (days)', 'FontSize', 10);
        xlim([min(tARTend_vec)-1, max(tARTend_vec)+1]);
        ylim([min(tARTstart_vec)-0.5, max(tARTstart_vec)+0.5]);
        grid on;
    end
end

sgtitle('Влияние максимальной скорости пролиферации c', ...
    'FontSize', 14, 'FontWeight', 'bold');


% ГРАФИК 2: Влияние dH (helper cell death rate)
figure('Position', [100, 100, 1200, 500]);

for dH_idx = 1:length(dH_values)
    dH_current = dH_values(dH_idx);
    
    for w_idx = 1:length(w_values)
        w_current = w_values(w_idx);
        
        % Пересчитываем зависимые параметры
        b_cur = T0 * dT;
        p_cur = R0 * dT * dI / b_cur;
        k_cur = dI / E0;
        y0 = [T0, Iin, Lin, Hin, Ein];
        
        params.R0 = R0; params.I0H = I0H; params.I0E = I0E;
        params.time = time; params.r0 = r0; params.pL0 = pL0; params.E0L = E0L;
        params.c = c; params.dT = dT; params.dI = dI; 
        params.dE = dE; params.dH = dH_current;
        params.rmax = rmax; params.b = b_cur; params.p = p_cur;
        params.H0 = H0; params.E0 = E0; params.k = k_cur;
        params.Em = Em; params.Iin = Iin; params.y0 = y0; params.effic = effic;
        
        points = find_control_points(params, w_current, tARTstart_vec, ...
                                      tARTend_vec, w_spontaneous, H_threshold);
        
        subplot(2, 2, (dH_idx-1)*2 + w_idx);
        if ~isempty(points)
            scatter(points(:,2), points(:,1), 40, 'filled', ...
                'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');
            title(sprintf('dH = %.2f, w = %d\n%d точек', dH_current, w_current, ...
                size(points,1)), 'FontSize', 11, 'FontWeight', 'bold');
        else
            text(0.5, 0.5, 'Нет контроля', 'HorizontalAlignment', 'center', ...
                'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
            title(sprintf('dH = %.2f, w = %d', dH_current, w_current), ...
                'FontSize', 11, 'FontWeight', 'bold');
        end
        
        xlabel('tARTend (days)', 'FontSize', 10);
        ylabel('tARTstart (days)', 'FontSize', 10);
        xlim([min(tARTend_vec)-1, max(tARTend_vec)+1]);
        ylim([min(tARTstart_vec)-0.5, max(tARTstart_vec)+0.5]);
        grid on;
    end
end

sgtitle('Влияние скорости смерти хелперов dH', ...
    'FontSize', 14, 'FontWeight', 'bold');

% ГРАФИК 3: Влияние dE (killer cell death rate)
figure('Position', [100, 100, 1200, 500]);

for dE_idx = 1:length(dE_values)
    dE_current = dE_values(dE_idx);
    
    for w_idx = 1:length(w_values)
        w_current = w_values(w_idx);
        
        b_cur = T0 * dT;
        p_cur = R0 * dT * dI / b_cur;
        k_cur = dI / E0;
        y0 = [T0, Iin, Lin, Hin, Ein];
        
        params.R0 = R0; params.I0H = I0H; params.I0E = I0E;
        params.time = time; params.r0 = r0; params.pL0 = pL0; params.E0L = E0L;
        params.c = c; params.dT = dT; params.dI = dI; 
        params.dE = dE_current; params.dH = dH;
        params.rmax = rmax; params.b = b_cur; params.p = p_cur;
        params.H0 = H0; params.E0 = E0; params.k = k_cur;
        params.Em = Em; params.Iin = Iin; params.y0 = y0; params.effic = effic;
        
        points = find_control_points(params, w_current, tARTstart_vec, ...
                                      tARTend_vec, w_spontaneous, H_threshold);
        
        subplot(2, 2, (dE_idx-1)*2 + w_idx);
        if ~isempty(points)
            scatter(points(:,2), points(:,1), 40, 'filled', ...
                'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');
            title(sprintf('dE = %.2f, w = %d\n%d точек', dE_current, w_current, ...
                size(points,1)), 'FontSize', 11, 'FontWeight', 'bold');
        else
            text(0.5, 0.5, 'Нет контроля', 'HorizontalAlignment', 'center', ...
                'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
            title(sprintf('dE = %.2f, w = %d', dE_current, w_current), ...
                'FontSize', 11, 'FontWeight', 'bold');
        end
        
        xlabel('tARTend (days)', 'FontSize', 10);
        ylabel('tARTstart (days)', 'FontSize', 10);
        xlim([min(tARTend_vec)-1, max(tARTend_vec)+1]);
        ylim([min(tARTstart_vec)-0.5, max(tARTstart_vec)+0.5]);
        grid on;
    end
end

sgtitle('Влияние скорости смерти киллеров dE', ...
    'FontSize', 14, 'FontWeight', 'bold');


% ГРАФИК 4: Влияние dT (target cell death rate)
figure('Position', [100, 100, 1200, 500]);

for dT_idx = 1:length(dT_values)
    dT_current = dT_values(dT_idx);
    
    for w_idx = 1:length(w_values)
        w_current = w_values(w_idx);
        
        % Пересчитываем все зависимые параметры
        b_cur = T0 * dT_current;
        p_cur = R0 * dT_current * dI / b_cur;  
        k_cur = dI / E0;
        y0 = [T0, Iin, Lin, Hin, Ein];
        
        params.R0 = R0; params.I0H = I0H; params.I0E = I0E;
        params.time = time; params.r0 = r0; params.pL0 = pL0; params.E0L = E0L;
        params.c = c; params.dT = dT_current; params.dI = dI; 
        params.dE = dE; params.dH = dH;
        params.rmax = rmax; params.b = b_cur; params.p = p_cur;
        params.H0 = H0; params.E0 = E0; params.k = k_cur;
        params.Em = Em; params.Iin = Iin; params.y0 = y0; params.effic = effic;
        
        points = find_control_points(params, w_current, tARTstart_vec, ...
                                      tARTend_vec, w_spontaneous, H_threshold);
        
        subplot(2, 2, (dT_idx-1)*2 + w_idx);
        if ~isempty(points)
            scatter(points(:,2), points(:,1), 40, 'filled', ...
                'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');
            title(sprintf('dT = %.2f, w = %d\n%d точек', dT_current, w_current, ...
                size(points,1)), 'FontSize', 11, 'FontWeight', 'bold');
        else
            text(0.5, 0.5, 'Нет контроля', 'HorizontalAlignment', 'center', ...
                'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
            title(sprintf('dT = %.2f, w = %d', dT_current, w_current), ...
                'FontSize', 11, 'FontWeight', 'bold');
        end
        
        xlabel('tARTend (days)', 'FontSize', 10);
        ylabel('tARTstart (days)', 'FontSize', 10);
        xlim([min(tARTend_vec)-1, max(tARTend_vec)+1]);
        ylim([min(tARTstart_vec)-0.5, max(tARTstart_vec)+0.5]);
        grid on;
    end
end

sgtitle('Влияние скорости смерти целевых клеток dT', ...
    'FontSize', 14, 'FontWeight', 'bold');

fprintf('\nАнализ завершен!\n');


function points = find_control_points(params, w, tARTstart_vec, tARTend_vec, ...
                                       w_spontaneous, H_threshold)
    points = [];
    
    for i = 1:length(tARTstart_vec)
        for j = 1:length(tARTend_vec)
            t_start = tARTstart_vec(i);
            t_end = tARTend_vec(j);
            if t_end <= t_start, continue; end
            
            [is_spontaneous, ~, ~] = check_spontaneous_control_extended(...
                t_start, t_end, params, w_spontaneous);
            if is_spontaneous, continue; end
            
            I1_min = compute_lower_bound(params, w);
            I1_max = I1_min * 1.1;
            I2 = compute_upper_bound(params);
            
            is_therapy = check_post_therapy_control_with_bounds(...
                t_start, t_end, params, w, I1_min, I1_max, I2, H_threshold);
            
            if is_therapy
                points = [points; t_start, t_end];
            end
        end
    end
end

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

function is_therapy_control = check_post_therapy_control_with_bounds(tARTstart, tARTend, params, w, I1, I1_max, I2, H_threshold)
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
    
    if c > dE
        I2 = I0E * log(c/(c - dE));
    else
        I2 = I0E * 10;
    end
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