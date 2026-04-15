% ПАРАМЕТРЫ МОДЕЛИ
R0 = 8.0;
I0E = 0.03;
I0H = 9e-6;
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

% Создаем структуру параметров
params.R0 = R0;
params.I0H = I0H;
params.I0E = I0E;
params.time = time;
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



% Значение w для спонтанного контроля
w_spontaneous = 1.0;
w_values = [10, 50, 100, 1000, 10000]; % значения для посттерапевтического контроля

% Диапазоны для tARTstart и tARTend
tARTstart_vec = 2.0:0.17:10;
tARTend_vec = 6:0.25:20;

% Порог для хелперов 
H_threshold = 1e-9;

spontaneous_points = []; % точки спонтанного контроля
therapy_points = cell(length(w_values), 1); % точки контроля для каждого w

for w_idx = 1:length(w_values)
    therapy_points{w_idx} = [];
end


fprintf('Анализ точек...\n');

for i = 1:length(tARTstart_vec)
    for j = 1:length(tARTend_vec)
        tARTstart_val = tARTstart_vec(i);
        tARTend_val = tARTend_vec(j);
        
        if tARTend_val <= tARTstart_val
            continue; 
        end
        
    
        % ШАГ 1: Проверка на спонтанный контроль (w=1, без терапии)
     
        [is_spontaneous, I_final, H_final] = check_spontaneous_control_extended(...
            tARTstart_val, tARTend_val, params, w_spontaneous);
        
        if is_spontaneous
            % Сохраняем точку спонтанного контроля
            spontaneous_points = [spontaneous_points; tARTstart_val, tARTend_val, I_final, H_final];
            
            % Переходим к следующей точке
            continue; 
        end
        
 
        % ШАГ 2: Анализ с терапией для разных w
        for w_idx = 1:length(w_values)
            w = w_values(w_idx);
            
            % Вычисляем границы из двух уравнений
            [I1_min, I1_max] = calculate_boundary_from_two_eqs(params.c, params.dH, params.p, w, params.I0H);
            I2 = compute_upper_bound(params);
            
            % Проверяем посттерапевтический контроль
            is_therapy_control = check_post_therapy_control_with_bounds(...
                tARTstart_val, tARTend_val, params, w, I1_min, I1_max, I2, H_threshold);
            
            if is_therapy_control
                therapy_points{w_idx} = [therapy_points{w_idx}; tARTstart_val, tARTend_val];
            end
        end
    end
end

fprintf('Анализ завершен.\n');


% АНАЛИЗ ШИРИНЫ ПОЛОСЫ
fprintf('\n=== АНАЛИЗ ШИРИНЫ ПОЛОСЫ ===\n');

% Анализируем ширину полосы для каждого w
bandwidth_results = cell(length(w_values), 1);

for w_idx = 1:length(w_values)
    w = w_values(w_idx);
    points = therapy_points{w_idx};
    
    fprintf('w = %d: всего точек = %d\n', w, size(points, 1));
    
    % Инициализация
    bandwidth = 0;
    start_min = 0;
    start_max = 0;
    is_detected = false;
    
    if ~isempty(points) && size(points, 1) >= 5
        % Простой анализ ширины полосы
        start_min = min(points(:,1));
        start_max = max(points(:,1));
        bandwidth = start_max - start_min;
        
        % Проверяем, достаточно ли точек для образования полосы
        unique_starts = unique(points(:,1));
        
        if bandwidth > 0.1 && length(unique_starts) >= 3
            % Проверяем плотность покрытия
            expected_count = round((start_max - start_min) / 0.17) + 1; 
            actual_count = length(unique_starts);
            density = actual_count / expected_count;
            
            if density > 0.3  % Если покрыто более 30%
                is_detected = true;
                fprintf('  Обнаружена полоса: d = %.2f дней (от %.2f до %.2f), плотность = %.1f%%\n', ...
                    bandwidth, start_min, start_max, density*100);
            else
                fprintf('  Плотность слишком низкая: %.1f%%\n', density*100);
            end
        else
            fprintf('  Слишком узкая полоса или мало уникальных значений\n');
        end
    else
        fprintf('  Недостаточно точек для анализа\n');
    end
    
    % Сохраняем результаты
    bandwidth_results{w_idx}.bandwidth = bandwidth;
    bandwidth_results{w_idx}.start_min = start_min;
    bandwidth_results{w_idx}.start_max = start_max;
    bandwidth_results{w_idx}.is_detected = is_detected;
end


% ВЫВОД 
fprintf('\n=== РЕЗУЛЬТАТЫ АНАЛИЗА ===\n');

% Точки спонтанного контроля
if ~isempty(spontaneous_points)
    % Сохраняем в файл
    save_spontaneous_therapy_points(spontaneous_points);
    fprintf('\nТочки спонтанного контроля сохранены в файл: spontaneous_control_points.txt\n');
else
    fprintf('\nНет точек спонтанного контроля\n');
end

% Сводка по ширине полос
fprintf('\n=== СВОДКА ПО ШИРИНЕ ПОЛОС ===\n');
for w_idx = 1:length(w_values)
    w = w_values(w_idx);
    result = bandwidth_results{w_idx};
    
    if result.is_detected
        fprintf('w = %d - d = %.2f дней (от %.2f до %.2f)\n', ...
            w, result.bandwidth, result.start_min, result.start_max);
    else
        fprintf('w = %d - четкая полоса не обнаружена\n', w);
    end
end


% ПОСТРОЕНИЕ ГРАФИКОВ
% Графики для каждого w
for w_idx = 1:length(w_values)
    w = w_values(w_idx);
    
    figure;
    hold on;
    grid on;
    title(sprintf('Посттерапевтический контроль\nw=%d R0 = %d I0E = %.1e IOH = %.1e', w, R0, I0E, I0H), 'FontSize', 14);
    xlabel('tARTend (days)', 'FontSize', 12);
    ylabel('tARTstart (days)', 'FontSize', 12);
    
    % Отображаем точки контроля
    if ~isempty(therapy_points{w_idx})
        plot(therapy_points{w_idx}(:,2), therapy_points{w_idx}(:,1), ...
            'r+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Контроль (+)');
        
        % Добавляем легенду 
        legend('Location', 'best');
        
        % Если обнаружена полоса, рисуем ее границы
        result = bandwidth_results{w_idx};
        if result.is_detected
            % Рисуем горизонтальные линии, ограничивающие полосу
            plot([min(tARTend_vec), max(tARTend_vec)], ...
                 [result.start_min, result.start_min], ...
                 'b--', 'LineWidth', 1, 'HandleVisibility', 'off');
            plot([min(tARTend_vec), max(tARTend_vec)], ...
                 [result.start_max, result.start_max], ...
                 'b--', 'LineWidth', 1, 'HandleVisibility', 'off');
            
            % Добавляем текст с шириной полосы
            text(mean(tARTend_vec), result.start_min - 0.2, ...
                sprintf('d=%.2f', result.bandwidth), ...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'blue');
        end
    end
    
    % Настраиваем график
    xlim([min(tARTend_vec)-0.5, max(tARTend_vec)+0.5]);
    ylim([min(tARTstart_vec)-0.25, max(tARTstart_vec)+0.25]);
    set(gca, 'FontSize', 10);
    
    hold off;
end


% ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ

function dxdt = odefun_no_therapy(t, x, R0, I0H, I0E, w, r0, pL0, E0L, c, dT, dI, dE, dH, rmax, b, p, H0, E0, k, Em, Iin)
    T = x(1); I = x(2); L = x(3); H = x(4); E = x(5);
    
    alphaH = 1 - exp(-I / I0H);
    alphaH = alphaH * max(0, (1 - H / Em))^(1/3);
    sigma = 1 - exp(-alphaH * H / H0 - I / I0E);
    sigma = sigma * max(0, (1 - E / Em))^(1/3);
    
    r = r0 + rmax * E / (E0L + E);
    pL = pL0 * E0L / (E0L + E);
    
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
    dIdt = ((1 - pL) * p1 * I * T - k * I * E - dI * I + r * L) * (I / (I + 1e-11));
    dLdt = pL * p1 * I * T - r * L;
    dHdt = (c * alphaH * H - p1 * w * I * H - dH * H) * (H / (H + 1e-11));
    dEdt = (c * sigma * E - dE * E) * (E / (E + 1e-11));
    
    dxdt = [dTdt; dIdt; dLdt; dHdt; dEdt];
end

function dxdt = odefun_with_therapy(t, x, R0, I0H, I0E, w, tARTstart, tARTend, effic, r0, pL0, E0L, c, dT, dI, dE, dH, rmax, b, p, H0, E0, k, Em, Iin)
    T = x(1); I = x(2); L = x(3); H = x(4); E = x(5);
    
    alphaH = 1 - exp(-I / I0H);
    alphaH = alphaH * max(0, (1 - H / Em))^(1/3);
    sigma = 1 - exp(-alphaH * H / H0 - I / I0E);
    sigma = sigma * max(0, (1 - E / Em))^(1/3);
    
    r = r0 + rmax * E / (E0L + E);
    pL = pL0 * E0L / (E0L + E);
    
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
    dIdt = ((1 - pL) * p1 * I * T - k * I * E - dI * I + r * L) * (I / (I + 1e-11));
    dLdt = pL * p1 * I * T - r * L;
    dHdt = (c * alphaH * H - p1 * w * I * H - dH * H) * (H / (H + 1e-11));
    dEdt = (c * sigma * E - dE * E) * (E / (E + 1e-11));
    
    dxdt = [dTdt; dIdt; dLdt; dHdt; dEdt];
end

function [is_spontaneous, I_final, H_final] = check_spontaneous_control_extended(tARTstart, tARTend, params, w)
    % ШАГ 1: Проверка спонтанного контроля (без терапии)
    % Моделируем до tARTstart без терапии
    
    R0 = params.R0; I0H = params.I0H; I0E = params.I0E;
    r0 = params.r0; pL0 = params.pL0; E0L = params.E0L;
    c = params.c; dT = params.dT; dI = params.dI; dE = params.dE; dH = params.dH;
    rmax = params.rmax; b = params.b; p = params.p; H0 = params.H0; E0 = params.E0;
    k = params.k; Em = params.Em; Iin = params.Iin;
    
    % Моделируем до tARTstart БЕЗ терапии
    tspan = [0, tARTstart];
    y0 = params.y0;
    
    % Решаем ОДУ БЕЗ терапии
    [t, x] = ode15s(@(t,x) odefun_no_therapy(t, x, R0, I0H, I0E, w, r0, pL0, E0L, c, dT, dI, dE, dH, rmax, b, p, H0, E0, k, Em, Iin), ...
                    tspan, y0, odeset('RelTol', 1e-3, 'AbsTol', 1e-14));
    
    I_final = x(end, 2);
    H_final = x(end, 4);
    
    % Критерий спонтанного контроля
    if H_final > 10 * I_final && I_final < 1e-3
        is_spontaneous = true;
    else
        is_spontaneous = false;
    end
end


function is_therapy_control = check_post_therapy_control_with_bounds(tARTstart, tARTend, params, w, I1_min, I1_max, I2, H_threshold)
    % ШАГ 2: Проверка посттерапевтического контроля (с терапией)
    
    R0 = params.R0; I0H = params.I0H; I0E = params.I0E;
    r0 = params.r0; pL0 = params.pL0; E0L = params.E0L;
    c = params.c; dT = params.dT; dI = params.dI; dE = params.dE; dH = params.dH;
    rmax = params.rmax; b = params.b; p = params.p; H0 = params.H0; E0 = params.E0;
    k = params.k; Em = params.Em; Iin = params.Iin;
    effic = params.effic; time = params.time;
    
    % Моделируем полный цикл С терапией
    tspan = [0, time];
    y0 = params.y0;
    
    % Решаем ОДУ С терапией
    [t, x] = ode15s(@(t,x) odefun_with_therapy(t, x, R0, I0H, I0E, w, tARTstart, tARTend, effic, r0, pL0, E0L, c, dT, dI, dE, dH, rmax, b, p, H0, E0, k, Em, Iin), ...
                    tspan, y0, odeset('RelTol', 1e-3, 'AbsTol', 1e-14));
    
    I_final = x(end, 2);
    H_final = x(end, 4);
    
 
    if H_final > H_threshold && I_final > I1_min && I_final < I1_max && I_final < I2
        is_therapy_control = true;
    else
        is_therapy_control = false;
    end
end

function I2 = compute_upper_bound(params)
    % Верхняя граница I2
    c = params.c;
    dE = params.dE;
    I0E = params.I0E;
    
    I2 = I0E * log(c/(c - dE));
end


function [I1_min, I1_max] = calculate_boundary_from_two_eqs(c, dH, p, w, I0H)
    % Расчет границ из двух уравнений 
    % Уравнение 1: I = I0H * log(c / (c - dH - p * w * I))
    % Уравнение 2: I = [c * (1 - exp(-I/I0H)) - dH] / (p * w)
    
    max_iter = 1000;
    tolerance = 1e-12;
    
    % Решение уравнения 1
    solutions_eq1 = [];
    initial_guesses = [I0H/1000, I0H/100, I0H/10, I0H, I0H*10, I0H*100];
    
    for init_guess = initial_guesses
        I_prev = init_guess;
        
        for iter = 1:max_iter
            denominator = c - dH - p * w * I_prev;
            
            if denominator <= 0 || denominator >= c
                break;
            end
            
            I_new = I0H * log(c / denominator);
            
            if abs(I_new - I_prev) < tolerance
                if I_new > 0 && I_new < 1e6
                    is_new = true;
                    for sol = solutions_eq1
                        if abs(I_new - sol) < tolerance * 10
                            is_new = false;
                            break;
                        end
                    end
                    if is_new
                        solutions_eq1 = [solutions_eq1, I_new];
                    end
                end
                break;
            end
            
            I_prev = I_new;
        end
    end
    
    % Решение уравнения 2
    solutions_eq2 = [];
    
    for init_guess = initial_guesses
        I_prev = init_guess;
        
        for iter = 1:max_iter
            numerator = c * (1 - exp(-I_prev/I0H)) - dH;
            
            if numerator <= 0
                I_new = 0;
            else
                I_new = numerator / (p * w);
            end
            
            if abs(I_new - I_prev) < tolerance
                if I_new > 0 && I_new < 1e6
                    is_new = true;
                    for sol = solutions_eq2
                        if abs(I_new - sol) < tolerance * 10
                            is_new = false;
                            break;
                        end
                    end
                    if is_new
                        solutions_eq2 = [solutions_eq2, I_new];
                    end
                end
                break;
            end
            
            I_prev = I_new;
        end
    end
    
    % Объединение решений
    all_solutions = unique([solutions_eq1, solutions_eq2]);
    
    if length(all_solutions) >= 2
        all_solutions = sort(all_solutions);
        I1_min = all_solutions(1);
        I1_max = all_solutions(end);
    else
        if ~isempty(all_solutions)
            I1_min = min(all_solutions);
            I1_max = max(all_solutions);
        else
            % Резервный расчет
            I1_min = calculate_simple_boundary(c, dH, p, w, I0H, true);
            I1_max = calculate_simple_boundary(c, dH, p, w, I0H, false);
        end
    end
end

function I_boundary = calculate_simple_boundary(c, dH, p, w, I0H, is_min)
    % Простой расчет границы
    max_iter = 100;
    tolerance = 1e-10;
    
    if is_min
        I_prev = I0H * 0.1;
    else
        I_prev = I0H * 10;
    end
    
    for iter = 1:max_iter
        numerator = c * (1 - exp(-I_prev/I0H)) - dH;
        
        if numerator <= 0
            I_new = I_prev;
        else
            I_new = numerator / (p * w);
        end
        
        if abs(I_new - I_prev) < tolerance
            I_boundary = I_new;
            return;
        end
        
        I_prev = I_new;
    end
    
    I_boundary = I_prev;
end

function save_spontaneous_therapy_points(spontaneous_points)
    % Сохранение точек спонтанного контроля
    filename = 'spontaneous_control_points.txt';
    
    fid = fopen(filename, 'w');
    if fid == -1
        error('Не удалось открыть файл для записи');
    end
    
    fprintf(fid, 'ТОЧКИ СПОНТАННОГО КОНТРОЛЯ (w=1, без терапии)\n\n');
    fprintf(fid, 'Всего точек: %d\n\n', size(spontaneous_points, 1));
    fprintf(fid, 'tARTstart\ttARTend\n');
    fprintf(fid, '----------------------\n');
    
    for k = 1:size(spontaneous_points, 1)
        fprintf(fid, '%.1f\t\t%.1f\n', ...
            spontaneous_points(k,1), spontaneous_points(k,2));
    end
    
    fclose(fid);
    fprintf('Точки спонтанного контроля сохранены в: %s\n', filename);
end
