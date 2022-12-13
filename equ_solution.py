"""
Equ_solution.
Процедуры подбора подходящих методов и решения с их использованием уравнений
"""
import math
from sympy import *

from equ_models import EPSILON, INF_cod, SCALE, STEP, SIGNS_named
from equ_models import f
import equ_models as m


# Округление результатов в соответствии с точностью расчета
def round_eps(data, eps=EPSILON):
    number_rounding = int(abs(math.log10(eps)) + 1)
    return round(data, number_rounding)                                         # корни округляем до точности eps+1


# Формирование списка интервалов между корнями по списку корней
def get_intervals(roots_l):
    roots_l.insert(0, f'-{chr(INF_cod)}')                                             # дополняем бесконечностями
    roots_l.append(chr(INF_cod))
    return [f'({rl} ... {rr})' for rl, rr in zip(roots_l[:-1], roots_l[1:])]    # формируем список интервалов


# Определение коэффициентов смежных алгебраических уравнений:
# Базовое уравнение Pn(x) = 0   Пример:     (-12, -18, 5, 10, 0) - базовые коэффициенты. Смежные уравнения:
# numb
#   0  P0(x) = Pn(x) = 0           Возвращаем: (12, 18, -5, -10, 0) - старший коэффициент должен быть положительным
#   1  P1(x) = x^n * Pn(1/x) = 0   Возвращаем: (10, 5, -18, -12)    - развернуть список, ведущие нули убрать
#   2  P2(x) = Pn(-x) = 0          Возвращаем: (12, -18, -5, 10, 0) - меняем знаки у коэффициентов нечетных степеней
#   3  P3(x) = x^n * Pn(-1/x) = 0  Возвращаем: (10, -5, -18, 12)    - разворачиваем список и меняем знаки
# Алгоритм:
# 1) для нечетных вариантов развернем список, ведущие нули уберем (степень полинома может уменьшиться)
# 2) для 3, 4 вариантов меняем знаки коэффициентов с нечетными индексами если степень уравнения четная и наоборот
# 3) старший коэффициент должен быть положительным, поэтому меняем знак на противоположный, если иначе
def get_kf_related(kfs, numb):
    if not kfs: return kfs
    ks = tuple(filter(lambda el: el, [kfs[-(1 + i)] for i in range(len(kfs))])) if numb % 2 else kfs
    n = len(ks) - 1
    ks = tuple([el if i % 2 == n % 2 else -el for i, el in enumerate(ks)]) if numb > 1 else ks
    ks = tuple(map(lambda el: -el, ks)) if ks[0] < 0 else ks
    return ks


# Расчет верхней границы положительных корней смежных уравнений - R
# В соответствии с теоремой Лагранжа о верхней границе положительных корней
#   R = 1 + (C/an) ^ 1/(n-i)
def get_upb_roots(kfs):
    C = max(map(lambda el: -el, kfs))
    n = len(kfs) - 1
    i = n - [i for i, el in enumerate(kfs) if el < 0][0]
    an = kfs[0]
    return 1 + (C/an) ** (1/(n-i))


# Отделение корней
# В соответствии с теоремой Лагранжа положительные корни +X*i и отрицательные корни -x*i
# уравнения удовлетворяют неравенствам 1/R1 <= +x*i <= R, -R2 <= -x*i <= -1/R3, соответственно,
#   где R, R1, R2, R3 = верхней границы положительных корней смежных уравнений (Rs - кортеж R, R1,...Ri)
#       br - граница корней = (0, 1, 2, 3) - (левая/правая граница положительных корней, л/п граница отрицательных)
def roots_sep(f_info, br=1):
    if not isinstance(f_info, tuple): return None               # другие виды уравнений не реализованы
    l_kfs   = [get_kf_related(f_info, i) for i in range(4)]     # списки коэффициентов смежных уравнений
    Rs      = tuple([get_upb_roots(el) for el in l_kfs])
    R, R1, R2, R3 = Rs
    roots_borders = [1/R1, R, -R2, -1/R3]
    return roots_borders[br]


# Поиск корня в интервале (..., b) - Метод секущих
def secant_method(function, eps, b):
    delta = eps                             # малое положительное число для расчета первой производной
    x0 = b                                  # Начальное приближение
    f0_div = (function(x0) - function(x0-delta)) / delta  # - крайне левое значение из интервала нахождения корней
    x1 = x0 - function(x0) / f0_div                # Первое приближение

    # Остальные приближения
    delta = x1 - x0
    xi_1 = x0
    xi = x1
    yi_1 = f(x0)
    i = 1                                   # число итераций
    while abs(delta) > eps:
        xi_2 = xi_1
        xi_1 = xi
        yi_2 = yi_1
        yi_1 = function(xi_1)
        delta = xi_1 - xi_2
        if yi_1 == yi_2:
            xi = None                       # больше корней нет
            break
        xi = xi_1 - yi_1 / (yi_1 - yi_2) * delta
        i += 1

    return xi, i


# Поиск всех корней
def roots_interval(function, eps, bi):
    save_Function = m.Function_us_def                           # сохраняем исходную функцию
    roots_self = []
    while True:
        bi -= eps
        root_c = secant_method(function, eps, bi)
        bi, i = root_c
        if bi is None:
            m.Function_us_def = save_Function                   # восстанавливаем исходную функцию
            break                                               # больше корней нет
        root_c = round_eps(bi), i                               # корни округляем до точности eps + 1
        roots_self.append(root_c)                               # корни округляем до точности eps + 1
        m.Function_us_def = f'({m.Function_us_def})/(x - {bi})'     # Отсекаем поиск уже найденных корней
    return roots_self


# Создание и группировка интервалов по оси X из его граничных точек - корней уравнения
def create_intervals(roots_ri, n_name_gr, eps=EPSILON):
    # Знак значения функции в точке левее первого корня ее уравнения
    # определяет чередующиеся интервалы + / возрастания и - / убывания
    roots_r = [r for r, i in roots_ri]
    roots_r.sort()
    root_min = roots_r[0]
    sgn = int(sign(f(root_min - eps)))

    # Список интервалов, располагающихся между корнями производной
    l_intervals = get_intervals(roots_r)

    # Формируем списки интервалов + / возрастания и - / убывания анализируемой функции
    intervals_show = dict()
    intervals_show[1] = intervals_show[-1] = None
    intervals_show[sgn]  = SIGNS_named[sgn][n_name_gr], l_intervals[::2]
    intervals_show[-sgn] = SIGNS_named[-sgn][n_name_gr], l_intervals[1::2]

    return roots_r[1:-1], intervals_show


# Определение границ графика для его комфортной визуализации
def get_comfortable_boundaries(y_top):

    save_Function = m.Function_us_def   # сохраняем исходную функцию
    _, kfs = f(None, what_ret=0)        # получим коэффициенты исходного уравнения

    # Расчет максимальной амплитуды графика для комфортной визуализации графика
    Y_max = abs(y_top * SCALE)
    Yl = f(roots_sep(kfs, br=2))
    Yr = f(roots_sep(kfs, br=1))
    Y_max_l = min(Y_max, abs(Yl)) * sign(Yl)
    Y_max_r = min(Y_max, abs(Yr)) * sign(Yr)

    x = Symbol('x')
    m.Function_us_def = str(eval(f'{save_Function} - {Y_max_l}'))   # составляем уравнение для нахождения левой границы
    Xl_lim, _ = secant_method(f, EPSILON, roots_sep(kfs, br=2))     # решение уравнения => левая граница графика
    m.Function_us_def = str(eval(f'{save_Function} - {Y_max_r}'))   # составляем уравнение для нахождения левой границы
    Xr_lim, _ = secant_method(f, EPSILON, roots_sep(kfs, br=1))     # решение уравнения => правая граница графика

    m.Function_us_def = save_Function
    return Xl_lim, Xr_lim, STEP
