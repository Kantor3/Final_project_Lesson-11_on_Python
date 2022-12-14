"""
Equ_controller.
Управление процессами расчета и вывода результатов в виде текстовой и графической информации
"""
from sympy import *

from equ_models import EPSILON
from equ_models import f
import equ_models as m
from equ_plot import show_chart
from equ_solution import get_comfortable_boundaries, round_eps, create_intervals, roots_interval, roots_sep
from equ_view import show_info


# 1. Определение корней уравнения
def get_roots(function, f_info, eps=EPSILON):
    roots_self = roots_interval(function, eps, f_info)
    if roots_self is None: return None                                   # Не известная ошибка
    return roots_self


# 2. Найти интервалы, на которых функция возрастает
# 3. Найти интервалы, на которых функция убывает.
# Интервалы, на которых функция монотонна (только возрастает или только убывает),
# соответствуют знакам ее производной: производная > 0 - функция возрастает; если < 0 - убывает.
# При этом, эти знаки чередуются. Поэтому достаточно определить знак значения функции
# производной слева от 1-го корня. Промежутки производной > 0 / < 0 лежат между корнями
# ее уравнения, на уровнях этих корней находятся также экстремумы самой функции.
def get_intervals_id(function, eps=EPSILON):

    save_Function = m.Function_us_def                                     # сохраняем исходную функцию

    # Корни производной
    x = Symbol('x')
    m.Function_us_def = str(eval(m.Function_us_def).diff(x, 1))             # Подменяем исходную функцию ее производной
    f_txt, f_info = f(None, what_ret=0)
    if f_txt is None: show_info(code=0, exit_prog=True)
    roots_self = get_roots(function, f_info)
    roots_d, intervals_show = create_intervals(roots_self, 0, eps=eps)  # Интервалы возрастания - убывания

    m.Function_us_def = save_Function                                     # восстанавливаем исходную функцию
    return intervals_show, roots_d


# 5. Вычислить вершину
# Если вершиной считать самый максимальный экстремум > 0, то он располагается в одной из точек, где производная = 0
def get_peak(roots_l):
    return round_eps(max(map(f, roots_l)))


# 6. Определить промежутки, на которых f > 0
# 7. Определить промежутки, на которых f < 0.
# Промежутки функции [ > 0 / < 0 ] лежат между корнями ее уравнения, при этом, эти промежутки чередуются.
# Поэтому достаточно определить знак значения этой функции слева от 1-го корня
def get_intervals_pn(roots_l, eps=EPSILON):
    return create_intervals(roots_l, 1, eps=eps)


# Мастер метод
def main():

    # Дана алгебраическая функция: '-12 * x**4 - 18 * x**3 + 5 * x**2 + 10 * x + 0'
    # Функция задана строкой в синтаксисе Phyton. В дальнейшем может рассматриваться как параметр,
    # задаваемый пользователем (вводом с консоли) или передаваемый при вызове метода определения корней
    # Условие корректной работы программы является указание всех членов многочлена, без пропуска.
    # Если к-л член отсутствует, его коэффициент устанавливается = 0
    m.Function_us_def = '-12 * x**4 - 18 * x**3 + 5 * x**2 + 10 * x + 0'

    txt_f, kf = f(None, what_ret=0)
    if txt_f is None: show_info(code=0, exit_prog=True)

    # 1. Определить корни уравнения
    roots_1 = get_roots(f, kf)
    if roots_1 is None: show_info(code=9, exit_prog=True)
    show_info(code=11, info=roots_1)

    # 2. Найти интервалы, на которых функция возрастает
    # 3. Найти интервалы, на которых функция убывает
    intervals_id, roots_df = get_intervals_id(f)
    show_info(code=12, info=intervals_id, add_info=True)

    # 4. Построить график
    x_arg = Symbol('x')
    f_diff = lambda p: f(None, what_ret=1).subs(x_arg, p)

    # 5. Вычислить вершину
    peak = get_peak(roots_df)
    show_info(code=15, info=peak, add_info=True)

    # 6. Определить промежутки, на которых f > 0
    # 7. Определить промежутки, на которых f < 0
    _, intervals_pn = get_intervals_pn(roots_1)
    show_info(code=16, info=intervals_pn, add_info=True)

    rang = get_comfortable_boundaries(peak, kf)
    show_chart(f, f_diff, rng=rang, txt_fx=txt_f, txt_other='Производная')   # график функции
