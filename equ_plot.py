"""
Equ_plot.
Графическое отражение функций уравнения и ее производной
"""
import numpy as np
from matplotlib import pyplot as plt

from equ_models import SCALE


# Масштабирование
#   scale - максимальное соотношение максимумов функций на интервале отображения
def scaling(Y, Yd, scale=SCALE):
    if scale is None: return 1
    max_Y = np.max(np.abs(Y))
    max_Yd = np.max(np.abs(Yd)) if isinstance(Yd, np.ndarray) else Yd
    if not max_Y: return 1
    rat = max_Yd / max_Y
    return max(1, rat / scale) if rat > 1 else min(1, rat * scale)


# Построение графика
def show_chart(*function_s, rng=None, txt_fx=None, txt_other=None, no_fill=None):

    rng = (-1, 1, 0.01) if rng is None else rng
    function, *functions_other = (*function_s, None, None)

    plt.figure(figsize=(9, 6))
    if not (txt_fx is None): plt.title(f'График функции f(x): {txt_fx}')
    plt.xlabel('Ось X')
    plt.ylabel('Ось Y')

    X = np.arange(*rng)
    Y = np.array([function(p) for p in X])
    chart = plt.plot(X, Y, '-', color='b', label=txt_fx, linewidth=2)
    x0 = y0 = np.zeros(len(X))
    plt.plot(X, y0, '-', color='k', linewidth=1)
    plt.grid()

    for f_other in functions_other:
        if f_other is None: break
        Y_d = np.array([f_other(p) for p in X])
        Y_d = Y_d / scaling(Y, Y_d)
        plt.plot(X, Y_d, '-', color='r', label=txt_other, linewidth=1)
        plt.plot(x0, Y_d, '-', color='k', linewidth=1)

    plt.plot(x0, Y, '-', color='k', linewidth=1)

    if no_fill is None:
        plt.fill_between(X, Y, where=(Y < 0), color='r', alpha=0.5)
        plt.fill_between(X, Y, where=(Y > 0), color='g', alpha=0.5)

    plt.legend()
    plt.show()

    return chart
