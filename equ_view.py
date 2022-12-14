"""
Equ_view.
Отражение результатов расчета и анализа
"""
from functools import reduce
from equ_models import f, Function_us_def


def print_interval(info, point, txt_info, sgn):
    print(f'\n{point}. {txt_info} {info[sgn][0]}')
    print(*info[sgn][1], sep='\n')


# Вывод результаты анализа функции, в т.ч. корни уравнения
# info_code - что выводим
#        11 - Вывод информацию о результатах расчета по 1-му пункту
#        12 - Вывод информацию о результатах расчета по 2 и 3-му пунктам
#        15 - Вывод информацию о результатах расчета по 5-му пункту
#        16 - Вывод информацию о результатах расчета по 6 и 7-му пунктам
def show_result(info, info_code):

    # 1. Корни уравнения
    if info_code == 11:
        print('1. Корни:')
        if info:
            print(*[f'{i+1}-й   {ri[0]}' for i, ri in enumerate(info)], sep='\n')
            print(f'Найдены за {reduce(lambda i, ri: i + ri[1], info, 0)} итераций')
        else:
            print('Нет корней')
        return

    # 2. Интервалы, на которых функция возрастает
    # 3. Интервалы, на которых функция убывает
    # {1: ('возрастает', (0, 2, 4)), -1: ('убывает', (1, 3, 5))}
    if info_code == 12:
        txt = 'Интервалы, на которых функция'
        print_interval(info, 2, txt, 1)
        print_interval(info, 3, txt, -1)
        return

    # 5. Вершина функции
    if info_code == 15:
        print('\n5. Вершина функции:')
        print(info, sep='\n')
        return

    # 6. Определить промежутки, на которых f > 0
    # 7. Определить промежутки, на которых f < 0
    if info_code == 16:
        txt = 'Промежутки, на которых'
        print_interval(info, 6, txt, 1)
        print_interval(info, 7, txt, -1)
        return


# Вывод информации, в т.ч. по завершении программы и собственно ее завершение
# code - код информации, если информация носит общий характер (ошибка, предупреждение и пр.):
#        None (по умолчанию) - показать результаты расчета
#        0 - Поддерживается поиск корней только алгебраических уравнений!
#        1 - Ошибка в задании выражения уравнения
#        9 - не известная ошибка и пр.
#        11...16 - Вывод информацию о результатах расчета по отдельным пунктам
#        12 - Вывод информацию о результатах расчета по 2 и 3-му пунктам
#        15 - Вывод информацию о результатах расчета по 5-му пункту
#        16 - Вывод информацию о результатах расчета по 6 и 7-му пунктам
# exit_prog - True - для немедленного завершения работы программы
# info - информация для вывода, например, результаты расчета
def show_info(code=None, info=None, exit_prog=None, add_info=None):
    if not add_info:
        title_info = f'\nАнализ уравнения "{f(None, what_ret=0)[0]}":'
        print(title_info)
        print('-' * len(title_info))
    if code == 0:
        print('\nПоддерживается поиск корней только алгебраических уравнений!')
    elif code == 1:
        print(f'\nОшибка в выражении уравнения {Function_us_def}. Исправьте выражение и запустите расчет вновь')
    elif code == 9:
        print('\nНе известная ошибка. Обратитесь к разработчику')
    elif code > 10:
        show_result(info, code)
    if exit_prog:
        print('\nРабота с программой завершена. До встречи!')
        exit()
