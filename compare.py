"""Основной модуль."""

import sys
import logging
from comparisons import solve_comparisons, parse_equation


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f'log/{__name__}.log')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)


def main():
    """Основная функция программы."""

    logger.info("Начало программы.")
    if len(sys.argv) < 3:
        logger.error('Неправильное количество аргументов.')
        print("Запускать: python3",
              sys.argv[0],
              "<файл_с_сравнениями> <файл для вывода результата>")
        return
    logger.info('Считывание с файла %s.', sys.argv[1])
    with open(sys.argv[1], encoding='utf-8') as fp:
        equations_lines = fp.readlines()
    equations = []
    for equation in equations_lines:
        logger.info('Парсинг %s', equation)
        if '#' in equation:
            continue
        equations.append(parse_equation(equation))
    if not equations:
        logger.info('Нет не закомментированных уравнений в %s.', sys.argv[1])
        print(f'Нет не закомментированных уравнений в {sys.argv[1]}.')
        return
    logger.info("Начинает решать.")
    result = solve_comparisons(equations)
    logger.info("Закончила решать. Начинает записывать ответ.")
    with open(sys.argv[2],"w", encoding='utf-8') as fp:
        for x in result:
            print(x, end=' ', file=fp)
    logger.info("Записала ответ.")


if __name__ == '__main__':
    main()
