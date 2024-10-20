import sys
import logging
from typing import List, Tuple, Set
import re
from polynomial import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f'log/{__name__}.log')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)

def solve_comparisons (equations: List[Tuple[Polynomial, int]]) -> Set[Tuple[int,int]]:
    logger.debug(f'solve_comparisons. equations={equations}')
    result = equations[0][0].roots(equations[0][1])
    for poly, module in equations[1:]:
        logger.debug(f'solve_comparisons. equations={equations}, result={result}')
        new_result = set()
        for a,m in result:
            new_poly = Polynomial()
            for exp, coeff in poly.monomials():
                new_poly += coeff * binomial_theorem(m,a,exp)
                logger.debug(f'solve_comparisons. result={result}, new_poly={new_poly}')
            for a1,m1 in new_poly.roots(module):
                new_result.add((m*a1+a,m*m1))
        result = new_result
    logger.debug(f'solve_comparisons. equations={equations}, result={result}')
    return result    

def parse_equation(equation: str) -> Tuple[Polynomial, int]:
    logger.debug(f'parse_equation. equation={equation}')
    pattern = r'(.+?)\s*=\s*(.+?)\s*\(\s*mod\s*(\d+)\s*\)'
    match = re.match(pattern, equation)
    if not match:
        logger.critical(f'parse_equation. equation={equation}. Неверный формат уравнения.')
        raise ValueError("Неверный формат строки. Ожидается формат 'P(x) = Q(x) (mod m)'")
    module = int(match.group(3).strip())
    poly = Polynomial(match.group(1).strip()) - Polynomial(match.group(2).strip())
    logger.debug(f'parse_equation. equation={equation}, poly={poly}, module={module}')
    return poly, module

def main():
    logger.info("Начало программы.")
    if len(sys.argv) < 3:
        logger.error(f'Неправильное количество аргументов.')
        print("Запускать: python3", sys.argv[0], "<файл_с_сравнениями> <файл для вывода результата>")
        return
    logger.info(f'Считывание с файла {sys.argv[1]}.')
    with open(sys.argv[1]) as fp: equations_lines = fp.readlines()
    equations = []
    for equation in equations_lines:
        logger.info(f'Парсинг {equation}')
        if '#' in equation: continue
        equations.append(parse_equation(equation))
    if not equations:
        logger.info(f'Нет незакомментированных уравнений в {sys.argv[1]}.')
        print(f'Нет незакомментированных уравнений в {sys.argv[1]}.')
        return
    logger.info("Начинает решать.")
    result = solve_comparisons(equations)
    logger.info("Закончила решать. Начинает записывать ответ.")
    with open(sys.argv[2],"w") as fp:
        for x,m in result:
            print(f'[{x}]_{m}',end=' ', file=fp)
    logger.info("Записала ответ.")

if __name__ == '__main__':
    main()
    handler.close()