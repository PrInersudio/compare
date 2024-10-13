import sys
import re
from math import gcd

def solve_linear_comparison(equation):
    a,b,m = equation
    d = gcd(a,m)
    if b % d != 0: return []
    a0,b0,m0 = a//d,b//d,m//d
    x0 = b0 * pow(a0,-1,m0)
    return [((x0 + i) % m, m) for i in range(0,m,m0)]

def solve_linear_comparisons (equations):
    result = solve_linear_comparison(equations[0])
    for a,b,m in equations[1:]:
        new_equations = [(a*m0,b-a0*a,m) for a0,m0 in result]
        temp_results = [solve_linear_comparison(equations) for equations in new_equations]
        new_result = []
        for (a0,m0),temp_result in zip(result,temp_results):
            for a1,m1 in temp_result:
                new_result.append((m0*a1+a0,m0*m1))
    result = new_result
    return result

if len(sys.argv) < 3:
    print("Запускать: python3", sys.argv[0], "<файл_с_сравнениями> <файл для вывода результата>")
    exit()

with open(sys.argv[1]) as fp:
    equations_lines = fp.readlines()

equation_regex = re.compile(r'(\d+)\s*\*?\s*x\s*=\s*(\d+)\s*\(\s*mod\s*(\d+)\s*\)')
equations = []
for equation in equations_lines:
    mo = equation_regex.search(equation)
    equations.append((int(mo.group(1)),int(mo.group(2)),int(mo.group(3))))

result = solve_linear_comparisons(equations)
with open(sys.argv[2],"w") as fp:
    for x,m in result:
        print(f'[{x}]_{m}',end=' ', file=fp)