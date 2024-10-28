all:
	gcc -shared -o libmul.so mul.c
remove:
	rm -r *.so log/*.log __pycache__