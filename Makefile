all:
	gcc -shared -o libmul.so mul.c
remove:
	rm *.so log/*.log