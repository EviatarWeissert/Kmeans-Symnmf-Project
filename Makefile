
symnmf: symnmf.o symnmf.h
	@echo "Building symnmf"
	@gcc -o symnmf symnmf.o -lm

symnmf.o: symnmf.c
	@echo "Compiling source file"
	@gcc -ansi -Wall -Wextra -Werror -pedantic-errors -c symnmf.c 

clean:
	@echo "Removing compiled files..."
	@rm -f *.o symnmf
