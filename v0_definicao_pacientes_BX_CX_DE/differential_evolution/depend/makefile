obj = $(wildcard %.o)
%.o:	%.cpp 
	@g++ -c $< 
main:	main.o 
	@g++ *.o -o main	
run: main
	@./main
clean:
	@rm *.o
	@rm main