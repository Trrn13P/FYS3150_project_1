all: execute compile execute



compile:
	c++ -o ./main.out main.cpp -larmadillo

execute:
	rm -r ./data; mkdir ./data
	./main.out
	rm -r ./figures; mkdir ./figures
	#python3 reader.py
