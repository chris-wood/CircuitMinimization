linear:
	javac -classpath .:pj.jar MatrixOptimize.java Matrix.java MatrixState.java Pair.java Vector.java BasePair.java SLP.java MatrixOptimizeExperiment.java

search:
	javac CircuitSearch.java Matrix.java Vector.java 

gs:
	javac SignalGraph.java Matrix.java SLP.java Signal.java

nonlinear:
	javac NonlinearOptimize.java SLP.java

mine:
	g++ -o matrixOptimize matrixOptimize.cpp

canright:
	gcc -o bestboth bestboth.c

clean:
	rm -f bestboth 
	rm -f matrixOptimize
	rm *class

