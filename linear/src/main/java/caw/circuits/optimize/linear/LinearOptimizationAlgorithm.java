package caw.circuits.optimize.linear;

public interface LinearOptimizationAlgorithm
{
	// Each row in this matrix is a linear form in the system to optimize
	public SLP optimize(Matrix forms); 

	// TODO: implement paar, paar exhaustive, peralta recursive distance, peralta brute distance, peralta graph distance, exhaustive factor
}
