import java.util.ArrayList;
import java.util.Scanner;
import edu.rit.pj.Comm;
import edu.rit.util.Random;
import java.util.concurrent.*;

public class FactorPerformanceTest
{

	// java -classpath .;pj.jar FactorPerformanceTest

	static void error(String s) { System.err.println(s); }
	static void disp(String s) { System.out.println(s); }

	public static void main(String[] args) throws Exception
	{
		Comm.init(args);
		if (args.length != 5) 
		{
			System.err.println("usage: java FactorPerformanceTest n m trials bias seed");
		}
		else
		{
			try
			{
				int n = Integer.parseInt(args[0]);
				int m = Integer.parseInt(args[1]);
				if (n < 2 || m < 2)
				{
					System.err.println("Sqaure matrix dimension must be > 1");
				}
				long trials = Long.parseLong(args[2]);
				if (trials < 1)
				{
					System.err.println("Number of trials must be > 0");	
				}
				double bias = Double.parseDouble(args[3]);
				if (bias < 0.0 || bias > 1.0)
				{
					System.err.println("Bias must be in the range [0.0,1.0]");
				}
				long seed = Long.parseLong(args[4]);

				// PRNG
				Random r = Random.getInstance(seed);

				// Storage containers for the number of gates
				ArrayList<ArrayList<Integer>> results = new ArrayList<ArrayList<Integer>>();
				results.add(new ArrayList<Integer>()); // parallel factor
				results.add(new ArrayList<Integer>()); // exhaustive factor

				// Create the matrix
				int[][] mat = new int[n][m];
				boolean allZero = true;
				for (int i = 0; i < n; i++)
				{
					for (int j = 0; j < m; j++)
					{
						if (r.nextDouble() <= bias)
						{
							allZero = false;
							mat[i][j] = 1;
						}
						else
						{
							mat[i][j] = 0;
						}
					}
				}
				Matrix sample = new Matrix(mat);

				int optimalGates = -1;
				for (long itr = 0; itr < trials; itr++)
				{
					ParallelSolution solution = new ParallelSolution(sample.getGateCount(), sample, 0, 0, n, m - 1);
					ParallelMatrixOptimize optimizer = new ParallelMatrixOptimize(solution);
					long start = System.currentTimeMillis();
					ForkJoinPool pool = new ForkJoinPool();
					pool.invoke(optimizer);
					long end = System.currentTimeMillis();

					if (optimalGates == -1) optimalGates = optimizer.solution.gOptimal;
					else 
					{
						if (optimalGates != optimizer.solution.gOptimal)
						{
							throw new Exception("Results don't match...");
						}
					}

					System.out.println("0 p " + (end - start) + " ms");
				}

				for (long itr = 0; itr < trials; itr++)
				{
					long start = System.currentTimeMillis();
					ArrayList<MatrixState> history3 = MatrixOptimize.exhaustiveOptimize(sample, 0, 0, n, m - 1); // 5x5 matrix at this point...
					long end = System.currentTimeMillis();
					int gc = history3.get(history3.size() - 1).matrix.getGateCount();

					if (optimalGates != gc)
					{
						throw new Exception("Results don't match...");
					}

					System.out.println("0 s " + (end - start) + " ms");
				}
			} 
			catch (Exception e)
			{
				e.printStackTrace();
				error(e.getMessage());
			}
		}
	}
}


