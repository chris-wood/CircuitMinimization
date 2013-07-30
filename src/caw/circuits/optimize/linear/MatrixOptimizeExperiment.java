import java.util.ArrayList;
import java.util.Scanner;
import java.util.Random;

public class MatrixOptimizeExperiment
{

	static void error(String s) { System.err.println(s); }
	static void disp(String s) { System.out.println(s); }

	public static void main(String[] args) throws Exception
	{
		if (args.length != 5) 
		{
			System.err.println("usage: java MatrixOptimizeExperiment n m trials bias seed");
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
				Random r = new Random(seed);

				// Storage containers for the number of gates
				ArrayList<ArrayList<Integer>> results = new ArrayList<ArrayList<Integer>>();
				results.add(new ArrayList<Integer>()); // paar
				results.add(new ArrayList<Integer>()); // paar exhaustive
				results.add(new ArrayList<Integer>()); // exhaustive factor
				results.add(new ArrayList<Integer>()); // peraltaOptimize 0
				results.add(new ArrayList<Integer>()); // peraltaOptimize 1
				results.add(new ArrayList<Integer>()); // peraltaOptimize 2
				results.add(new ArrayList<Integer>()); // peraltaOptimize 3
				results.add(new ArrayList<Integer>()); // bernstein

				for (long itr = 0; itr < trials; itr++)
				{
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
					Matrix sample = new Matrix(mat, m);
					error("" + sample.toString());

					// Run the experiments
					if (!allZero)
					{
						// TODO: create matrices, run optimizations, store result in a list
						// disp(sample.toString());
						error("STARTING PAAR TEST");
						ArrayList<MatrixState> history = MatrixOptimize.paarOptimize(sample, sample.getLength()); 
						results.get(0).add(MatrixOptimize.numGates(sample.getDimension(), history));
						// dispStrings(buildSLP(history));

						// error("STARTING PAAR PSEUDO-EXHAUSTIVE TEST");
						// ArrayList<MatrixState> history1 = MatrixOptimize.paarOptimizeExhaustive(sample, sample.getLength(), 0, 0, sample.getLength() - 1); // 5x5 matrix at this point...
						// results.get(1).add(MatrixOptimize.numGates(sample.getDimension(), history1));
						// dispStrings(buildSLP(history1));

						// disp("STARTING EXHAUSTIVE FACTORIZATION TEST");
						// ArrayList<MatrixState> history3 = MatrixOptimize.exhaustiveOptimize(sample, 0, 0, sample.getLength(), sample.getLength() - 1, new ArrayList<Pair>()); // 5x5 matrix at this point...
						// results.get(2).add(MatrixOptimize.numGates(sample.getDimension(), history3));
						// dispStrings(buildSLP(history3));

						error("STARTING PERALTA TEST - TIE 0");
						SLP slp4 = MatrixOptimize.parallelPeraltaOptimize_v2(sample, sample.getDimension(), sample.getLength(), 0);
						results.get(3).add(slp4.xc);
						error("" + slp4.xc);
						// dispStrings(slp4.lines);

						// error("STARTING PERALTA TEST - TIE 1");
						// slp4 = MatrixOptimize.parallelPeraltaOptimize_v2(sample, sample.getDimension(), sample.getLength(), 1);
						// results.get(4).add(slp4.xc);
						// error("" + slp4.xc);
						// // dispStrings(slp4.lines);

						// error("STARTING PERALTA TEST - TIE 2");
						// slp4 = MatrixOptimize.parallelPeraltaOptimize_v2(sample, sample.getDimension(), sample.getLength(), 2);
						// results.get(5).add(slp4.xc);
						// error("" + slp4.xc);
						// // dispStrings(slp4.lines);

						// error("STARTING PERALTA TEST - TIE 3");
						// slp4 = MatrixOptimize.parallelPeraltaOptimize_v2(sample, sample.getDimension(), sample.getLength(), 3);
						// results.get(6).add(slp4.xc);
						// error("" + slp4.xc);
						// dispStrings(slp4.lines);

						error("STARTING BERNSTEIN TEST");
						ArrayList<String> slp5 = MatrixOptimize.bernsteinOptimize(sample.getVectors(), sample.getLength(), sample.getLength(), sample.getDimension());
						results.get(7).add(slp5.size());
						error("" + slp5.size());
						// dispStrings(slp5);
					}
				}

				// Compute and display the averages
				int[] totals = new int[results.size()];
				for (int i = 0; i < results.size(); i++)
				{
					for (int j = 0; j < results.get(i).size(); j++)
					{
						totals[i] += results.get(i).get(j);
					}
				}
				for (int i = 0; i < totals.length; i++)
				{
					// double avg = (double)totals[i] / (double)trials;
					disp(i + ": " + totals[i]);
				}
			} 
			catch (Exception e)
			{
				error(e.getMessage());
			}
		}
	}
}
