package caw.circuits.optimize.linear;

import java.util.Random;

public class TestBinarySearch
{
	public static void main(String[] args) throws Exception
	{
		if (args.length != 4)
		{
			System.err.println("usage: java TestBinarySearch n m trials bias");
			System.exit(-1);
		}

		int n = Integer.parseInt(args[0]);
		int m = Integer.parseInt(args[1]);
		if (n < 2 || m < 2)
		{
			System.err.println("Square matrix dimension must be > 1");
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
		long seed = System.currentTimeMillis();

		// PRNG
		Random r = new Random(seed);
		for (int t = 0; t < trials; t++)
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
			System.err.println("" + sample.toString());

			// Testing against our original sequential code...
			int[][] b = new int[m][m];
			for (int i = 0; i < m; i++)
			{
				for (int j = 0; j < m; j++)
				{
					if (i == j)
					{
						b[i][j] = 1;
					}
				}
			}
			Matrix base = new Matrix(b);
			int[] f = new int[m];
			Random prng = new Random();
			for (int i = 0; i < m; i++) f[i] = prng.nextInt(2); 
			boolean unique = false;
			while (!unique)
			{
				boolean match = false;
				for (int row = 0; row < base.getDimension(); row++)
				{
					if (MatrixOptimize.areEqual(base.getRow(row), f))
					{
						match = true;
						break;
					}
				}
				if (match == false)
				{
					unique = true;
				}
				else
				{
					// repopulate the thing
					for (int i = 0; i < m; i++) f[i] = prng.nextInt() % 2; 
				}
			}

			System.out.println("Target:");
			System.out.print("[ ");
			for (int i = 0; i < m - 1; i++) System.out.print(f[i] + " ");
			System.out.println(f[f.length - 1] + "]");

			// compute the distances
			System.out.println("Brute: " + MatrixOptimize.peraltaDistance(base, f));

			// public static int peraltaDistance(final Matrix base, final int[] f)

			// Run both distance computations
			System.err.println("STARTING PERALTA HEURISTIC TEST - TIE 0, RECURSIVE DISTANCE");
			long start4 = System.currentTimeMillis();
			SLP slp4 = MatrixOptimize.peraltaOptimize(sample, sample.getDimension(), sample.getLength(), 0, 0);
			long end4 = System.currentTimeMillis();
			System.err.println("Elapsed time: " + (end4 - start4));

			System.err.println("STARTING PERALTA HEURISTIC TEST - TIE 0, GRAPH DISTANCE WITH BIT ARRAYS");
			long start5 = System.currentTimeMillis();
			SLP slp5 = MatrixOptimize.peraltaOptimize(sample, sample.getDimension(), sample.getLength(), 0, 1);
			long end5 = System.currentTimeMillis();
			System.err.println("Elapsed time: " + (end5 - start5));

			// System.err.println("STARTING PERALTA HEURISTIC TEST - TIE 0, GRAPH DISTANCE WITH UNSIGNED INTEGERS");
			// long start6 = System.currentTimeMillis();
			// SLP slp6 = MatrixOptimize.peraltaOptimize(sample, sample.getDimension(), sample.getLength(), 0, 2);
			// long end6 = System.currentTimeMillis();
			// System.err.println("Elapsed time: " + (end6 - start6));

			// Make sure we actually have the same result for both
			System.out.println("Recursive: " + slp4.xc);
			System.out.println("Graph: " + slp5.xc);
			if (slp4.lines.equals(slp5.lines) == false)
			{
				System.err.println("SLPs did not match. Terminating.");
				System.err.println("" + slp4.lines.size() + "," + slp5.lines.size());
				System.err.println("SLP #1");
				System.err.println(slp4.lines);
				System.err.println("SLP #2");
				System.err.println(slp5.lines);
				System.exit(-1);
			}

			// Display results
			// System.out.println(n + "," + m + "," + (end4 - start4) + "," + (end5 - start5));
		}
	}
}