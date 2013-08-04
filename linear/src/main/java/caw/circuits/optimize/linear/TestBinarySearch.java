package caw.circuits.optimize.linear;

import java.util.Random;
import java.io.*;

public class TestBinarySearch
{
	public static void main(String[] args) throws Exception
	{
		if (args.length != 4 && args.length != 1)
		{
			System.err.println("usage: java TestBinarySearch n m trials bias");
			System.exit(-1);
		}

		if (args.length == 4)
		{
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
		else
		{
			BufferedReader in = new BufferedReader(new FileReader(file)); 
			String line = "";
			while ((line = in.readLine()) != null && !line.isEmpty())
			{
				String params = line.split(" ");
				int n = Integer.parseInt(params[0]);
				int m = Integer.parseInt(params[1]);
				int trials = Integer.parseInt(params[2]);
				double bias = Double.parseDouble(params[3]);

				// PRNG
				Random r = new Random(System.currentTimeMillis(););
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
	}
}