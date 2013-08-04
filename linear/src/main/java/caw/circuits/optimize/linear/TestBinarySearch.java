package caw.circuits.optimize.linear;

import java.util.Random;
import java.io.*;

public class TestBinarySearch
{
	public static error(String s) { System.err.println(s); }
	public static disp(String s) { System.out.println(s); }

	public static void main(String[] args) throws Exception
	{
		if (args.length != 4 && args.length != 1)
		{
			error("usage: java TestBinarySearch n m trials bias");
			System.exit(-1);
		}

		if (args.length == 4)
		{
			int n = Integer.parseInt(args[0]);
			int m = Integer.parseInt(args[1]);
			if (n < 2 || m < 2)
			{
				error("Square matrix dimension must be > 1");
			}
			long trials = Long.parseLong(args[2]);
			if (trials < 1)
			{
				error("Number of trials must be > 0");	
			}
			double bias = Double.parseDouble(args[3]);
			if (bias < 0.0 || bias > 1.0)
			{
				error("Bias must be in the range [0.0,1.0]");
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
				error("" + sample.toString());

				// Run both distance computations
				error("STARTING PERALTA HEURISTIC TEST - TIE 0, RECURSIVE DISTANCE");
				long start4 = System.currentTimeMillis();
				SLP slp4 = MatrixOptimize.peraltaOptimize(sample, sample.getDimension(), sample.getLength(), 0, 0);
				long end4 = System.currentTimeMillis();
				error("Elapsed time: " + (end4 - start4));

				error("STARTING PERALTA HEURISTIC TEST - TIE 0, GRAPH DISTANCE WITH BIT ARRAYS");
				long start5 = System.currentTimeMillis();
				SLP slp5 = MatrixOptimize.peraltaOptimize(sample, sample.getDimension(), sample.getLength(), 0, 1);
				long end5 = System.currentTimeMillis();
				error("Elapsed time: " + (end5 - start5));

				// Make sure we actually have the same result for both
				disp("Recursive: " + slp4.xc);
				disp("Graph: " + slp5.xc);
				if (slp4.lines.equals(slp5.lines) == false)
				{
					error("SLPs did not match. Terminating.");
					error("" + slp4.lines.size() + "," + slp5.lines.size());
					error("SLP #1");
					error(slp4.lines);
					error("SLP #2");
					error(slp5.lines);
					System.exit(-1);
				}
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
					try 
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

						// public static int peraltaDistance(final Matrix base, final int[] f)

						// Run both distance computations
						error("STARTING PERALTA HEURISTIC TEST - TIE 0, RECURSIVE DISTANCE");
						long start4 = System.currentTimeMillis();
						SLP slp4 = MatrixOptimize.peraltaOptimize(sample, sample.getDimension(), sample.getLength(), 0, 0);
						long end4 = System.currentTimeMillis();
						error("Elapsed time: " + (end4 - start4));

						error("STARTING PERALTA HEURISTIC TEST - TIE 0, GRAPH DISTANCE WITH BIT ARRAYS");
						long start5 = System.currentTimeMillis();
						SLP slp5 = MatrixOptimize.peraltaOptimize(sample, sample.getDimension(), sample.getLength(), 0, 1);
						long end5 = System.currentTimeMillis();
						error("Elapsed time: " + (end5 - start5));

						// display the times in CSV format
						disp(itr + "," + n + "," + m + "," + bias + "," + (end4 - start4) + "," + (end5 - start5));

						// Make sure we actually have the same result for both
						if (slp4.lines.equals(slp5.lines) == false)
						{
							error("SLPs did not match. Terminating.");
							error("" + slp4.lines.size() + "," + slp5.lines.size());
							error("SLP #1");
							error(slp4.lines);
							error("SLP #2");
							error(slp5.lines);
							System.exit(-1);
						}
					}
					catch (Exception e)
					{
						t--;
						continue;
					}
				}	
			}
		}
	}
}