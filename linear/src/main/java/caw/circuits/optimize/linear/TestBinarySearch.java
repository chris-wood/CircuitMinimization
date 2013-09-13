package caw.circuits.optimize.linear;

import java.util.Random;
import java.io.BufferedReader;
import java.io.FileReader;
import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedLong;
import edu.rit.pj.reduction.SharedInteger;
import edu.rit.pj.reduction.IntegerOp;
import edu.rit.pj.Comm;
import edu.rit.pj.BarrierAction;

public class TestBinarySearch
{
	public static void error(String s) { System.err.println(s); }
	public static void disp(String s) { System.out.println(s); }

	public static void main(String[] args) throws Exception
	{
		Comm.init(args);

		if (args.length != 4 && args.length != 1)
		{
			error("usage: java TestBinarySearch [n m trials bias | filename]");
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
				SLP slp4 = MatrixOptimize.peralta_optimize(sample, sample.getDimension(), sample.getLength(), 0, 0);
				long end4 = System.currentTimeMillis();
				error("Elapsed time: " + (end4 - start4));

				error("STARTING PERALTA HEURISTIC TEST - TIE 0, GRAPH DISTANCE WITH BIT ARRAYS");
				long start5 = System.currentTimeMillis();
				SLP slp5 = MatrixOptimize.peralta_optimize(sample, sample.getDimension(), sample.getLength(), 0, 1);
				long end5 = System.currentTimeMillis();
				error("Elapsed time: " + (end5 - start5));

				error("STARTING PERALTA HEURISTIC TEST - TIE 0, SPLIT SUBSET SUM");
				long start7 = System.currentTimeMillis();
				SLP slp7 = MatrixOptimize.peralta_optimize(sample, sample.getDimension(), sample.getLength(), 0, 2);
				long end7 = System.currentTimeMillis();
				error("Elapsed time: " + (end7 - start7));

				error("STARTING PERALTA HEURISTIC TEST - TIE 0, PARALLEL RECURSIVE DISTANCE");
				long start6 = System.currentTimeMillis();
				SLP slp6 = ParallelMatrixOptimize.parallel_peralta_optimize(sample, sample.getDimension(), sample.getLength(), 0);
				long end6 = System.currentTimeMillis();
				error("Elapsed time: " + (end6 - start6));

				// Make sure we actually have the same result for both
				error("Recursive: " + slp4.xc);
				error(end4 - start4 + "ms");
				error("Graph: " + slp5.xc);
				error(end5 - start5 + "ms");
				error("Split: " + slp7.xc);
				error(end7 - end7 + "ms");
				// error("GRAPH NUMTRIES: " + MatrixOptimize.numTries);
				error("Parallel: " + slp6.xc);
				error(end6 - start6 + "ms");
				if (slp4.lines.equals(slp5.lines) == false || slp4.lines.equals(slp6.lines) == false || slp4.lines.equals(slp7.lines) == false)
				{
					error("SLPs did not match. Terminating.");
					error("" + slp4.lines.size() + "," + slp5.lines.size() + "," + slp7.lines.size() + "," + slp6.lines.size());
					error("SLP #1");
					error(slp4.lines.toString());
					error("SLP #2");
					error(slp5.lines.toString());
					error("SLP #3");
					error(slp7.lines.toString());
					error("SLP #4");
					error(slp6.lines.toString());
					System.exit(-1);
				}
			}	
		}
		else
		{
			BufferedReader in = new BufferedReader(new FileReader(args[0])); 
			String line = "";
			while ((line = in.readLine()) != null && !line.isEmpty())
			{
				String[] params = line.split(" ");
				int n = Integer.parseInt(params[0]);
				int m = Integer.parseInt(params[1]);
				int trials = Integer.parseInt(params[2]);
				double bias = Double.parseDouble(params[3]);

				// PRNG
				Random r = new Random(System.currentTimeMillis());
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

						// Run both distance computations
						error("STARTING PERALTA HEURISTIC TEST - TIE 0, RECURSIVE DISTANCE");
						long start4 = System.currentTimeMillis();
						SLP slp4 = MatrixOptimize.peralta_optimize(sample, sample.getDimension(), sample.getLength(), 0, 0);
						long end4 = System.currentTimeMillis();
						error("Elapsed time: " + (end4 - start4));

						error("STARTING PERALTA HEURISTIC TEST - TIE 0, GRAPH DISTANCE WITH BIT ARRAYS");
						long start5 = System.currentTimeMillis();
						SLP slp5 = MatrixOptimize.peralta_optimize(sample, sample.getDimension(), sample.getLength(), 0, 1);
						long end5 = System.currentTimeMillis();
						error("Elapsed time: " + (end5 - start5));

						error("STARTING PERALTA HEURISTIC TEST - TIE 0, SPLIT SUBSET SUM");
						long start7 = System.currentTimeMillis();
						SLP slp7 = MatrixOptimize.peralta_optimize(sample, sample.getDimension(), sample.getLength(), 0, 2);
						long end7 = System.currentTimeMillis();
						error("Elapsed time: " + (end7 - start7));

						error("STARTING PERALTA HEURISTIC TEST - TIE 0, PARALLEL RECURSIVE DISTANCE");
						long start6 = System.currentTimeMillis();
						SLP slp6 = ParallelMatrixOptimize.parallel_peralta_optimize(sample, sample.getDimension(), sample.getLength(), 0);
						long end6 = System.currentTimeMillis();
						error("Elapsed time: " + (end6 - start6));

						// Make sure we actually have the same result for both
						error("Recursive: " + slp4.xc);
						error("Graph: " + slp5.xc);
						error("Split: " + slp7.xc);
						error("Parallel: " + slp6.xc);
						disp(t + "," + n + "," + m + "," + bias + "," + (end4 - start4) + "," + (end5 - start5) + "," + (end7 - start7) + "," + (end6 - start6));
						if (slp4.lines.equals(slp5.lines) == false || slp4.lines.equals(slp6.lines) == false || slp4.lines.equals(slp7.lines) == false)
						{
							error("SLPs did not match. Terminating.");
							error("" + slp4.lines.size() + "," + slp5.lines.size() + "," + slp7.lines.size() + "," + slp6.lines.size());
							error("SLP #1");
							error(slp4.lines.toString());
							error("SLP #2");
							error(slp5.lines.toString());
							error("SLP #3");
							error(slp7.lines.toString());
							error("SLP #4");
							error(slp6.lines.toString());
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