import java.util.ArrayList;
import java.util.Scanner;
import edu.rit.pj.Comm;
import edu.rit.util.Random;

public class PerformanceTest
{

	// java -Xmx2000m -Dpj.schedule=guided PerformanceTest 18 18 1 0.5 123
	// java -Xmx2000m -Dpj.schedule=guided PerformanceTest 17 17 1 0.5 123

	static void error(String s) { System.err.println(s); }
	static void disp(String s) { System.out.println(s); }

	public static boolean inRows(ArrayList<int[]> rows, int[] row)
	{
		for (int[] r : rows)
		{
			boolean match = true;
			for (int i = 0; i < r.length; i++)
			{
				if (r[i] != row[i]) 
				{
					match = false;
					break;
				}
			}
			if (match) return true;
		}
		return false;
	}

	public static void main(String[] args) throws Exception
	{
		Comm.init(args);
		if (args.length != 5) 
		{
			System.err.println("usage: java PerformanceTest n m trials bias seed");
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
				results.add(new ArrayList<Integer>()); // peraltaOptimize 0
				results.add(new ArrayList<Integer>()); // peraltaOptimize 1
				results.add(new ArrayList<Integer>()); // peraltaOptimize 2
				results.add(new ArrayList<Integer>()); // peraltaOptimize 3

				// Create the matrix
				ArrayList<int[]> rows = new ArrayList<int[]>();
				boolean allZero = true;
				for (int i = 0; i < n; i++)
				{
					int[] row = new int[m];
					do 
					{
						row = new int[m];
						for (int j = 0; j < m; j++)
						{
							if (r.nextDouble() <= bias)
							{
								allZero = false;
								row[j] = 1;
							}
						}
					} while (inRows(rows, row));
					rows.add(row);
				}
				Matrix sample = new Matrix(rows);
				SLP slp1;
				
				// slp1 = null;
				// for (long itr = 0; itr < trials; itr++)
				// {
				// 	long start1 = System.currentTimeMillis();
				// 	SLP slp4 = MatrixOptimize.parallelPeraltaOptimize_v2(sample, sample.getDimension(), sample.getLength(), 0);
				// 	// results.get(3).add(slp4.xc);
				// 	long end1 = System.currentTimeMillis();

				// 	if (slp1 == null) slp1 = slp4;
				// 	else 
				// 	{
				// 		// disp("stuff");
				// 		// disp(slp4.toString());
				// 		// if (!(slp1.equals(slp4)))
				// 		if (!(slp1.equals(slp4)))
				// 		{
				// 			throw new Exception("Results from the same code don't match...");
				// 		}
				// 	}

				// 	System.out.println("0 p " + (end1 - start1) + " ms");
				// }
				
				// // sample = new Matrix(rows);
				// for (long itr = 0; itr < trials; itr++)
				// {
				// 	long sstart1 = System.currentTimeMillis();
				// 	SLP slp4 = MatrixOptimize.peraltaOptimize(sample, sample.getDimension(), sample.getLength(), 0);
				// 	// results.get(3).add(slp4.xc);
				// 	long send1 = System.currentTimeMillis();

				// 	// disp("sequential stuff");
				// 	// disp(slp4.toString());
				// 	if (!(slp1.equals(slp4)))
				// 	{
				// 		throw new Exception("Parallel and sequential results don't match...");
				// 	}

				// 	System.out.println("0 s " + (send1 - sstart1) + " ms");
				// }
				
				// slp1 = null;
				// for (long itr = 0; itr < trials; itr++)
				// {
				// 	long start2 = System.currentTimeMillis();
				// 	SLP slp4 = MatrixOptimize.parallelPeraltaOptimize_v2(sample, sample.getDimension(), sample.getLength(), 1);
				// 	// results.get(4).add(slp4.xc);
				// 	long end2 = System.currentTimeMillis();

				// 	if (slp1 == null) slp1 = slp4;
				// 	else 
				// 	{
				// 		if (!(slp1.equals(slp4)))
				// 		{
				// 			// disp("old stuff");
				// 			// disp(slp1.toString());
				// 			// disp("new stuff");
				// 			// disp(slp4.toString());
				// 			throw new Exception("Results from the same code don't match...");
				// 		}
				// 	}

				// 	System.out.println("1 p " + (end2 - start2) + " ms");
				// }
				
				// for (long itr = 0; itr < trials; itr++)
				// {
				// 	long sstart2 = System.currentTimeMillis();
				// 	SLP slp4 = MatrixOptimize.peraltaOptimize(sample, sample.getDimension(), sample.getLength(), 1);
				// 	// results.get(4).add(slp4.xc);
				// 	long send2 = System.currentTimeMillis();

				// 	if (!(slp1.equals(slp4)))
				// 	{
				// 		// disp("parallel stuff");
				// 		// disp(slp1.toString());
				// 		// disp("sequential stuff");
				// 		// disp(slp4.toString());
				// 		throw new Exception("Parallel and sequential results don't match...");
				// 	}

				// 	System.out.println("1 s " + (send2 - sstart2) + " ms");
				// }
				
				slp1 = null;
				for (long itr = 0; itr < trials; itr++)
				{
					long start3 = System.currentTimeMillis();
					SLP slp4 = MatrixOptimize.parallelPeraltaOptimize_v2(sample, sample.getDimension(), sample.getLength(), 2);
					// results.get(5).add(slp4.xc);
					long end3 = System.currentTimeMillis();

					if (slp1 == null) slp1 = slp4;
					else 
					{
						if (!(slp1.equals(slp4)))
						{
							throw new Exception("Results from the same code don't match...");
						}
					}

					System.out.println("2 p " + (end3 - start3) + " ms");
				}
				
				
				for (long itr = 0; itr < trials; itr++)
				{
					long sstart3 = System.currentTimeMillis();
					SLP slp4 = MatrixOptimize.peraltaOptimize(sample, sample.getDimension(), sample.getLength(), 2);
					// results.get(5).add(slp4.xc);
					long send3 = System.currentTimeMillis();

					if (!(slp1.equals(slp4)))
					{
						throw new Exception("Parallel and sequential results don't match...");
					}

					System.out.println("2 s " + (send3 - sstart3) + " ms");
				}
				
				// // slp1 = null;
				// for (long itr = 0; itr < trials; itr++)
				// {
				// 	long start4 = System.currentTimeMillis();
				// 	SLP slp4 = MatrixOptimize.parallelPeraltaOptimize_v2(sample, sample.getDimension(), sample.getLength(), 3);
				// 	// results.get(6).add(slp4.xc);
				// 	long end4 = System.currentTimeMillis();

				// 	// if (slp1 == null) slp1 = slp4;
				// 	// else 
				// 	// {
				// 	// 	if (!(slp1.equals(slp4)))
				// 	// 	{
				// 	// 		throw new Exception("REsults don't match...");
				// 	// 	}
				// 	// }

				// 	System.out.println("3 p " + (end4 - start4) + " ms");
				// }
				
				
				// for (long itr = 0; itr < trials; itr++)
				// {
				// 	long sstart4 = System.currentTimeMillis();
				// 	SLP slp4 = MatrixOptimize.peraltaOptimize(sample, sample.getDimension(), sample.getLength(), 3);
				// 	// results.get(6).add(slp4.xc);
				// 	long send4 = System.currentTimeMillis();

				// 	// don't compare - this is the randomized heuristic
				// 	// if (!(slp1.equals(slp4)))
				// 	// {
				// 	// 	throw new Exception("REsults don't match...");
				// 	// }

				// 	System.out.println("3 s " + (send4 - sstart4) + " ms");
				// }
				
			} 
			catch (Exception e)
			{
				e.printStackTrace();
				error(e.getMessage());
			}
		}
	}
}


