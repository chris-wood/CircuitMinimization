package caw.circuits.optimize.linear;

import java.util.ArrayList;
import java.util.Random;
import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedLong;
import edu.rit.pj.reduction.SharedInteger;
import edu.rit.pj.reduction.IntegerOp;
import edu.rit.pj.Comm;
import edu.rit.pj.BarrierAction;

public class ParallelMatrixOptimize
{
	private static SolutionCircuit solutionCircuit = new SolutionCircuit();

	static int p_xorCount = 0;
	static int[] p_dist;
	static ArrayList<String> p_slp;
	static Matrix p_base;
	static int p_n;
	static int[] p_newDist;
	static int p_i;
	static int p_ii;
	static int p_jj;
	static int p_d;
	static int[] p_optimalBase;
	static ArrayList<int[]> p_optimalBases;
	static ArrayList<int[]> p_optimalDistances;
	static ArrayList<Pair> p_pairs;
	static ParallelTeam team = new ParallelTeam();

	public static SLPStrings parallel_peralta_optimize(final Matrix m, final int r, final int c, final int tieBreaker) throws Exception
	{
		p_xorCount = 0;
		solutionCircuit.reset();
		p_dist = null;
		p_slp = null;
		p_base = null;
		p_n = 0;
		p_newDist = null;
		p_optimalBases = null;
		p_optimalDistances = null;
		p_pairs = null;

		p_slp = new ArrayList<String>();
		// Create the initial base
		int[][] b = new int[c][c];
		for (int i = 0; i < c; i++)
		{
			for (int j = 0; j < c; j++)
			{
				if (i == j)
				{
					b[i][j] = 1;
				}
			}
		}
		p_base = new Matrix(b, c);
		
		// Initialize the distance array and then get the ball rolling
		p_dist = MatrixOptimize.compute_distance(p_base, m);
		for (int f = 0; f < p_dist.length; f++)
		{
			if (p_dist[f] == 0)
			{
				int[] row = m.getRow(f);
				for (int cc = 0; cc < row.length; cc++)
				{
					if (row[cc] == 1)
					{
						p_slp.add("y_" + f + " = " + "x_" + cc);
					}
				}
			}
		}
		p_newDist = new int[m.getDimension()];
		p_i = m.getLength();
		final int nVars = m.getLength();

		p_n = p_base.getDimension(); // each row is a base...
		p_d = MatrixOptimize.sum(p_dist);
		p_optimalBase = new int[p_base.getLength()]; // was dimension
		p_optimalBases = new ArrayList<int[]>();
		p_optimalDistances = new ArrayList<int[]>();
		p_pairs = new ArrayList<Pair>();

		// disp("Entering: " + p_n);

		// Create the parallel team so everything isn't anonymous
		ParallelTeam team = new ParallelTeam();

		while (MatrixOptimize.is_zero(p_dist) == false)
		{
			// disp("Here: " + p_n);
			// error("Another go round");
			team.execute(new ParallelRegion() // linear combinations get computed in parallel - but can we push parallelism to the outer loop?
			{
				public void run() throws Exception // There are many, many combinations (potentially...)
				{
					execute (0, p_n - 1, new IntegerForLoop()
					{
						int t_d = p_d;
						ArrayList<int[]> t_optimalBases = new ArrayList<int[]>();
						ArrayList<int[]> t_optimalDistances = new ArrayList<int[]>();
						ArrayList<Pair> t_pairs = new ArrayList<Pair>();

						public IntegerSchedule schedule()
						{
							return IntegerSchedule.guided();
						}

						public void run (int first, int last) throws Exception
						{
							for (int t_ii = first; t_ii <= last; t_ii++)
							{
								for (int t_jj = t_ii + 1; t_jj < p_n; t_jj++)
								{
									if (t_ii != t_jj)
									{
										// t_newDist = new int[m.getDimension()];

										int[] ssum = MatrixOptimize.add_mod(p_base.getRow(t_ii), p_base.getRow(t_jj), 2);
										if (!(p_base.containsRow(ssum)) && !MatrixOptimize.is_zero(ssum))
										{												
											// int[] innersum = add_mod(p_base.getRow(t_ii), p_base.getRow(t_jj), 2);
											// // Matrix newBase = p_base.copy();
											// // newBase.appendRow(innersum);
											// for (int k = 0; k < last; k++)
											// {
											// 	int[] t_fi = m.getRow(k);
											// 	// p_newDist[k] = optimizedPeraltaDistance(newBase, t_fi);
											// 	t_newDist[k] = optimizedPeraltaDistance(p_base, innersum, t_fi, p_dist[k]);
											// }
											int[] t_newDist = MatrixOptimize.optimized_compute_distance(p_base, m, ssum, p_dist);
											int newDistSum = MatrixOptimize.sum(t_newDist);
											// disp("Sum for " + t_ii + "," + t_jj + " = " + newDistSum);
											if (newDistSum < t_d)
											{
												t_d = newDistSum;
												t_optimalBases = new ArrayList<int[]>();
												t_optimalDistances = new ArrayList<int[]>();
												t_pairs = new ArrayList<Pair>();
												t_optimalBases.add(ssum);
												t_optimalDistances.add(t_newDist);
												t_pairs.add(new Pair(t_ii, t_jj));
											}
											else if (newDistSum == t_d)
											{
												t_optimalBases.add(ssum);
												t_optimalDistances.add(t_newDist);
												t_pairs.add(new Pair(t_ii, t_jj));
											}
											ssum = null;
										}	
									}
								}
							}
						}

						public void finish() throws Exception
						{
							if (t_d < p_d)
							{
								solutionCircuit.record(t_d, t_optimalBases, t_optimalDistances, t_pairs);
							}
						}
					});
				}
			});
	
			// Copy over the solution...
			p_d = solutionCircuit.d;	
			p_optimalBases = solutionCircuit.optimalBases;
			p_optimalDistances = solutionCircuit.optimalDistances;
			p_pairs = solutionCircuit.pairs;
	
			// disp("Executing the barrier action for tie breaking and round updates");
			// Resolve using norms...
			int mi = 0;
			int mj = 0;
			double maxNorm = Double.MAX_VALUE;
			switch (tieBreaker)
			{
				case 0: // norm
					maxNorm = MatrixOptimize.euclidean_norm(p_optimalDistances.get(0));
					p_optimalBase = p_optimalBases.get(0);
					p_newDist = p_optimalDistances.get(0);
					mi = p_pairs.get(0).i;
					mj = p_pairs.get(0).j;
					for (int i = 1; i < p_optimalDistances.size(); i++)
					{
						double norm = MatrixOptimize.euclidean_norm(p_optimalDistances.get(i));
						if (norm > maxNorm)
						{
							maxNorm = norm;
							p_optimalBase = p_optimalBases.get(i);
							p_newDist = p_optimalDistances.get(i);
							mi = p_pairs.get(i).i;
							mj = p_pairs.get(i).j;
						}
						else if (norm == maxNorm && MatrixOptimize.compare_pairs(p_pairs.get(i).i, p_pairs.get(i).j, mi, mj) <= 0)
						{
							maxNorm = norm;
							p_optimalBase = p_optimalBases.get(i);
							p_newDist = p_optimalDistances.get(i);
							mi = p_pairs.get(i).i;
							mj = p_pairs.get(i).j;	
						}
					}
					break;
				case 1: // norm-largest
					int largestDist = MatrixOptimize.largest(p_optimalDistances.get(0));
					// for (int i = 0; i < p_optimalDistances.get(0).length; i++)
					// {
					// 	int tdist = p_optimalDistances.get(0)[i];
					// 	largestDist = tdist > largestDist ? tdist : largestDist;
					p_optimalBase = p_optimalBases.get(0);
					p_newDist = p_optimalDistances.get(0);
					mi = p_pairs.get(0).i;
					mj = p_pairs.get(0).j;
					// }
					maxNorm = Math.pow(MatrixOptimize.euclidean_norm(p_optimalDistances.get(0)), 2);
					maxNorm = maxNorm - (double)largestDist;
					for (int i = 1; i < p_optimalDistances.size(); i++)
					{
						double norm = Math.pow(MatrixOptimize.euclidean_norm(p_optimalDistances.get(i)), 2);
						norm = norm - (double)MatrixOptimize.largest(p_optimalDistances.get(i));
						if (norm > maxNorm)
						{
							maxNorm = norm;
							p_optimalBase = p_optimalBases.get(i);
							p_newDist = p_optimalDistances.get(i);
							mi = p_pairs.get(i).i;
							mj = p_pairs.get(i).j;
						}
						else if (norm == maxNorm && MatrixOptimize.compare_pairs(p_pairs.get(i).i, p_pairs.get(i).j, mi, mj) <= 0)
						{
							maxNorm = norm;
							p_optimalBase = p_optimalBases.get(i);
							p_newDist = p_optimalDistances.get(i);
							mi = p_pairs.get(i).i;
							mj = p_pairs.get(i).j;	
						}
					}
					break;
				case 2: // norm-diff
					int firstLargestDist = 0;
					int secondLargestDist = -1;
					mi = p_pairs.get(0).i;
					mj = p_pairs.get(0).j;
					for (int i = 0; i < p_optimalDistances.get(0).length; i++)
					{
						int tdist = p_optimalDistances.get(0)[i];
						p_optimalBase = p_optimalBases.get(0);
						p_newDist = p_optimalDistances.get(0);
						if (tdist >= firstLargestDist)
						{
							secondLargestDist = firstLargestDist;
							firstLargestDist = tdist;
						}
					}

					maxNorm = Math.pow(MatrixOptimize.euclidean_norm(p_optimalDistances.get(0)), 2);
					maxNorm = maxNorm - (double)(firstLargestDist - secondLargestDist);
					for (int i = 1; i < p_optimalDistances.size(); i++)
					{
						double norm = Math.pow(MatrixOptimize.euclidean_norm(p_optimalDistances.get(i)), 2);
						norm = norm - (double)(firstLargestDist - secondLargestDist);
						if (norm > maxNorm)
						{
							maxNorm = norm;
							p_optimalBase = p_optimalBases.get(i);
							p_newDist = p_optimalDistances.get(i);
							mi = p_pairs.get(i).i;
							mj = p_pairs.get(i).j;
						}
						else if (norm == maxNorm && MatrixOptimize.compare_pairs(p_pairs.get(i).i, p_pairs.get(i).j, mi, mj) <= 0)
						{
							maxNorm = norm;
							p_optimalBase = p_optimalBases.get(i);
							p_newDist = p_optimalDistances.get(i);
							mi = p_pairs.get(i).i;
							mj = p_pairs.get(i).j;	
						}
					}
					break;
				case 3: // random
					Random rng = new Random(System.currentTimeMillis());
					int i = rng.nextInt(p_optimalDistances.size());
					p_optimalBase = p_optimalBases.get(i);
					p_newDist = p_optimalDistances.get(i);
					mi = p_pairs.get(i).i;
					mj = p_pairs.get(i).j;
					break;
				default: // just pick the first one
					p_optimalBase = p_optimalBases.get(0);
					p_newDist = p_optimalDistances.get(0);
					break;
			}

			// disp("picked " + mi + " " + mj + ": " + maxNorm);

			// return new BasePair(optimalBase, new Pair(mi, mj));
			BasePair newBase = new BasePair(p_optimalBase, p_newDist, new Pair(mi, mj));
			// disp("OPTIMAL BASE");
			// disp(p_optimalBase);
			// disp("TRYING TO APPEND WITH i = " + p_i);
			p_base.appendRow(newBase.base);

			// Insert the new line
			String c1 = newBase.p.i < nVars ? "x" : "t";
			String c2 = newBase.p.j < nVars ? "x" : "t";
			p_slp.add("t_" + p_i + " = " + c1 + newBase.p.i + " XOR " + c2 + newBase.p.j);
			p_xorCount++;
			for (int f = 0; f < m.getDimension(); f++)
			{
				if (MatrixOptimize.are_equal(newBase.base, m.getRow(f)))
				{
					p_slp.add("y_" + f + " = " + "t_" + p_i);
				}
			}
			p_i++;

			// Update variables for the next round...
			p_n = p_base.getDimension(); // each row is a base...
			// disp("UPDATED P_N: " + p_n);
			// BasePair newBase = pickBase(base, m, dist, tieBreaker);
			p_d = MatrixOptimize.sum(p_dist);
			p_optimalBase = null;
			p_optimalBase = new int[p_base.getLength()]; // was dimension
			p_optimalBases = null;
			p_optimalBases = new ArrayList<int[]>();
			p_optimalDistances = null;
			p_optimalDistances = new ArrayList<int[]>();
			p_pairs = null;
			p_pairs = new ArrayList<Pair>();

			for (int k = 0; k < p_newDist.length; k++)
			{
				p_dist[k] = p_newDist[k];
			}
			// disp("NEW BASE: ");
			// disp(p_dist);
		}

		return new SLPStrings(p_slp, p_xorCount, 0); // no ANDs required
	}
}