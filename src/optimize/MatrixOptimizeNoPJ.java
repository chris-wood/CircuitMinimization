import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.util.Deque;
import java.util.ArrayDeque;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;

public class MatrixOptimizeNoPJ
{
	static void error(String s) { System.err.println(s); }
	static void disp(String s) { System.out.println(s); }
	static void disp(int n) { System.out.println(n); }
	static void disp(int[] v) {
		System.out.print("[");
		for (int i = 0; i < v.length - 1; i++) System.out.print(v[i] + " ");
		System.out.println(v[v.length - 1] + "]");
	}
	static void disp(ArrayList<int[]> vs) {
		for (int i = 0; i < vs.size(); i++) disp(vs.get(i));
	}
	static void dispStrings(ArrayList<String> vs) {
		for (int i = 0; i < vs.size(); i++) disp(vs.get(i));
	}
	static void dispVectors(ArrayList<Vector> vs) {
		for (int i = 0; i < vs.size(); i++) disp(vs.get(i).row);
	}
	static ArrayList<Pair> copy(ArrayList<Pair> pairs) {
		ArrayList<Pair> np = new ArrayList<Pair>();
		for (int i = 0; i < pairs.size(); i++) np.add(pairs.get(i));
		return np;
	}

	public static double taxiNorm(int[] v)
	{
		return (double)sum(v);
	}

	public static double euclideanNorm(int[] v)
	{
		double sum = 0.0;
		for (int i = 0; i < v.length; i++)
		{
			sum += Math.pow(v[i], 2);
		}
		return Math.sqrt(sum);
	}

	public static int largest(int[] v)
	{
		int l = v[0];
		for (int i = 1; i < v.length; i++)
		{
			l = v[i] > l ? v[i] : l;
		}
		return l;
	}

	public static int sum(int[] v)
	{
		int s = 0;
		for (int i = 0; i < v.length; i++)
		{
			s += v[i];
		}
		return s;
	}

	public static int[] zeros(int n)
	{
		int[] z = new int[n];
		for (int i = 0; i < n; i++)
		{
			z[i] = 0;
		}
		return z;
	}

	public static int weight(int[] v) 
	{
		int wt = 0;
		for (int i = 0; i < v.length; i++) 
		{
			if (v[i] != 0)
			{
				wt++;
			}
		}
		return wt;
	}

	public static boolean isZero(int[] v)
	{
		for (int i = 0; i < v.length; i++)
		{
			if (v[i] != 0) return false;
		}
		return true;
	}

	public static boolean areEqual(int[] a, int[] b)
	{
		for (int i = 0; i < a.length; i++)
		{
			if (a[i] != b[i]) return false;
		}
		return true;
	}

	public static int[] add(int[] a, int[] b) throws Exception
	{
		int[] c = new int[a.length];
		if (a.length != b.length) throw new Exception("Vector dimensions don't match for addition.");
		for (int i = 0; i < a.length; i++)
		{
			c[i] = a[i] + b[i];
		}
		return c;
	}

	public static int[] addMod(int[] a, int[] b, int mod) throws Exception
	{
		int[] c = new int[a.length];
		if (a.length != b.length) throw new Exception("Vector dimensions don't match for addition.");
		for (int i = 0; i < a.length; i++)
		{
			c[i] = (a[i] + b[i]) % mod;
		}
		return c;
	}

	public static int compareReverse(int[] x, int[] y) throws Exception
	{
		if (x.length != y.length) throw new Exception("Error: vector dimensions do not match");
		for (int i = x.length - 1; i >= 0; i--) 
		{
			if (x[i] < y[i]) return 1;
			if (x[i] > y[i]) return -1;
		}
		return 0;
	}

	public static int compare(int[] x, int[] y) throws Exception
	{
		if (x.length != y.length) throw new Exception("Error: vector dimensions do not match");
		for (int i = x.length - 1; i >= 0; i--) 
		{
			if (x[i] < y[i]) return -1;
			if (x[i] > y[i]) return 1;
		}
		return 0;
	}

	public static int[] XOR(int[] x, int[] y) throws Exception
	{
		if (x.length != y.length) throw new Exception("Invalid x/y dimensions.");
		int[] and = new int[x.length];
		for (int i = 0; i < x.length; i++) 
		{
			and[i] = x[i] ^ y[i];
		}
		return and;
	}

	public static int[] AND(int[] x, int[] y) throws Exception
	{
		if (x.length != y.length) throw new Exception("Invalid x/y dimensions.");
		int[] and = new int[x.length];
		for (int i = 0; i < x.length; i++) 
		{
			and[i] = x[i] & y[i];
		}
		return and;
	}

	public static int[] NEGATE(int[] x) 
	{
		int[] nx = new int[x.length];
		for (int i = 0; i < x.length; i++) 
		{
			nx[i] = x[i] == 0 ? 1 : 0; // enforce only two possibilities
		}
		return nx;
	}

	public static ArrayList<MatrixState> paarOptimize(Matrix currMatrix, int n) throws Exception
	{
		ArrayList<MatrixState> history = new ArrayList<MatrixState>();
		int lastCol = n - 1;
		int hmax = 0;
		int wt = 0;
		int maxi = 0;
		int maxj = 0;
		int i = 0;
		int j = 0;

		// Column containers
		int[] coli;
		int[] colj;
		do
		{
			hmax = 0;
			for (i = 0; i <= lastCol - 1; i++)
			{
				for (j = i + 1; j <= lastCol; j++)
				{	
					coli = currMatrix.getColumn(i);
					colj = currMatrix.getColumn(j);
					wt = weight(AND(coli, colj));
					if (wt > hmax)
					{
						hmax = wt;
						maxi = i;
						maxj = j;
					}
				}
			}

			history.add(new MatrixState(currMatrix, maxi, maxj, lastCol));
			if (hmax > 1)
			{
				int[] maxcoli = currMatrix.getColumn(maxi);
				int[] maxcolj = currMatrix.getColumn(maxj);
				int[] newcol = AND(maxcoli, maxcolj);

				Matrix newMatrix = new Matrix(currMatrix);
				newMatrix.appendColumn(newcol);
				// Drop the variables from the other two entries...
				for (i = 0; i < currMatrix.getDimension(); i++)
				{
					if (newcol[i] == 1)
					{
						newMatrix.setEntry(i, maxi, 0);
						newMatrix.setEntry(i, maxj, 0);
					}
				}

				currMatrix = newMatrix;
				lastCol++;
				MatrixState newState = new MatrixState(currMatrix, maxi, maxj, lastCol);
			}
			//disp("Another go round: " + maxi + " " + maxj);
			//disp(currMatrix.toString());
		} while (hmax > 1);

		return history;
	}

	public static ArrayList<MatrixState> paarOptimizeExhaustive(Matrix currMatrix, int n, int oldi, int oldj, int lastcol) throws Exception
	{
		ArrayList<MatrixState> history = new ArrayList<MatrixState>();
		MatrixState state = new MatrixState(currMatrix, oldi, oldj, lastcol - 1);
		history.add(state);

		int hmax = 0;
		int wt = 0;
		int maxi = 0;
		int maxj = 0;
		int i = 0;
		int j = 0;
		int optimalGates = currMatrix.getGateCount();

		// Column containers
		int[] coli;
		int[] colj;

		hmax = 0;
		ArrayList<Pair> pairs = new ArrayList<Pair>();
		for (i = 0; i <= lastcol - 1; i++)
		{
			for (j = i + 1; j <= lastcol; j++)
			{	
				coli = currMatrix.getColumn(i);
				colj = currMatrix.getColumn(j);
				wt = weight(AND(coli, colj));
				if (wt > hmax)
				{
					pairs = new ArrayList<Pair>();
					pairs.add(new Pair(i, j));
					hmax = wt;
					maxi = i;
					maxj = j;
				} 
				else if (wt == hmax)
				{
					pairs.add(new Pair(i, j));
				}
			}
		}

		if (hmax > 1)
		{
			for (Pair p : pairs) 
			{
				MatrixState newState = new MatrixState(currMatrix, maxi, maxj, lastcol);
				history.add(newState);

				int[] maxcoli = currMatrix.getColumn(maxi);
				int[] maxcolj = currMatrix.getColumn(maxj);
				int[] newcol = AND(maxcoli, maxcolj);

				Matrix newMatrix = new Matrix(currMatrix);
				newMatrix.appendColumn(newcol);
				for (i = 0; i < newMatrix.getDimension(); i++)
				{
					if (newcol[i] == 1)
					{
						newMatrix.setEntry(i, maxi, 0);
						newMatrix.setEntry(i, maxj, 0);
					}
				}

				// Recurse
				ArrayList<MatrixState> rHistory = paarOptimizeExhaustive(newMatrix, i, j, n, lastcol + 1);
				if (rHistory.size() > 0)
				{
					int g = rHistory.get(rHistory.size() - 1).matrix.getGateCount();
					if (g <= optimalGates) 
					{
						optimalGates = g;
						history = rHistory;
						history.add(0, state);
					}
				}
			}
		}

		return history;
	}

	public static int crossCorrelation(int[] a, int[] b) throws Exception
	{
		return weight(AND(a, b));
	}

	public static ArrayList<Pair> maxCrossCorrelation(ArrayList<Pair> pairs, Matrix m) throws Exception
	{
		ArrayList<Pair> filtered = new ArrayList<Pair>();

		int maxCc = 0;
		for (Pair p : pairs)
		{
			int cc = crossCorrelation(m.getColumn(p.i), m.getColumn(p.j));
			if (cc > maxCc) maxCc = cc;
		}
		for (Pair p : pairs)
		{
			int cc = crossCorrelation(m.getColumn(p.i), m.getColumn(p.j));
			if (cc == maxCc) filtered.add(p);
		}

		return filtered;
	}

	public static ArrayList<Pair> maxExclusionCoefficient(ArrayList<Pair> pairs, Matrix m) throws Exception
	{
		ArrayList<Pair> filtered = new ArrayList<Pair>();

		int minCc = m.getDimension() * m.getDimension();
		for (int i = 0; i < pairs.size(); i++)
		{
			Pair pi = pairs.get(i);
			int cc = 0;
			for (int j = 0; j < pairs.size(); j++)
			{
				if (i != j)
				{
					Pair pj = pairs.get(j);
					int[] pairi = AND(m.getColumn(pi.i), m.getColumn(pi.j));
					int[] pairj = AND(m.getColumn(pj.i), m.getColumn(pj.j));
					int paircc = weight(AND(pairi, pairj));
					cc += paircc;
				}
			}
			if (cc < minCc) 
			{
				minCc = cc;
			}
		}

		for (int i = 0; i < pairs.size(); i++)
		{
			Pair pi = pairs.get(i);
			int cc = 0;
			for (int j = 0; j < pairs.size(); j++)
			{
				if (i != j)
				{
					Pair pj = pairs.get(j);
					int[] pairi = AND(m.getColumn(pi.i), m.getColumn(pi.j));
					int[] pairj = AND(m.getColumn(pj.i), m.getColumn(pj.j));
					int paircc = weight(AND(pairi, pairj));
					cc += paircc;
				}
			}

			if (cc == minCc) filtered.add(pi);
		}

		return filtered;
	}

	public static ArrayList<Pair> smallestPairWeight(ArrayList<Pair> pairs, Matrix m) throws Exception
	{
		ArrayList<Pair> filtered = new ArrayList<Pair>();

		int minWeight = m.getDimension() * m.getDimension();
		for (Pair p : pairs)
		{
			int sum = weight(m.getColumn(p.i)) + weight(m.getColumn(p.j));
			if (sum < minWeight) minWeight = sum;
		}
		for (Pair p : pairs)
		{
			int sum = weight(m.getColumn(p.i)) + weight(m.getColumn(p.j));
			if (sum == minWeight) filtered.add(p);
		}

		return filtered;
	}

	public static Pair nextLargestCorrelation(ArrayList<Pair> pairs, Matrix m, int lc) throws Exception
	{
		int maxCc = 0;
		Pair optimalPair = null;
		int n = m.getDimension();
		for (Pair p : pairs)
		{
			// Perform the factoring to look one step ahead in the process
			Matrix currMatrix = m.copy();
			int[] newCol = new int[n];
			int maxi = p.i;
			int maxj = p.j;
			for (int r = 0; r < n; r++)
			{
				int[] row = currMatrix.getRow(r);
				if (row[maxi] == 1 && row[maxj] == 1)
				{
					currMatrix.setEntry(r, maxi, 0);
					currMatrix.setEntry(r, maxj, 0);
					newCol[r] = 1;
				}
			}
			currMatrix.appendColumn(newCol);
			int lastcol = lc + 1;

			// Determine the largest coefficient from the next iteration (factoring)
			int maxCount = 0;
			ArrayList<Pair> newPairs = new ArrayList<Pair>();
			for (int i = 0; i <= lastcol - 1; i++)
			{
				for (int j = 0; j <= lastcol; j++)
				{
					if (i != j) // i and j are the columns, but really the variable IDs
					{
						int count = 0;
						for (int r = 0; r < n; r++)
						{
							int[] row = currMatrix.getRow(r);
							if (row[i] == 1 && row[j] == 1)
							{
								count++;
							}
						}
						if (count > maxCount)
						{
							newPairs = new ArrayList<Pair>();
							newPairs.add(new Pair(i, j));
							maxCount = count;
						}
					}
				}
			}

			// Compute the largest cross correlation and reset the optimal pair (this is the last heuristic in this example)
			ArrayList<Pair> ccList = maxCrossCorrelation(newPairs, currMatrix);
			int[] col1 = currMatrix.getColumn(ccList.get(0).i);
			int[] col2 = currMatrix.getColumn(ccList.get(0).j);
			int newCc = crossCorrelation(col1, col2);
			if (newCc > maxCc)
			{
				optimalPair = p;
			} 
		}

		return optimalPair;
	}

	public static ArrayList<MatrixState> cseOptimize(Matrix currMatrix, int n) throws Exception
	{
		ArrayList<MatrixState> history = new ArrayList<MatrixState>();
		int index = 0;
		int lastcol = n - 1;
		int maxCount = 0;
		int maxi = 0;
		int maxj = 0;

		do 
		{
			maxCount = maxi = maxj = 0;
			ArrayList<Pair> pairs = new ArrayList<Pair>();
			for (int i = 0; i <= lastcol - 1; i++)
			{
				for (int j = 0; j <= lastcol; j++)
				{
					if (i != j) // i and j are the columns, but really the variable IDs
					{
						int count = 0;
						for (int r = 0; r < n; r++)
						{
							int[] row = currMatrix.getRow(r);
							if (row[i] == 1 && row[j] == 1)
							{
								count++;
							}
						}
						if (count > maxCount)
						{
							pairs = new ArrayList<Pair>();
							pairs.add(new Pair(i, j));
							maxCount = count;
						}
						if (count == maxCount)
						{
							pairs.add(new Pair(i, j));
						}
					}
				}
			}

			// Apply the heuristic to break ties...
			if (pairs.size() == 1)
			{
				maxi = pairs.get(0).i;
				maxj = pairs.get(0).j;
			}
			else
			{
				ArrayList<Pair> step1 = maxCrossCorrelation(pairs, currMatrix);
				if (step1.size() == 0)
				{
					error("Step 1 failed.");
					System.exit(-1);
				}
				if (step1.size() == 1)
				{
					maxi = step1.get(0).i;
					maxj = step1.get(0).j;
				}
				else
				{
					ArrayList<Pair> step2 = maxExclusionCoefficient(step1, currMatrix);
					if (step2.size() == 0)
					{
						error("Step 2 failed.");
						System.exit(-1);
					}
					if (step2.size() == 1)
					{
						maxi = step2.get(0).i;
						maxj = step2.get(0).j;
					}
					else
					{
						ArrayList<Pair> step3 = smallestPairWeight(step2, currMatrix);
						if (step3.size() == 0)
						{
							error("Step 3 failed.");
							System.exit(-1);
						}
						if (step3.size() == 1)
						{
							maxi = step3.get(0).i;
							maxj = step3.get(0).j;
						}
						else
						{
							Pair optimal = nextLargestCorrelation(step3, currMatrix, lastcol);	
							maxi = optimal.i;
							maxj = optimal.j;
						}
					}
				}
			}

			if (maxCount > 1)
			{
				// replace xi + xj with xu, and add the new variable (a new column) to the set of equations
				MatrixState newState = new MatrixState(currMatrix, maxi, maxj, lastcol);
				history.add(newState);

				int[] newCol = new int[n];
				for (int r = 0; r < n; r++)
				{
					int[] row = currMatrix.getRow(r);
					if (row[maxi] == 1 && row[maxj] == 1)
					{
						currMatrix.setEntry(r, maxi, 0);
						currMatrix.setEntry(r, maxj, 0);
						newCol[r] = 1;
					}
				}
				currMatrix.appendColumn(newCol);
			}

			// Extend the matrix...
			lastcol++;
			disp(currMatrix.toString());
		} while (maxCount > 1);

		return history;
	}

	// static int combs; // number of linear combinations...
	// static SharedInteger d; // worst case...
	public static int peraltaDistance(final Matrix base, final int[] f) throws Exception
	{
		final int n = base.getDimension(); // number of rows - this is correct
		final int vars = f.length;
		int combs = 1;
		int d = vars;
		// d = new SharedInteger(vars);

		for (int i = 0; i <= n; i++) combs *= 2; // compute 2^b

		// ORIGINAL SEQUENTIAL CODE
		for (int c = 0; c < combs; c++) // we need at least one element from the basis...
		{
			int[] sum = zeros(vars);
			int adds = -1; 
			for (int i = 0; i < n; i++)
			{
				if (((1 << i) & c) != 0)
				{
					adds++;
					sum = addMod(sum, base.getRow(i), 2);
					// sum = add(sum, base.getRow(i));
				}
			}

			if (areEqual(f, sum) == true && adds < d)
			{
				d = adds;
			} 
		}
		return d;
		// return d.get();
	}

	public static int optimizedPeraltaDistance(Matrix base, int[] newBase, int[] f, int od) throws Exception
	{
		for (int i = 0; i < base.getDimension(); i++)
		{
			if (areEqual(base.getRow(i), f))
			{
				return 0;
			}
		}
		if (areEqual(newBase, f))
		{
			return 0;
		}

		if (reachable(base, XOR(newBase, f), od - 1, 0))
		{
			return od - 1; // distance can decrease by at most one since we're introducing a single new base
		}
		else
		{
			return od; // keep the old distance
		}
	}

	public static boolean reachable(Matrix base, int[] f, int targetDist, int index) throws Exception
	{
		if (base.getDimension() - index < targetDist) return false;
		if (targetDist == 0) return false;
		if (targetDist == 1)
		{
			for (int i = index; i < base.getDimension(); i++)
			{
				if (areEqual(f, base.getRow(i)))
				{
					return true;
				}
			}
			return false;
		}

		// Three cases from Peralta's implementation
		if (reachable(base, XOR(f, base.getRow(index)), targetDist - 1, index + 1)) 
		{
			return true;
		}
		if (reachable(base, f, targetDist, index + 1))
		{
			return true;
		}
		return false;
	}

	public static int[] computeDistance(Matrix base, Matrix m) throws Exception
	{
		int n = m.getDimension();
		int[] distance = new int[n];
		for (int i = 0; i < n; i++)
		{
			int[] fi = m.getRow(i);
			distance[i] = peraltaDistance(base, fi);
		}
		return distance;
	}

	public static int[] optimizedComputeDistance(Matrix base, Matrix m, int[] newBase, int[] dist) throws Exception
	{
		int n = m.getDimension();
		int[] distance = new int[n];
		for (int i = 0; i < n; i++)
		{
			int[] fi = m.getRow(i);
			distance[i] = optimizedPeraltaDistance(base, newBase, fi, dist[i]);
		}
		return distance;
	}

	public static BasePair pickBase(Matrix base, Matrix m, int[] dist, int tieBreaker) throws Exception
	{
		int n = base.getDimension(); // each row is a base...
		int d = sum(dist);
		int[] optimalBase = new int[base.getLength()]; // was dimension
		int[] optimalDist = new int[dist.length];
		ArrayList<int[]> optimalBases = new ArrayList<int[]>();
		ArrayList<int[]> optimalDistances = new ArrayList<int[]>();
		ArrayList<Pair> pairs = new ArrayList<Pair>();

		for (int i = 0; i < n; i++)
		{
			for (int j = i + 1; j < n; j++)
			{
				if (i != j)
				{
					int[] sum = addMod(base.getRow(i), base.getRow(j), 2);
					// Reduce modulo 2
					// int[] sum = add(base.getRow(i), base.getRow(j));
					// for (int k = 0; k < sum.length; k++) sum[k] = sum[k] % 2;

					if (!(base.containsRow(sum)) && !isZero(sum)) // bases can't be zero and we don't repeat bases
					{
						// Matrix newBase = base.copy();
						// newBase.appendRow(sum);
						int[] newDist = optimizedComputeDistance(base, m, sum, dist);
						int newDistSum = sum(newDist);
						if (newDistSum < d)
						{
							d = newDistSum;
							optimalBases = new ArrayList<int[]>();
							optimalDistances = new ArrayList<int[]>();
							pairs = new ArrayList<Pair>();
							optimalBases.add(sum);
							optimalDistances.add(newDist);
							pairs.add(new Pair(i, j));
						}
						else if (newDistSum == d)
						{
							optimalBases.add(sum);
							optimalDistances.add(newDist);
							pairs.add(new Pair(i, j));
						}
					}
				}
			}
		}

		// Resolve using norms...
		int mi = 0;
		int mj = 0;
		double maxNorm = Double.MAX_VALUE;
		switch (tieBreaker)
		{
			case 0: // norm
				maxNorm = euclideanNorm(optimalDistances.get(0));
				optimalBase = optimalBases.get(0);
				optimalDist = optimalDistances.get(0);
				mi = pairs.get(0).i;
				mj = pairs.get(0).j;
				for (int i = 1; i < optimalDistances.size(); i++)
				{
					double norm = euclideanNorm(optimalDistances.get(i));
					if (norm > maxNorm)
					{
						maxNorm = norm;
						optimalBase = optimalBases.get(i);
						optimalDist = optimalDistances.get(i);
						mi = pairs.get(i).i;
						mj = pairs.get(i).j;
					}
					else if (norm == maxNorm && comparePairs(pairs.get(i).i, pairs.get(i).j, mi, mj) <= 0)
					{
						maxNorm = norm;
						optimalBase = optimalBases.get(i);
						optimalDist = optimalDistances.get(i);
						mi = pairs.get(i).i;
						mj = pairs.get(i).j;	
					}
				}
				break;
			case 1: // norm-largest
				int largestDist = largest(optimalDistances.get(0));
				// for (int i = 0; i < optimalDistances.get(0).length; i++)
				// {
				// 	int tdist = optimalDistances.get(0)[i];
				// 	largestDist = tdist > largestDist ? tdist : largestDist;
				optimalBase = optimalBases.get(0);
				optimalDist = optimalDistances.get(0);
				mi = pairs.get(0).i;
				mj = pairs.get(0).j;
				// }
				maxNorm = Math.pow(euclideanNorm(optimalDistances.get(0)), 2);
				maxNorm = maxNorm - (double)largestDist;
				for (int i = 1; i < optimalDistances.size(); i++)
				{
					double norm = Math.pow(euclideanNorm(optimalDistances.get(i)), 2);
					norm = norm - (double)largest(optimalDistances.get(i));
					if (norm > maxNorm)
					{
						maxNorm = norm;
						optimalBase = optimalBases.get(i);
						optimalDist = optimalDistances.get(i);
						mi = pairs.get(i).i;
						mj = pairs.get(i).j;
					}
					else if (norm == maxNorm && comparePairs(pairs.get(i).i, pairs.get(i).j, mi, mj) <= 0)
					{
						maxNorm = norm;
						optimalBase = optimalBases.get(i);
						optimalDist = optimalDistances.get(i);
						mi = pairs.get(i).i;
						mj = pairs.get(i).j;	
					}
				}
				break;
			case 2: // norm-diff
				int firstLargestDist = 0;
				int secondLargestDist = -1;
				mi = pairs.get(0).i;
				mj = pairs.get(0).j;
				for (int i = 0; i < optimalDistances.get(0).length; i++)
				{
					int tdist = optimalDistances.get(0)[i];
					optimalBase = optimalBases.get(0);
					optimalDist = optimalDistances.get(0);
					if (tdist > firstLargestDist)
					{
						secondLargestDist = firstLargestDist;
						firstLargestDist = tdist;
					}
				}

				maxNorm = Math.pow(euclideanNorm(optimalDistances.get(0)), 2);
				maxNorm = maxNorm - (double)(firstLargestDist - secondLargestDist);
				for (int i = 1; i < optimalDistances.size(); i++)
				{
					double norm = Math.pow(euclideanNorm(optimalDistances.get(i)), 2);
					norm = norm - (double)(firstLargestDist - secondLargestDist);
					if (norm > maxNorm)
					{
						maxNorm = norm;
						optimalBase = optimalBases.get(i);
						optimalDist = optimalDistances.get(i);
						mi = pairs.get(i).i;
						mj = pairs.get(i).j;
					}
					else if (norm == maxNorm && comparePairs(pairs.get(i).i, pairs.get(i).j, mi, mj) <= 0)
					{
						maxNorm = norm;
						optimalBase = optimalBases.get(i);
						optimalDist = optimalDistances.get(i);
						mi = pairs.get(i).i;
						mj = pairs.get(i).j;	
					}
				}
				break;
			case 3: // random
				Random r = new Random(System.currentTimeMillis());
				int i = r.nextInt(optimalDistances.size());
				optimalBase = optimalBases.get(i);
				optimalDist = optimalDistances.get(i);
				mi = pairs.get(i).i;
				mj = pairs.get(i).j;
				break;
			default: // just pick the first one
				optimalBase = optimalBases.get(0); 
				optimalDist = optimalDistances.get(0);
				break;
		}

		//disp("picked " + mi + " " + mj + ": " + maxNorm);

		return new BasePair(optimalBase, optimalDist, new Pair(mi, mj));
	}

	public static SLP peraltaOptimize(Matrix m, int r, int c, int tieBreaker) throws Exception
	{
		ArrayList<String> slp = new ArrayList<String>();

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
		Matrix base = new Matrix(b, c);
		
		// Initialize the distance array and then get the ball rolling
		int[] dist = computeDistance(base, m);
		for (int f = 0; f < dist.length; f++)
		{
			if (dist[f] == 0)
			{
				int[] row = m.getRow(f);
				for (int cc = 0; cc < row.length; cc++)
				{
					if (row[cc] == 1)
					{
						slp.add("y_" + f + " = " + "x_" + cc);
					}
				}
			}
		}
		int i = m.getLength();
		int nVars = m.getLength();
		int xorCount = 0;
		while (isZero(dist) == false)
		{
			BasePair newBase = pickBase(base, m, dist, tieBreaker);
			base.appendRow(newBase.base);

			// Insert the new line
			String c1 = newBase.p.i < nVars ? "x" : "t";
			String c2 = newBase.p.j < nVars ? "x" : "t";
			slp.add("t_" + i + " = " + c1 + newBase.p.i + " XOR " + c2 + newBase.p.j);
			xorCount++;
			for (int f = 0; f < m.getDimension(); f++)
			{
				if (areEqual(newBase.base, m.getRow(f)))
				{
					slp.add("y_" + f + " = " + "t_" + i);
				}
			}
			i++;

			// Update the distance vector with this new base
			// dist = computeDistance(base, m);
			dist = newBase.dist;
		}

		return new SLP(slp, xorCount, 0); // no ANDs required
	}

	public static ArrayList<MatrixState> exhaustiveOptimize(Matrix currMatrix, int oldi, int oldj, int n, int lastcol) throws Exception
	{
		ArrayList<MatrixState> history = new ArrayList<MatrixState>();
		MatrixState state = new MatrixState(currMatrix, oldi, oldj, lastcol - 1);
		history.add(state);

		int[] coli;
		int[] colj;
		int wt = 0;
		int hmax = 0;
		int optimalGates = currMatrix.getGateCount();

		for (int i = 0; i <= lastcol - 1; i++) // lastcol - 1 = nm - 1 = n - 2, lastcol = nm
		{
			for (int j = i + 1; j <= lastcol; j++) // lastcol = n - 1
			{
				if (!(i < oldi && j != oldi && j != oldj && j < lastcol)) // prune
				{
					coli = currMatrix.getColumn(i);
					colj = currMatrix.getColumn(j);
					wt = weight(AND(coli, colj));
					Pair p = new Pair(i, j);
					if (wt > 1) // if the weight of the AND is > 1, then 
					{
						int[] newcol = AND(coli, colj);
						// disp(newcol);
						Matrix newMatrix = new Matrix(currMatrix);
						newMatrix.appendColumn(newcol);
						for (int k = 0; k < n; k++)
						{
							if (newcol[k] == 1)
							{
								newMatrix.setEntry(k, i, 0);
								newMatrix.setEntry(k, j, 0);
							}
						}

						// Recurse and only save the output if it's better than what we have so far
						ArrayList<MatrixState> rHistory = exhaustiveOptimize(newMatrix, i, j, n, lastcol + 1);
						if (rHistory.size() > 0)
						{
							int g = rHistory.get(rHistory.size() - 1).matrix.getGateCount();
							if (g < optimalGates)
							{
								optimalGates = g;
								history = rHistory;
								history.add(0, state);
							}
						}
					}
				}
			}
		}

		return history;
	}

	public static int comparePairs(int i1, int j1, int i2, int j2)
	{
		if (i1 < i2)
		{
			return -1;
		}
		else if (i1 > i2)
		{ 
			return 1;
		}
		else
		{
			if (j1 < j2)
			{
				return -1;
			}
			else if (j1 > j2)
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}
	}

	public static ArrayList<Vector> sortRowsReverseOrder(ArrayList<Vector> rows) throws Exception
	{
		boolean swapped = true;
		while (swapped) // simple bubble sort to put in reverse order
		{
			swapped = false;
			for (int i = 0; i < rows.size(); i++)
			{
				for (int j = i + 1; j < rows.size(); j++)
				{
					Vector vi = rows.get(i);
					Vector vj = rows.get(j);
					int[] ri = vi.row;;
					int[] rj = vj.row;
					if (compareReverse(ri, rj) == -1) // then swap
					{
						rows.remove(i);
						rows.add(i, vj);
						rows.remove(j);
						rows.add(j, vi);
						swapped = true;
					}
				}
			}
		}
		return rows;	
	}

	// this builds an SLP...
	public static ArrayList<String> bernsteinOptimize(ArrayList<Vector> rows, int n, int p, int q) throws Exception
	{
		ArrayList<String> slp = new ArrayList<String>();
		rows = sortRowsReverseOrder(rows); // sort in reverse order - this is redundant at times

		if (q == 0) return slp;
		if (p == 0)
		{
			return slp;
		}
		if (rows.get(q - 1).row[p - 1] == 0) // x_{p-1} is unused, so don't consider next time around
		{
			slp.addAll(bernsteinOptimize(rows, n, p - 1, q));
			return slp;
		}
		else if (q >= 2) // q - 1 == j
		{
			if (rows.get(q - 2).row[p - 1] == 1) // q - 2 == i
			{
				int ti = rows.get(q - 2).id;
				int tj = rows.get(q - 1).id; // newRowTag[q - 1];
				int[] xor = XOR(rows.get(q - 1).row, rows.get(q - 2).row);
				rows.remove(q - 1);
				rows.add(q - 1, new Vector(xor, tj));
				slp.add("y_" + tj + " = y_" + ti + " XOR y_" + tj);
				slp.addAll(bernsteinOptimize(rows, n, p, q));
				return slp;
			}
			else
			{
				rows.get(q - 1).row[p - 1] = 0; // set the row to 0...
				slp.add("y_" + rows.get(q - 1).id + " = x_" + (p - 1) + " XOR y_" + rows.get(q - 1).id);
				slp.addAll(bernsteinOptimize(rows, n, p, q));
				return slp;
			}
		}

		return slp;
	}

	public static int numGates(int n, ArrayList<MatrixState> history) throws Exception
	{
		if (history.size() == 0)
		{
			throw new Exception("Invalid history.");
		}
		return history.get(history.size() - 1).matrix.getGateCount();
	}

	public static ArrayList<String> buildSLP(ArrayList<MatrixState> history, int vars) throws Exception
	{
		ArrayList<String> slp = new ArrayList<String>();

		Matrix full = history.get(history.size() - 1).matrix;
		int cutoff = full.getDimension();
		// int vars = cutoff + history.size() - 1;
		Matrix opt = full.subMatrix(cutoff);

		// Compute the new variables...
		for (int i = history.size() - 1; i > 0; i--)
		{
			int ii = history.get(i).i;
			int jj = history.get(i).j;
			String c1 = ii >= cutoff ? "t" : "x";
			String c2 = jj >= cutoff ? "t" : "x";
			String line = "t_" + (cutoff + i - 1) + " = " + c1 + "_" + ii + " XOR " + c2 + "_" + jj;
			slp.add(line);
		}

		// Now that we have the values of each and every one of these thingys, compute the output for each bit
		for (int i = 0; i < cutoff; i++)
		{
			int[] row = full.getRow(i);
			int rw = weight(row);
			if (rw == 0)
			{
				String line = "y_" + i + " = 0";
				slp.add(line);
			}
			else if (rw == 1)
			{
				for (int j = 0; j < vars; j++)
				{
					if (row[j] == 1)
					{
						String c = j < cutoff ? "x" : "t";
						String line = "y_" + i + " = " + c + "_" + j;
						slp.add(line);
					}
				}
			}
			else // at least two variables in this linear form
			{
				int op1 = -1;
				int op2 = -1;
				for (int j = 0; j < vars; j++)
				{
					if (row[j] == 1 && op2 == -1) // setting up first operand
					{
						op2 = j;
					}
					else if (row[j] == 1 && op1 == -1) // first time with two operands
					{
						op1 = op2;
						op2 = j;
						String c1 = op1 >= cutoff ? "t" : "x";
						String c2 = op2 >= cutoff ? "t" : "x";
						String line = "y_" + i + " = " + c1 + "_" + op1 + " XOR " + c2 + "_" + op2;
						slp.add(line);
					}
					else if (row[j] == 1)
					{
						op1 = op2;
						op2 = j;
						String c2 = op2 >= cutoff ? "t" : "x";
						String line = "y_" + i + " = y_" + i + " XOR " + c2 + "_" + op2;
						slp.add(line);
					}
				}	
			}
		}

		return slp;	
	}

	public static void test() throws Exception
	{
		// disp("STARTING PAAR TEST");
		// // int[][] test = {{1,1,0,1,1,1,0,0},{0,1,1,1,1,1,1,0},{0,1,0,1,1,1,0,0},{1,1,0,1,1,0,1,0},{1,0,0,1,0,0,1,0},{1,1,0,1,1,0,1,0},{1,1,0,1,1,0,0,0},{1,1,1,1,1,1,1,1}};
		// int[][] test = {{1,0,0,0,1},{0,1,1,0,1},{0,1,0,1,0},{0,1,1,0,1},{0,0,1,1,1}};
		// Matrix m = new Matrix(test, 5);
		// disp(m.toString());
		// ArrayList<MatrixState> history = paarOptimize(m, 5); // 5x5 matrix at this point...
		// disp("" + history.size());
		// disp(history.get(history.size() - 1).matrix.toString());
		// disp("Number of gates: " + numGates(5, history));

		// disp("STARTING CSE TEST");
		// int[][] test2 = {{1,1,0,1,1,1,0,0},{0,1,1,1,1,1,1,0},{0,1,0,1,1,1,0,0},{1,1,0,1,1,0,1,0},{1,0,0,1,0,0,1,0},{1,1,0,1,1,0,1,0},{1,1,0,1,1,0,0,0},{1,1,1,1,1,1,1,1}};
		// Matrix m2 = new Matrix(test2, 8);
		// disp(m2.toString());
		// ArrayList<MatrixState> history2 = cseOptimize(m2, 8);
		// disp("Number of gates: " + numGates(8, history2));

		// disp("STARTING EXHAUSTIVE TEST");
		// // int[][] test3 = {{1,0,0,0,1},{0,1,1,0,1},{0,1,0,1,0},{0,1,1,0,1},{0,0,1,1,1}};
		// int[][] test3 = {{1,1,0,1,1,1,0,0},{0,1,1,1,1,1,1,0},{0,1,0,1,1,1,0,0},{1,1,0,1,0,0,1,0},{1,0,0,1,1,0,1,0},{1,1,0,1,1,0,1,0},{1,1,0,1,0,0,0,0},{1,1,1,1,1,1,1,1}};
		// Matrix m3 = new Matrix(test3, 8);
		// disp(m3.toString());
		// long start = System.currentTimeMillis();
		// ArrayList<MatrixState> history3 = exhaustiveOptimize(m3, 0, 0, 8, 7); // 5x5 matrix at this point...
		// long end = System.currentTimeMillis();
		// // ArrayList<MatrixState> history3 = exhaustiveOptimize(m3, 0, 0, 5, 4); // 5x5 matrix at this point...
		// disp("" + history3.size());
		// disp("" + history3.get(history3.size() - 1).matrix.getGateCount());
		// disp(history3.get(history3.size() - 1).matrix.toString());
		// disp("" + (end - start) + " ms");
		// buildSLP(history3);

		// disp("STARTING EXHAUSTIVE TEST");
		// int[][] test3 = {{1,1,1,1},{1,0,1,1},{1,1,1,0},{1,0,0,1}};
		// Matrix m3 = new Matrix(test3, 4);
		// disp(m3.toString());
		// ArrayList<MatrixState> history3 = exhaustiveOptimize(m3, 0, 0, 4, 3, new ArrayList<Pair>()); // 5x5 matrix at this point...
		// disp("" + history3.size());
		// disp(history3.get(history3.size() - 1).matrix.toString());
		// for (int i = 0; i < history3.size(); i++)
		// {
		// 	disp("Pair: " + history3.get(i).i + " " + history3.get(i).j);
		// }
		// disp("Number of gates: " + numGates(4, history3));

		disp("STARTING PERALTA TEST");
		int[][] test4 = {{1,1,1,0,0},{0,1,0,1,1},{1,0,1,1,1},{0,1,1,1,0},{1,1,0,1,0},{0,1,1,1,1}};
		Matrix m4 = new Matrix(test4, 5);
		disp(m4.toString());
		SLP slp4 = peraltaOptimize(m4, 6, 5, 0);
		// SLP slp4 = parallelPeraltaOptimize(m4, 6, 5, 0);
		dispStrings(slp4.lines);
		// long start = System.currentTimeMillis();
		// SLP slp5 = parallelPeraltaOptimize(m4, 6, 5, 0);
		// long end = System.currentTimeMillis();
		// disp("" + (end - start));
		// dispStrings(slp5.lines);
		// start = System.currentTimeMillis();
		// SLP slp6 = parallelPeraltaOptimize_v2(m4, 6, 5, 0);
		// end = System.currentTimeMillis();
		// disp("" + (end - start));
		// dispStrings(slp6.lines);

		// disp("Number of gates: " + numGates(5, history4));

		// disp("STARTING BERNSTEIN TEST");
		// int[][] test5 = {{1,1,1,1},{1,0,1,1},{1,1,1,0},{1,0,0,1}}; // tags make up the last column, for easy purposes
		// // int[][] test5 = {{1,0,1},{1,1,1},{0,0,1}};
		// // Matrix m5 = new Matrix(test5, 4);
		// Matrix m5 = new Matrix(test5, 3);
		// disp(m5.toString());
		// ArrayList<String> slp = bernsteinOptimize(m5.getVectors(), 4, 4, 4);
		// // ArrayList<MatrixState> history5 = bernsteinOptimize(m5.getVectors(), 3, 3, 3);
	}

	public static ArrayList<Matrix> buildMatrices(String file) throws Exception
	{
		ArrayList<Matrix> matrices = new ArrayList<Matrix>();
		BufferedReader in = new BufferedReader(new FileReader(file)); 
		String line = "";
		while ((line = in.readLine()) != null && !line.isEmpty())
		{
			String[] data = line.split(" ");
			int p = Integer.parseInt(data[0]);
			int q = Integer.parseInt(data[1]);
			line = in.readLine();
			
			ArrayList<int[]> rows = new ArrayList<int[]>();
			int index = 0;
			for (int i = 0; i < p; i++)
			{
				int[] row = new int[q];
				for (int j = 0; j < q; j++)
				{
					row[j] = Integer.parseInt("" + line.charAt(index++));
				}
				rows.add(row);
			}
			Matrix m = new Matrix(rows, q); // include the initial number of columns (variables)
			matrices.add(m);
		}
		return matrices;
	}

	public static int optimize(ArrayList<Matrix> matrices, boolean debug) throws Exception
	{
		int total = 0;
		for (Matrix m : matrices)
		{
			int gc = Integer.MAX_VALUE;
			int tgc = 0;

			error(m.toString());
			error("STARTING PAAR TEST");
			ArrayList<MatrixState> history = paarOptimize(m, m.getLength()); 
			tgc = numGates(m.getDimension(), history);
			gc = tgc < gc ? tgc : gc;
			// disp("Number of gates: " + numGates(m.getDimension(), history));
			if (debug) dispStrings(buildSLP(history, m.getLength()));

			// disp("STARTING PAAR PSEUDO-EXHAUSTIVE TEST");
			// ArrayList<MatrixState> history1 = paarOptimizeExhaustive(m, m.getLength(), 0, 0, m.getLength() - 1); // 5x5 matrix at this point...
			// disp("Number of gates: " + numGates(m.getDimension(), history1));
			// dispStrings(buildSLP(history1));

			// disp("STARTING EXHAUSTIVE FACTORIZATION TEST");
			// ArrayList<MatrixState> history3 = exhaustiveOptimize(m, 0, 0, m.getLength(), m.getLength() - 1); // 5x5 matrix at this point...
			// disp("Number of gates: " + numGates(m.getDimension(), history3));
			// dispStrings(buildSLP(history3));

			error("STARTING PERALTA TEST - TIE 0");
			SLP slp4 = peraltaOptimize(m, m.getDimension(), m.getLength(), 0);
			tgc = slp4.xc;
			gc = tgc < gc ? tgc : gc;
			// error("Number of gates: " + slp4.xc);
			if (debug) dispStrings(slp4.lines);

			// error("STARTING PERALTA TEST - TIE 1");
			// slp4 = parallelPeraltaOptimize_v2(m, m.getDimension(), m.getLength(), 1);
			// tgc = slp4.xc;
			// gc = tgc < gc ? tgc : gc;
			// // error("Number of gates: " + slp4.xc);
			// if (debug) dispStrings(slp4.lines);

			// error("STARTING PERALTA TEST - TIE 2");
			// slp4 = parallelPeraltaOptimize_v2(m, m.getDimension(), m.getLength(), 2);
			// tgc = slp4.xc;
			// gc = tgc < gc ? tgc : gc;
			// // error("Number of gates: " + slp4.xc);
			// if (debug) dispStrings(slp4.lines);

			// error("STARTING PERALTA TEST - TIE 3");
			// slp4 = parallelPeraltaOptimize_v2(m, m.getDimension(), m.getLength(), 3);
			// tgc = slp4.xc;
			// gc = tgc < gc ? tgc : gc;
			// // disp("Number of gates: " + slp4.xc);
			// if (debug) dispStrings(slp4.lines);

			// error("STARTING BERNSTEIN TEST");
			// ArrayList<String> slp5 = bernsteinOptimize(m.getVectors(), m.getLength(), m.getLength(), m.getDimension());
			// tgc = slp5.size();
			// gc = tgc < gc ? tgc : gc;
			// // disp("Number of gates: " + slp5.size());
			// if (debug) dispStrings(slp5);

			total += gc;
		}

		return total;
	}

	public static void main(String[] args) throws Exception
	{
		if (args.length == 0)
		{
			test(); // just run the hard-coded test cases
		}
		else if (args.length == 1)
		{
			ArrayList<Matrix> matrices = buildMatrices(args[0]);
			disp(optimize(matrices, true)); 
		}
		else if (args.length == 3)
		{
			int lb = Integer.parseInt(args[0]);
			int ub = Integer.parseInt(args[1]);
			String header = args[2];

			long start = System.currentTimeMillis();
			for (int i = lb; i < ub; i++)
			{
				// Forward
				String fname = header + "_" + i + "_forward.txt";
				ArrayList<Matrix> matrices = buildMatrices(fname);
				int ftotal = optimize(matrices, false);

				// Inverse
				fname = header + "_" + i + "_inverse.txt";
				matrices = buildMatrices(fname);
				int itotal = optimize(matrices, false);

				// Merged
				fname = header + "_" + i + "_merged.txt";
				matrices = buildMatrices(fname);
				int mtotal = optimize(matrices, false);

				// Write the results to the file header_i_count.txt
				// First line contains the forward and merged counts, separately, and then combined
				// Second line contains merged count
				disp("" + ftotal + " " + itotal + " " + (ftotal + itotal));
				disp("" + mtotal);
			}
			long end = System.currentTimeMillis();
			error("Time: " + (end - start) + "ms");
			error("Memory: " + ((double)(Runtime.getRuntime().freeMemory() / (double)Runtime.getRuntime().maxMemory())) + ", " + Runtime.getRuntime().totalMemory());
		}
		else
		{	
			error("Usage: java MatrixOptimize [file]");
			error("   - [0] no argument will cause the hard-coded test cases to run");
			error("   - [1] file contains a list of file names that contain matrices to optimize");
			error("   - [2] lower bound, upper bound, and header for all files - will write a _count file for each index");
		}
	}
}

