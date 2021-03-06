package caw.circuits.optimize.linear;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Deque;
import java.util.ArrayDeque;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;

public class MatrixOptimize
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
	static void disp_strings(ArrayList<String> vs) {
		for (int i = 0; i < vs.size(); i++) disp(vs.get(i));
	}
	static void disp_vectors(ArrayList<Vector> vs) {
		for (int i = 0; i < vs.size(); i++) disp(vs.get(i).row);
	}
	static ArrayList<Pair> copy(ArrayList<Pair> pairs) {
		ArrayList<Pair> np = new ArrayList<Pair>();
		for (int i = 0; i < pairs.size(); i++) np.add(pairs.get(i));
		return np;
	}

	public static double taxi_norm(int[] v)
	{
		return (double)sum(v);
	}

	public static double euclidean_norm(int[] v)
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

	public static boolean is_zero(int[] v)
	{
		for (int i = 0; i < v.length; i++)
		{
			if (v[i] != 0) return false;
		}
		return true;
	}

	public static boolean are_equal(int[] a, int[] b)
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

	public static int[] add_mod(int[] a, int[] b, int mod) throws Exception
	{
		int[] c = new int[a.length];
		if (a.length != b.length) throw new Exception("Vector dimensions don't match for addition.");
		for (int i = 0; i < a.length; i++)
		{
			c[i] = (a[i] + b[i]) % mod;
		}
		return c;
	}

	public static int compare_reverse(int[] x, int[] y) throws Exception
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

	public static int[] FLIP(int[] x, int i) throws Exception
	{
		if (i < 0 || i >= x.length) throw new Exception("Invalid index into the vector.");
		x[i] = x[i] == 1 ? 0 : 1; // flip
		return x;
	}

	public static int[] CLONE(int[] x)
	{
		int[] nx = new int[x.length];
		for (int i = 0; i < x.length; i++) nx[i] = x[i];
		return nx;
	}

	public static ArrayList<MatrixState> paar_optimize(Matrix currMatrix, int n) throws Exception
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

	public static ArrayList<MatrixState> paar_optimize_exhaustive(Matrix currMatrix, int n, int oldi, int oldj, int lastcol) throws Exception
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
				ArrayList<MatrixState> rHistory = paar_optimize_exhaustive(newMatrix, i, j, n, lastcol + 1);
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

	// static int combs; // number of linear combinations...
	// static SharedInteger d; // worst case...
	public static int peralta_distance(final Matrix base, final int[] f) throws Exception
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
					sum = add_mod(sum, base.getRow(i), 2);
					// sum = add(sum, base.getRow(i));
				}
			}

			if (are_equal(f, sum) == true && adds < d)
			{
				d = adds;
			} 
		}
		return d;
		// return d.get();
	}

	public static ArrayList<Vector> walk_weighted_sequences(int n, int k) throws Exception
	{
		ArrayList<Vector> visited = new ArrayList<Vector>();
		ArrayList<int[]> queue = new ArrayList<int[]>();
	
		// Build the max node (root)
		int[] max = new int[n];
		for (int i = 0; i < k; i++) max[i] = 1;
		queue.add(max);

		// Walk...
		while (queue.size() > 0)
		{
			int[] seq = queue.remove(0); // pop 
			ArrayList<Integer> nz = find_nonzeros(seq);
			ArrayList<Integer> z = find_zeros(seq);
			for (int i = 0; i < nz.size(); i++)
			{
				for (int j = 0; j < z.size(); j++)
				{
					int[] mod = CLONE(seq);
					mod = FLIP(FLIP(mod, nz.get(i)), z.get(j)); 
					Vector v = new Vector(mod);
					if (!(visited.contains(v)))
					{
						visited.add(v);
						queue.add(mod);
					}
				}
			}
		}
		return visited;
	}

// void HamiltonianPath(char direction, int n, int H, string pre, string post)
	public enum WALK_DIR {LEFT, RIGHT};
	public static String lastString;
	public static HashMap<String, int[]> xorMap = new HashMap<String, int[]>();

	// public static int[] APPEND(int[] l, int[] r)
	// {
	// 	int ii = 0;
	// 	int[] x = new int[l.length + r.length];
	// 	for (int i = 0; i < l.length; i++) x[ii] = l[ii++];
	// 	for (int i = 0; i < r.length; i++) x[ii] = r[ii++];
	// 	return x;
	// }

	// public static int[] APPEND(int[] l, int[] r, int bit)
	// {
	// 	int ii = 0;
	// 	int[] x = new int[l.length + r.length + 1];
	// 	for (int i = 0; i < l.length; i++) x[ii] = l[ii++];
	// 	x[ii++] = bit;
	// 	for (int i = 0; i < r.length; i++) x[ii] = r[ii++];
	// 	return x;
	// }

	public static ArrayList<Vector> walk_weighted_sequences_v2(WALK_DIR dir, int n, int k, String pre, String post) throws Exception
	{
		ArrayList<Vector> seqs = new ArrayList<Vector>();
		if (n == 1)
		{
			String bit = k == 0 ? "0" : "1";
			String newSeq = pre + bit + post;
			System.out.println("seq: " + newSeq);
			seqs.add(new Vector(newSeq));
			lastString = newSeq;
			return seqs;
		}
		else if (k == 0)
		{
			String tmp = pre;
			for (int i = 0; i < n; i++) tmp = tmp + "0";
			tmp = tmp + post;
			System.out.println("seq: " + tmp);
			seqs.add(new Vector(tmp));
			lastString = tmp;
			return seqs;
		}
		else if (k == n)
		{
			String tmp = pre;
			for (int i = 0; i < n; i++) tmp = tmp + "1";
			tmp = tmp + post;
			System.out.println("seq: " + tmp);
			seqs.add(new Vector(tmp));
			lastString = tmp;
			return seqs;
		}
		else if (dir == WALK_DIR.LEFT)
		{
			System.out.println("left to right: " + n + " : " + k + " " + pre + " 1 " + post);
			seqs.addAll(walk_weighted_sequences_v2(WALK_DIR.RIGHT, n - 1, k - 1, pre + "1", post));
			System.out.println("left to left: " + n + " : " + k + " " + pre + " 0 " + post);
			seqs.addAll(walk_weighted_sequences_v2(WALK_DIR.LEFT, n - 1, k, pre + "0", post));
			return seqs;
		}
		else
		{
			System.out.println("right to left: " + n + " : " + k + " " + pre + " 1 " + post);
			seqs.addAll(walk_weighted_sequences_v2(WALK_DIR.LEFT, n - 1, k - 1, pre + "1", post));
			System.out.println("right to right: " + n + " : " + k + " " + pre + " 0 " + post);
			seqs.addAll(walk_weighted_sequences_v2(WALK_DIR.RIGHT, n - 1, k, pre + "0", post));
			return seqs;
		}
	}

	// To avoid throwing stuff on the stack
	public static int[] walk_ft;
	public static int[] walk_newBase;
	public static int walk_od;
	public static Matrix walk_base;

	public static int walk_weighted_sequences_and_change_with_graph_v2(WALK_DIR dir, int n, 
		int k, int bit, int[] acc) throws Exception
	{
		if (n == 1)
		{
			int[] target = acc;
			if (k == 1)
			{	
				target = XOR(walk_base.getRow(bit), acc);
			}
			if (are_equal(XOR(walk_newBase, target), walk_ft))
			{
				return walk_od;
			}
			else
			{
				return walk_od + 1;	
			}
		}
		else if (k == 0)
		{
			if (are_equal(XOR(walk_newBase, acc), walk_ft))
			{
				return walk_od;
			}
			else
			{
				return walk_od + 1;
			}
		}
		else if (k == n)
		{
			int[] target = acc;
			for (int i = bit; i < (n + bit); i++) 
			{
				target = XOR(target, walk_base.getRow(i));
			}
			if (are_equal(XOR(walk_newBase, target), walk_ft))
			{
				return walk_od;
			}
			else
			{
				return walk_od + 1;
			}
		}
		else if (dir == WALK_DIR.LEFT)
		{
			if (walk_weighted_sequences_and_change_with_graph_v2(WALK_DIR.RIGHT, n - 1, k - 1, bit + 1, 
				XOR(acc, walk_base.getRow(bit))) == walk_od)
			{
				return walk_od;
			}
			if (walk_weighted_sequences_and_change_with_graph_v2(WALK_DIR.LEFT, n - 1, k, bit + 1, acc) == walk_od)
			{
				return walk_od;
			}
			return walk_od + 1;
		}
		else
		{
			if (walk_weighted_sequences_and_change_with_graph_v2(WALK_DIR.LEFT, n - 1, k - 1, bit + 1, 
				XOR(walk_base.getRow(bit), acc)) == walk_od)
			{
				return walk_od;
			}
			if (walk_weighted_sequences_and_change_with_graph_v2(WALK_DIR.RIGHT, n - 1, k, bit + 1, acc) == walk_od)
			{
				return walk_od;
			}
			return walk_od + 1;
		}
	}

	// To avoid throwing stuff on the stack
	public static void kweight_sums(ArrayList<int[]> sums, WALK_DIR dir, int n, int k, int bit, int[] acc) throws Exception
	{

		// NOTE: sums assumed to be non-null from caller

		if (n == 1)
		{
			int[] target = acc;
			if (k == 1)
			{	
				target = XOR(walk_base.getRow(bit), acc);
			}
			sums.add(XOR(walk_newBase, target));
			// if (are_equal(XOR(walk_newBase, target), walk_ft))
			// {
			// 	return walk_od;
			// }
			// else
			// {
			// 	return walk_od + 1;	
			// }
		}
		else if (k == 0)
		{
			sums.add(XOR(walk_newBase, acc));
			// if (are_equal(XOR(walk_newBase, acc), walk_ft))
			// {
			// 	return walk_od;
			// }
			// else
			// {
			// 	return walk_od + 1;
			// }
		}
		else if (k == n)
		{
			int[] target = acc;
			for (int i = bit; i < (n + bit); i++) 
			{
				target = XOR(target, walk_base.getRow(i));
			}
			sums.add(XOR(walk_newBase, target));
			// if (are_equal(XOR(walk_newBase, target), walk_ft))
			// {
			// 	return walk_od;
			// }
			// else
			// {
			// 	return walk_od + 1;
			// }
		}
		else if (dir == WALK_DIR.LEFT)
		{
			kweight_sums(sums, WALK_DIR.RIGHT, n - 1, k - 1, bit + 1, XOR(acc, walk_base.getRow(bit)));
			kweight_sums(sums, WALK_DIR.LEFT, n - 1, k, bit + 1, acc);
		}
		else
		{
			kweight_sums(sums, WALK_DIR.LEFT, n - 1, k - 1, bit + 1, XOR(acc, walk_base.getRow(bit)));
			kweight_sums(sums, WALK_DIR.RIGHT, n - 1, k, bit + 1, acc);
		}
	}

	// compute the sums of all k-combinations of rows between start and end index of the base rows
	// that is... pick k base signals in [start,end] and compute their sums
	public static ArrayList<int[]> hamiltonian_sums(Matrix base, int[] newBase, int start, int end, int k) throws Exception
	{
		ArrayList<int[]> sums = new ArrayList<int[]>(); // end is inclusive

		// disp("start/end = " + start + "," + end);
		// disp(base.toString());
		walk_base = base.subMatrix(start, end + 1, 0, base.getLength());
		int[] acc = new int[base.getLength()];
		for (int j = 0; j < acc.length; j++) acc[j] = 0;
		walk_newBase = newBase;
		// distance[i] = walk_weighted_sequences_and_change_with_graph_v2(WALK_DIR.RIGHT, base.getDimension(), dist[i] - 1, 0, ft1);
		kweight_sums(sums, WALK_DIR.RIGHT, end - start + 1, k, 0, acc); // n = start - end + 1 
		return sums;
	}

	public static int walk_weighted_sequencesAndCheck(int n, int k, Matrix base, int[] f, int[] newBase) throws Exception
	{
		//ArrayList<Integer> nz = find_nonzeros(seq);
		//ArrayList<Integer> z = find_zeros(seq);
		//ArrayList<Integer> seqs = new ArrayList<Integer>();
		ArrayList<Vector> visited = new ArrayList<Vector>();
		ArrayList<int[]> queue = new ArrayList<int[]>();
		ArrayList<Integer> starts = new ArrayList<Integer>();
		ArrayList<Integer> ends = new ArrayList<Integer>();
	
		// Build the max node (root)
		int[] max = new int[n];
		for (int i = 0; i < k; i++) max[i] = 1;
		queue.add(max);

		// Walk the tree...
		boolean firstRun = true;
		while (queue.size() > 0)
		{
			int[] seq = queue.remove(0); // pop out the first one...
			int start = firstRun ? -1 : starts.remove(0);
			int end = firstRun ? -1 : ends.remove(0);
			for (int i = 0; i < n; i++)
			{
				if (seq[i] == 1 && (i > start || start == -1))
				{
					for (int j = 0; j < n; j++)
					{
						if (i != j && seq[j] == 0 && (j > end || end == -1))
						{
							// int[] mod = CLONE(seq);
							// mod = FLIP(FLIP(mod, i), j); // flip both spots (keeps the same weight but produces a new vector)
							Vector v = new Vector(seq);
							v.row[i] = v.row[i] == 0 ? 1 : 0;
							v.row[j] = v.row[j] == 0 ? 1 : 0;
							
							if (!(visited.contains(v))) // this won't happen with pruning... right?
							{
								// firstRun = false;
								starts.add(i);
								ends.add(j);
								visited.add(v);
								queue.add(v.row);

								// Check to see if we have a match...
								if (are_equal(XOR(newBase, base.xorRows(v.row)), f))
								{
									return k;
								}
							}
						}
					}
				}
			}
		}

		return k + 1;
	}

	public static int walk_weighted_sequences_and_check_v2(int n, int k, Matrix base, int[] f, int[] newBase) throws Exception
	{
		LinkedList<Integer> queue = new LinkedList<Integer>();
		HashSet<Integer> visited = new HashSet<Integer>();

		String ms = "";
		for (int i = 0; i < k; i++) ms = ms + "1";
		for (int i = k; i < n; i++) ms = ms + "0";
		queue.add(Integer.parseInt(ms, 2));

		while (queue.size() > 0)
		{
			int seq = queue.removeFirst(); // pop
			for (int i = 0; i < n; i++)
			{ 
				if ((seq & (1 << i)) > 0) // ones
				{
					for (int j = 0; j < n; j++)
					{
						if ((seq & (1 << j)) == 0) // zeros
						{
							// int mod = seq.flipBit(i).flipBit(j);
							// if (!(visited.contains(mod)))
							int mod = ((seq ^ (1 << i)) ^ (1 << j));
							if (!(visited.contains(mod)))
							{
								visited.add(mod);
								queue.addLast(mod);

								// Check to see if we have a match...
								if (are_equal(XOR(newBase, base.xorRows(mod)), f))
								{
									return k;
								}
							}
						}
					}
				}
			}
		}

		return k + 1;
	}

	public static int walk_weighted_sequences_and_check_v3(int n, int k, Matrix base, int[] f, int[] newBase) throws Exception
	{
		LinkedList<Integer> queue = new LinkedList<Integer>();
		HashSet<Integer> visited = new HashSet<Integer>();

		String ms = "";
		for (int i = 0; i < k; i++) ms = ms + "1";
		for (int i = k; i < n; i++) ms = ms + "0";
		queue.add(Integer.parseInt(ms, 2));

		while (queue.size() > 0)
		{
			int seq = queue.removeFirst(); // pop
			for (int i = 0; i < n; i++)
			{ 
				if ((seq & (1 << i)) > 0) // ones
				{
					for (int j = 0; j < n; j++)
					{
						if ((seq & (1 << j)) == 0) // zeros
						{
							// int mod = seq.flipBit(i).flipBit(j);
							// if (!(visited.contains(mod)))
							int mod = ((seq ^ (1 << i)) ^ (1 << j));
							if (!(visited.contains(mod)))
							{
								visited.add(mod);
								queue.addLast(mod);

								// Check to see if we have a match...
								if (are_equal(XOR(newBase, base.xorRows(mod)), f))
								{
									return k;
								}
							}
						}
					}
				}
			}
		}

		return k + 1;
	}

	public static ArrayList<Integer> find_nonzeros(int[] x)
	{
		ArrayList<Integer> nz = new ArrayList<Integer>();
		for (int i = 0; i < x.length; i++)
		{
			if (x[i] != 0) nz.add(i);
		}	
		return nz;
	}

	public static ArrayList<Integer> find_zeros(int[] x)
	{
		ArrayList<Integer> z = new ArrayList<Integer>();
		for (int i = 0; i < x.length; i++)
		{
			if (x[i] == 0) z.add(i);
		}
		return z;
	}

	public static int optimized_peralta_distance(Matrix base, int[] newBase, int[] f, int od) throws Exception
	{
		for (int i = 0; i < base.getDimension(); i++)
		{
			if (are_equal(base.getRow(i), f))
			{
				return 0;
			}
		}
		if (are_equal(newBase, f))
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
				if (are_equal(f, base.getRow(i)))
				{
					return true;
				}
			}
			return false;
		}

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

	public static int[] compute_distance(Matrix base, Matrix m) throws Exception
	{
		int n = m.getDimension();
		int[] distance = new int[n];
		for (int i = 0; i < n; i++)
		{
			int[] fi = m.getRow(i);
			distance[i] = peralta_distance(base, fi);
		}
		return distance;
	}

	public static int[] optimized_compute_distance(Matrix base, Matrix m, int[] newBase, int[] dist) throws Exception
	{
		int n = m.getDimension();
		int[] distance = new int[n];
		for (int i = 0; i < n; i++)
		{
			int[] fi = m.getRow(i);
			distance[i] = optimized_peralta_distance(base, newBase, fi, dist[i]);
		}
		return distance;
	}

	public static int[] compute_distance_graph_walk(Matrix base, Matrix m, int[] newBase, int[] dist) throws Exception
	{
		int n = m.getDimension();
		int[] distance = new int[n];
		for (int i = 0; i < n; i++)
		{
			int[] fi = m.getRow(i);
			boolean match = false;
			for (int j = 0; j < base.getDimension(); j++)
			{
				if (are_equal(base.getRow(j), fi))
				{
					distance[i] = 0;
					match = true;
				}
			}
			if (!match)
			{
				if (are_equal(newBase, fi))
				{
					distance[i] = 0;
				}
				else
				{
					int[] ft1 = new int[fi.length];
					for (int j = 0; j < ft1.length; j++) ft1[j] = 0;
					walk_ft = fi;
					walk_newBase = newBase;
					walk_od = dist[i] - 1;
					walk_base = base;
					distance[i] = walk_weighted_sequences_and_change_with_graph_v2(WALK_DIR.RIGHT, base.getDimension(), dist[i] - 1, 0, ft1);
				}
			}
		}
		return distance;
	}

	// the heuristic(T, A, B) works as follows:
	// 1. Iniitlaize array for T with 0 elements (unchosen)
	// 2. for each unchosen index i
	// 3. split A and B into two lists each (A0/A1 and B0/B1)
	// 4. set u = T[i]
	// 5. set pickedPositions(T,i) = 1 (set ith bit as chosen)
	// 6. if (u == 0) recursively solve (T, A0, B0) and (T, A1, B1), else solve (T, A0, B1) and (T, A1, B0)
	// ... this splits all the way down to single-element lists, which become the elements alpha and beta that are sought after

	// heuristic:
	//  Pick an unchosen position i from T.

	//  Split A into A0 and A1, where A0 are those elements
	//  of A that contain a 0 at position i and A1 are those
	//  elements of A that contain a 1 at position i.
	//  Do the same with B.

	//  u = T[i];

	//  declare position i "chosen";

	//  if (u == 0) recursively solve (A0,B0,T) and (A1,B1,T);
	//  else recursively solve (A0,B1,T) and (A1,B0,T);
	public static VectorPair split_subset_sum(int[] T, ArrayList<int[]> A, ArrayList<int[]> B) throws Exception
	{
		int[] tbits = new int[T.length]; // automatically initialized to zero when malloc'd

		// Check the base case - the size of A/B is one entry and the XOR is T
		// we can't split the list that has one element

		// System.out.println(A.toString());
		// System.out.println(B.toString());

		System.err.println(A.size() + ", " + B.size());

		if (A.size() == 1 || B.size() == 1)
		{
			for (int i = 0; i < A.size(); i++)
			{
				for (int j = 0; j < B.size(); j++)
				{
					if (are_equal(XOR(A.get(i), B.get(j)), T))
					{
						return new VectorPair(A.get(i), B.get(j)); // this pair added up to T
					}
				}
			}

			// there wasn't a pair that added up
			return null; // this pair didn't add up to T
		}

		// if either one of the lists is empty, then... whoops!
		if (A.size() == 0 || B.size() == 0) return null;

		for (int i = 0; i < tbits.length; i++)
		{
			if (tbits[i] == 0) // unchosen
			{
				ArrayList<int[]> A0 = new ArrayList<int[]>();
				ArrayList<int[]> A1 = new ArrayList<int[]>();
				ArrayList<int[]> B0 = new ArrayList<int[]>();
				ArrayList<int[]> B1 = new ArrayList<int[]>();

				// Split as per comments above
				for (int[] aentry : A)
				{
					if (aentry[i] == 0)
					{
						A0.add(aentry);
					}
					else
					{
						A1.add(aentry);	
					}
				}
				for (int[] bentry : B)
				{
					if (bentry[i] == 0)
					{
						B0.add(bentry);
					}
					else
					{
						B1.add(bentry);	
					}
				}

				// Set the ith bit as chosen
				int u = T[i];
				tbits[i] = 1; 

				// Recursively solve the two splits
				VectorPair vp = null;
				if (u == 0)
				{
					// (A0,B0,T) and (A1,B1,T);
					vp = split_subset_sum(T, A0, B0);
					if (vp == null)
					{
						vp = split_subset_sum(T, A1, B1);
					}

				}
				else
				{
					// (A0,B1,T) and (A1,B0,T);
					vp = split_subset_sum(T, A0, B1);
					if (vp == null)
					{
						vp = split_subset_sum(T, A1, B0);
					}
				}
				return vp;
			}
		}

		// Should not get here... 
		throw new Exception("Should NOT GET HERE!");
	}

	public static int[] compute_distance_target_split(Matrix base, Matrix m, int[] newBase, int[] dist) throws Exception
	{
		int n = m.getDimension();
		int[] distance = new int[n];
		for (int i = 0; i < n; i++)
		{
			// Input target T is each row of m (the set of targets)
			// k value is the corresponding index in the distance array
			int[] T = m.getRow(i);
			int k = dist[i];
			distance[i] = k;

			// #1 p1 and p2 are pointers to the starts of each split
			int sp1 = 0;                             // 0 - 3/2 = 1, 3/2 + 1 = 2 - 3
			int sp2 = (base.getDimension() / 2) + 1; // 0 - 4/2 = 2, 5/2 + 1 = 3 - 4

			// #2 pick a random a (and corresponding b...)
			// Random prng = new Random(System.currentTimeMillis());
			VectorPair pair = null;
			int a = 0;
			int b = k;
			while (pair == null && b >= 0)
			{
				error("Trying a/b = " + a + "/" + b);

				// #3 use hamiltonian path method to commpute sums of all a and b elements from each list and store in another list A/B - resulting in A = C(n/2,a) and B = C(n/2,b) sized lists
				ArrayList<int[]> Alist = hamiltonian_sums(base, newBase, 0, sp2 - 1, a); // |Alist| = C(n/2,a)
				ArrayList<int[]> Blist = hamiltonian_sums(base, newBase, sp2, base.getDimension() - 1, b); // |Blist| = C(n/2,b)

				// #4 invoke heuristic with Alist and Blist
				// TODO: see split_subset_sum method above
				pair = split_subset_sum(T, Alist, Blist);
				if (pair == null)
				{
					error("GAH! Didn't find a matching summation...");
				}
				else // found a matching sum, so the distance can decrease by at most one
				{
					dist[i] = dist[i] - 1;
					break;
				}

				a++;
				b--;
				// int[] T, ArrayList<int[]> A, ArrayList<int[]> B
			}


			// int a = prng.nextInt(k);
			// int b = k - a;

			// // #3 use hamiltonian path method to commpute sums of all a and b elements from each list and store in another list A/B - resulting in A = C(n/2,a) and B = C(n/2,b) sized lists
			// ArrayList<int[]> Alist = hamiltonian_sums(base, newBase, 0, sp2 - 1, a); // |Alist| = C(n/2,a)
			// ArrayList<int[]> Blist = hamiltonian_sums(base, newBase, sp2, base.getDimension() - 1, b); // |Blist| = C(n/2,b)

			// // #4 invoke heuristic with Alist and Blist
			// // TODO: see split_subset_sum method above
			// VectorPair pair = split_subset_sum(T, Alist, Blist);
			// if (pair == null)
			// {
			// 	System.err.println("GAH! Didn't find a matching summation...");
			// 	throw new Exception("halt and examine");
			// }
			// int[] T, ArrayList<int[]> A, ArrayList<int[]> B
		}

		return distance;

		// compute the distance as follows:
		// 1. split the rows of base into two lists 
		// 2. guess a and b (# elements from first half and second half, resp) such that a+b = k and a,b <= k
		// 3. With hamiltonian path method, commpute sums of all a and b elements from each list and store in another list A/B - resulting in A = C(n/2,a) and B = C(n/2,b) sized lists
		// 4. invoke heuristic with T, A and B

		// the heuristic(T, A, B) works as follows:
		// 1. Iniitlaize array for T with 0 elements (unchosen)
		// 2. for each unchosen index i
		// 3. split A and B into two lists each (A0/A1 and B0/B1)
		// 4. set u = T[i]
		// 5. set pickedPositions(T,i) = 1 (set ith bit as chosen)
		// 6. if (u == 0) recursively solve (T, A0, B0) and (T, A1, B1), else solve (T, A0, B1) and (T, A1, B0)
		// ... this splits all the way down to single-element lists, which become the elements alpha and beta that are sought after
	}

	public static BasePair pick_base(Matrix base, Matrix m, int[] dist, int tieBreaker, int distanceMethod) throws Exception
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
					int[] sum = add_mod(base.getRow(i), base.getRow(j), 2);
					if (!(base.containsRow(sum)) && !is_zero(sum)) // bases can't be zero and we don't repeat bases
					{
						int[] newDist = null;
						switch (distanceMethod)
						{
							case 1: 
								// newDist = compute_distance_graph_walk(base, m, sum, dist);
								newDist = compute_distance_graph_walk(base, m, sum, dist);
								break;
							case 2:
								newDist = compute_distance_target_split(base, m, sum, dist);
								break;
							default: // bank on recursive method
								// just FYI: sum is the new base
								newDist = optimized_compute_distance(base, m, sum, dist);
								break;
						}
						int newDistSum = sum(newDist);
						// disp("sum = " + newDistSum);
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
				maxNorm = euclidean_norm(optimalDistances.get(0));
				optimalBase = optimalBases.get(0);
				optimalDist = optimalDistances.get(0);
				mi = pairs.get(0).i;
				mj = pairs.get(0).j;
				for (int i = 1; i < optimalDistances.size(); i++)
				{
					double norm = euclidean_norm(optimalDistances.get(i));
					if (norm > maxNorm)
					{
						maxNorm = norm;
						optimalBase = optimalBases.get(i);
						optimalDist = optimalDistances.get(i);
						mi = pairs.get(i).i;
						mj = pairs.get(i).j;
					}
					else if (norm == maxNorm && compare_pairs(pairs.get(i).i, pairs.get(i).j, mi, mj) <= 0)
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
				maxNorm = Math.pow(euclidean_norm(optimalDistances.get(0)), 2);
				maxNorm = maxNorm - (double)largestDist;
				for (int i = 1; i < optimalDistances.size(); i++)
				{
					double norm = Math.pow(euclidean_norm(optimalDistances.get(i)), 2);
					norm = norm - (double)largest(optimalDistances.get(i));
					if (norm > maxNorm)
					{
						maxNorm = norm;
						optimalBase = optimalBases.get(i);
						optimalDist = optimalDistances.get(i);
						mi = pairs.get(i).i;
						mj = pairs.get(i).j;
					}
					else if (norm == maxNorm && compare_pairs(pairs.get(i).i, pairs.get(i).j, mi, mj) <= 0)
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

				maxNorm = Math.pow(euclidean_norm(optimalDistances.get(0)), 2);
				maxNorm = maxNorm - (double)(firstLargestDist - secondLargestDist);
				for (int i = 1; i < optimalDistances.size(); i++)
				{
					double norm = Math.pow(euclidean_norm(optimalDistances.get(i)), 2);
					norm = norm - (double)(firstLargestDist - secondLargestDist);
					if (norm > maxNorm)
					{
						maxNorm = norm;
						optimalBase = optimalBases.get(i);
						optimalDist = optimalDistances.get(i);
						mi = pairs.get(i).i;
						mj = pairs.get(i).j;
					}
					else if (norm == maxNorm && compare_pairs(pairs.get(i).i, pairs.get(i).j, mi, mj) <= 0)
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

	public static SLPStrings peralta_optimize(Matrix m, int r, int c, int tieBreaker, int distanceMethod) throws Exception
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
		int[] dist = compute_distance(base, m);
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
		while (is_zero(dist) == false)
		{
			BasePair newBase = pick_base(base, m, dist, tieBreaker, distanceMethod);
			base.appendRow(newBase.base);

			// Insert the new line
			String c1 = newBase.p.i < nVars ? "x" : "t";
			String c2 = newBase.p.j < nVars ? "x" : "t";
			slp.add("t_" + i + " = " + c1 + newBase.p.i + " XOR " + c2 + newBase.p.j);
			xorCount++;
			for (int f = 0; f < m.getDimension(); f++)
			{
				if (are_equal(newBase.base, m.getRow(f)))
				{
					slp.add("y_" + f + " = " + "t_" + i);
				}
			}
			i++;

			// Update the distance vector with this new base
			// dist = compute_distance(base, m);
			dist = newBase.dist;
		}

		return new SLPStrings(slp, xorCount, 0); // no ANDs required
	}

	public static ArrayList<MatrixState> exhaustive_optimize(Matrix currMatrix, int oldi, int oldj, int n, int lastcol) throws Exception
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
						ArrayList<MatrixState> rHistory = exhaustive_optimize(newMatrix, i, j, n, lastcol + 1);
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

	public static int compare_pairs(int i1, int j1, int i2, int j2)
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

	public static ArrayList<Vector> sort_rows_reverse_order(ArrayList<Vector> rows) throws Exception
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
					if (compare_reverse(ri, rj) == -1) // then swap
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

	// this builds an SLPStrings...
	public static ArrayList<String> bernstein_optimize(ArrayList<Vector> rows, int n, int p, int q) throws Exception
	{
		ArrayList<String> slp = new ArrayList<String>();
		rows = sort_rows_reverse_order(rows); // sort in reverse order - this is redundant at times

		if (q == 0) return slp;
		if (p == 0)
		{
			return slp;
		}
		if (rows.get(q - 1).row[p - 1] == 0) // x_{p-1} is unused, so don't consider next time around
		{
			slp.addAll(bernstein_optimize(rows, n, p - 1, q));
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
				slp.addAll(bernstein_optimize(rows, n, p, q));
				return slp;
			}
			else
			{
				rows.get(q - 1).row[p - 1] = 0; // set the row to 0...
				slp.add("y_" + rows.get(q - 1).id + " = x_" + (p - 1) + " XOR y_" + rows.get(q - 1).id);
				slp.addAll(bernstein_optimize(rows, n, p, q));
				return slp;
			}
		}

		return slp;
	}

	public static int num_gates(int n, ArrayList<MatrixState> history) throws Exception
	{
		if (history.size() == 0)
		{
			throw new Exception("Invalid history.");
		}
		return history.get(history.size() - 1).matrix.getGateCount();
	}

	public static ArrayList<String> build_slp(ArrayList<MatrixState> history, int vars) throws Exception
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
		// disp("STARTING PERALTA TEST");
		// int[][] test4 = {{1,1,1,0,0},{0,1,0,1,1},{1,0,1,1,1},{0,1,1,1,0},{1,1,0,1,0},{0,1,1,1,1}};
		// Matrix m4 = new Matrix(test4, 5);
		// disp(m4.toString());
		// SLPStrings slp4 = peralta_optimize(m4, 6, 5, 0);
		// SLPStrings slp4 = parallelperalta_optimize(m4, 6, 5, 0);
		// disp_strings(slp4.lines);
		// long start = System.currentTimeMillis();
		// SLPStrings slp5 = parallelperalta_optimize(m4, 6, 5, 0);
		// long end = System.currentTimeMillis();
		// disp("" + (end - start));
		// disp_strings(slp5.lines);
		// start = System.currentTimeMillis();
		// SLPStrings slp6 = parallelperalta_optimize_v2(m4, 6, 5, 0);
		// end = System.currentTimeMillis();
		// disp("" + (end - start));
		// disp_strings(slp6.lines);

		Vector v1 = new Vector(new int[] {1,1,0,0});
		Vector v2 = new Vector(new int[] {1,1,0,0});
		Vector v3 = new Vector(new int[] {0,0,0,0});
		disp("" + v1.equals(v2));
		disp("" + v1.equals(v3));

//		ArrayList<Vector> seqs = new ArrayList<Vector>();
//		int[] max = {1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		int n = 5;
		int k = 3;
		disp("Finding sequences for: " + n + "," + k);
		long start = System.currentTimeMillis();
		ArrayList<Vector> seqs = walk_weighted_sequences(n, k); 
		long end = System.currentTimeMillis();
		disp("Finding sequences for: " + n + "," + k + " with BigIntegers");
		// ArrayList<BigInteger> ints = walk_weighted_sequences_v2(n, k); 
		disp("Result #1: " + (end - start) + "ms");
		for (Vector s : seqs)
		{
			disp(s.row);
		}

		lastString = "";
		for (int i = 0; i < n; i++) lastString = lastString + "0";
		start = System.currentTimeMillis();
		ArrayList<Vector> seqs2 = walk_weighted_sequences_v2(WALK_DIR.RIGHT, n, k, "", "");
		end = System.currentTimeMillis();
		disp("Result #2: " + (end - start) + "ms");
		for (Vector s : seqs)
		{
			disp(s.row);
		}	

		disp("" + (seqs.size() == seqs2.size()));
	}

	public static ArrayList<Matrix> build_matrices(String file) throws Exception
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

			// error(m.toString());
			// error("STARTING PAAR TEST");
			// ArrayList<MatrixState> history = paar_optimize(m, m.getLength()); 
			// tgc = num_gates(m.getDimension(), history);
			// gc = tgc < gc ? tgc : gc;
			// // disp("Number of gates: " + num_gates(m.getDimension(), history));
			// if (debug) disp_strings(build_slp(history, m.getLength()));

			error("STARTING PERALTA TEST - TIE 0, RECURSIVE DISTANCE");
			long start4 = System.currentTimeMillis();
			SLPStrings slp4 = peralta_optimize(m, m.getDimension(), m.getLength(), 0, 0);
			long end4 = System.currentTimeMillis();
			tgc = slp4.xc;
			gc = tgc < gc ? tgc : gc;
			// error("Number of gates: " + slp4.xc);
			if (debug) disp_strings(slp4.lines);
			error("Elapsed time: " + (end4 - start4));

			error("STARTING PERALTA TEST - TIE 0, GRAPH DISTANCE");
			long start5 = System.currentTimeMillis();
			SLPStrings slp5 = peralta_optimize(m, m.getDimension(), m.getLength(), 0, 1);
			long end5 = System.currentTimeMillis();
			tgc = slp5.xc;
			gc = tgc < gc ? tgc : gc;
			// error("Number of gates: " + slp4.xc);
			error("Elapsed time: " + (end5 - start5));
			if (debug) disp_strings(slp5.lines);

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
			ArrayList<Matrix> matrices = build_matrices(args[0]);
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
				ArrayList<Matrix> matrices = build_matrices(fname);
				int ftotal = optimize(matrices, false);

				// Inverse
				fname = header + "_" + i + "_inverse.txt";
				matrices = build_matrices(fname);
				int itotal = optimize(matrices, false);

				// Merged
				fname = header + "_" + i + "_merged.txt";
				matrices = build_matrices(fname);
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

