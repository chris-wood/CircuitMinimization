import java.util.HashSet;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

public class ParallelMatrixOptimize extends RecursiveAction 
{
	public static void disp(String s) { System.out.println(s); }
	public static void disp(int[] v) {
		System.out.print("[");
		for (int i = 0; i < v.length - 1; i++) System.out.print(v[i] + " ");
		System.out.println(v[v.length - 1] + "]");
	}

	public ParallelSolution solution;
	public Matrix currMatrix;
	public int n;
	public int gOptimal;
	public int lastcol;
	public int oldi;
	public int oldj;

	public ParallelMatrixOptimize(ParallelSolution solution) 
	{
		this.solution = solution;
		this.currMatrix = solution.currMatrix;
		this.oldi = solution.oldi;
		this.oldj = solution.oldj;
		this.n = solution.n;
		this.lastcol = solution.lastcol;
		this.gOptimal = solution.gOptimal;
		// this.seen = seen; // for pruning...
	}

	public int weight(int[] v) 
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

	public boolean isZero(int[] v)
	{
		for (int i = 0; i < v.length; i++)
		{
			if (v[i] != 0) return false;
		}
		return true;
	}

	public int[] add(int[] a, int[] b) throws Exception
	{
		int[] c = new int[a.length];
		if (a.length != b.length) throw new Exception("Vector dimensions don't match for addition.");
		for (int i = 0; i < a.length; i++)
		{
			c[i] = a[i] + b[i];
		}
		return c;
	}

	public int[] AND(int[] x, int[] y) throws Exception
	{
		if (x.length != y.length) throw new Exception("Invalid x/y dimensions.");
		int[] and = new int[x.length];
		for (int i = 0; i < x.length; i++) 
		{
			and[i] = x[i] & y[i];
		}
		return and;
	}

	public int[] NEGATE(int[] x) 
	{
		int[] nx = new int[x.length];
		for (int i = 0; i < x.length; i++) 
		{
			nx[i] = x[i] == 0 ? 1 : 0; // enforce only two possibilities
		}
		return nx;
	}

	protected void compute()
	{
		try
		{
			ArrayList<Pair> pairs = new ArrayList<Pair>();
			for (int i = 0; i <= lastcol - 1; i++)
			{
				for (int j = i + 1; j <= lastcol; j++)
				{
					if (!(i < oldi && j != oldi && j != oldj && j < lastcol)) // prune
					{
						int[] coli = currMatrix.getColumn(i);
						int[] colj = currMatrix.getColumn(j);
						int wt = weight(AND(coli, colj));
						if (wt > 1) 
						{
							pairs.add(new Pair(i, j));
						}
					}
				} 
				
			}
			if (pairs.size() > 0)
			{
				exhaustiveOptimize(pairs); // breadth first...
			}
		}
		catch (Exception e)
		{
			System.err.println("Error in compute method");
			e.printStackTrace();
			System.exit(-1);
		}
	}

	public void exhaustiveOptimize(ArrayList<Pair> pairs) throws Exception
	{
		int wt = 0;
		int hmax = 0;
		HashSet<ParallelMatrixOptimize> tasks = new HashSet<ParallelMatrixOptimize>();

		for (Pair p : pairs)
		{
			int i = p.i;
			int j = p.j;
			int[] coli = currMatrix.getColumn(i);
			int[] colj = currMatrix.getColumn(j);
			int[] newcol = AND(coli, colj);
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

			ParallelSolution sol = new ParallelSolution(newMatrix.getGateCount(), 
				newMatrix, i, j, n, lastcol + 1);
			tasks.add(new ParallelMatrixOptimize(sol));
		}

		invokeAll(tasks); // fork & join
		for (ParallelMatrixOptimize result : tasks)
		{
			if (result.solution.gOptimal < this.solution.gOptimal)
			{
				this.solution = result.solution;
			}
		}
	}

	public static void main(String[] args) throws Exception
	{
		// int[][] test3 = {{1,1,0,1},{0,0,0,0},{0,0,1,0},{1,1,0,1}};
		int[][] test2 = {{1,1,0,1,1,1,0,0},{0,1,1,1,1,1,1,0},{0,1,0,1,1,1,0,0},{1,1,0,1,0,0,1,0},{1,0,0,1,1,0,1,0},{1,1,0,1,1,0,1,0},{1,1,0,1,0,0,0,0},{1,1,1,1,1,1,1,1}};
		Matrix m3 = new Matrix(test2, 8);
		disp(m3.toString());
		ParallelSolution solution = new ParallelSolution(m3.getGateCount(), m3, 0, 0, 8, 7);
		// ParallelMatrixOptimize optimizer = new ParallelMatrixOptimize(m3, 0, 0, 8, 7, m3.getGateCount(), new ArrayList<Pair>());
		ParallelMatrixOptimize optimizer = new ParallelMatrixOptimize(solution);
		ForkJoinPool pool = new ForkJoinPool();
		long start = System.currentTimeMillis();
		pool.invoke(optimizer);
		long end = System.currentTimeMillis();
		disp(optimizer.solution.currMatrix.toString());
		disp("" + optimizer.solution.gOptimal);
		disp("" + (end - start) + " ms");
	}
}

