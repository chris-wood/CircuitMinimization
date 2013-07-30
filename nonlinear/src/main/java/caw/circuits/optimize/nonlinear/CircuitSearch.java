import java.util.Scanner;
import java.util.ArrayList;;

public class CircuitSearch
{
	static void error(String s) { System.err.println(s); }
	static void disp(String s) { System.out.println(s); }
	static void disp(int n) { System.out.println(n); }
	static void disp(int[] v) {
		System.out.print("[");
		for (int i = 0; i < v.length - 1; i++) System.out.print(v[i] + " ");
		System.out.println(v[v.length - 1] + "]");
	}
	static int[] copy(int[] v)
	{
		int[] c = new int[v.length];
		for (int i = 0; i < v.length; i++) c[i] = v[i];
		return c;
	}
	static boolean equal(int[] a, int[] b)
	{
		for (int i = 0; i < a.length; i++)
		{
			if (a[i] != b[i]) return false;
		}
		return true;
	}
	static boolean isZero(int[] v) 
	{
		for (int i = 0; i < v.length; i++)
		{
			if (v[i] != 0) return false;
		}
		return true;
	}
	public boolean isZero(int[][] m)
	{
		for (int i = 0; i < m.length; i++)
		{
			for (int j = 0; j < m[i].length; j++)
			{
				if (m[i][j] != 0) return false;
			}
		}
		return true;
	}

	public static int evaluate(Matrix f, int[] x) throws Exception
	{
		int sum = 0;

		// Evaluate each term and then XOR the results together...
		for (int t = 0; t < f.getDimension(); t++)
		{
			int[] ft = f.getRow(t);
			if (x.length != ft.length) throw new Exception("Evaluate input doesn't match function domain");
			int prod = 0;
			if (!isZero(ft)) 
			{
				prod = 1;
				for (int i = 0; i < x.length; i++)
				{
					if (x[i] == 0 && ft[i] == 1)  // #ft == #x (the number of variables)
					{
						prod = 0; // an x term in the prodct was 0, so the product immediately evaluates to 0
						break;
					}
				}
			}
			sum += prod;
		}

		// BF output is 1/0
		return sum % 2; 
	}

	public static boolean computesSignal(Matrix f, int[] signal, int n) throws Exception
	{
		int[] x = new int[n];
		for (int i = 0; i < signal.length; i++)
		{
			if (evaluate(f, x) == signal[i])
			{
				for (int j = 0; j < n; j++)
				{
					if (x[j] == 0) 
					{
						x[j] = 1;
						break;
					} 
					x[j] = 0;
				}
			}
			else
			{
				return false;
			}
		}
		return true;
	}

	// Compute the "cost" of a nonlinear predicate that contains addition and multiplication operations
	public static int cost(Matrix f, int n, int addCost, int multCost) throws Exception
	{
		int c = 0;
		int nonzeroRows = 0;
		for (int i = 0; i < f.getDimension(); i++)
		{
			int sum = 0;
			int[] r = f.getRow(i);
			if (!isZero(r)) 
			{
				nonzeroRows++;
				for (int j = 0; j < n; j++)
				{
					sum += r[j];
				}
				if (sum > 1) c += (sum * (multCost - 1));
			}
		}
		if (nonzeroRows > 1) c += nonzeroRows * (addCost - 1);
		
		return c;
	}

	public static String buildLinearForm(int id, Matrix f, int mode) throws Exception
	{
		String form = "";
		switch (mode)
		{
			case 0: form = "y" + id + " = "; break;
			case 1: form = "y_" + id + " = & "; break;
			case 2: break; // string is default
		}
		if (f.isZero()) 
		{
			form = form + "0";
			return form;
		}
		else
		{
			for (int i = 0; i < f.getDimension(); i++) // dimension == row
			{
				int[] row = f.getRow(i);
				if (!isZero(row))
				{
					String prod = "";
					for (int j = 0; j < row.length; j++)
					{
						if (row[j] == 1)
						{
							switch (mode)
							{
								case 0: prod = prod + "(x" + j + ")"; break;
								case 1: prod = prod + "x_" + j + ""; break;
								case 2: prod = prod + j + " ^ ";
							}
						}
					}
					switch (mode)
					{
						case 0: form = form + prod + " + "; break;
						case 1: form = form + prod + " + "; break;
						case 2: form = form + prod.substring(0, prod.lastIndexOf("^") - 1) + " + "; break;
					}
				}
			}
			switch (mode)
			{
				case 0: form = form.substring(0, form.lastIndexOf('+')); break;
				case 1: form = form.substring(0, form.lastIndexOf("+")) + "\\\\"; break;
				case 2: form = form.substring(0, form.lastIndexOf("+") - 1); break;
			}
			return form;
		}
	}

	// TODO
	public static Matrix findMinimalEquation(int[] signal, int n, int addCost, int multCost) throws Exception
	{
		Matrix eqn = null;
		int maxCost = (multCost * (n - 1)) + (addCost * (n - 1));

		// Zero function is special
		if (isZero(signal)) return null;

		// Initialize the products to choose from
		int np = 1 << n; // for 2 vars, this yields 4
		disp("Number of predicates: " + np);
		ArrayList<int[]> rows = new ArrayList<int[]>();
		for (int i = 0; i < np; i++)
		{
			int[] prod = new int[n];
			for (int j = 0; j < n; j++)
			{
				if (((1 << j) & i) > 0) prod[j] = 1;
			}
			rows.add(prod);
		}
		Matrix eqns = new Matrix(rows);

		// Now try all possibilities
		boolean foundOne = false;
		int nf = 1 << np;
		for (int i = 0; i < nf; i++)
		{
			ArrayList<int[]> preds = new ArrayList<int[]>();
			for (int f = 0; f < np; f++)
			{
				if (((1 << f) & i) > 0) preds.add(eqns.getRow(f));
			}
			if (preds.size() != 0)
			{
				Matrix tmpEqn = new Matrix(preds);
				if (computesSignal(tmpEqn, signal, n) == true)
				{
					error("Found a matching circuit...");
					int cost = cost(tmpEqn, n, addCost, multCost);
					if (cost < maxCost || !foundOne)
					{
						eqn = tmpEqn;
						maxCost = cost;
						foundOne = true;
					}
				}
			}
		}

		// This is an exhaustive search over every Boolean function - there will be at least one that works.
		// SANITY CHECK
		if (eqn == null) throw new Exception("Error: didn't find a matching predicate");

		return eqn;
	}

	// Read in the TT of an (n,m) circuit in canonical form (i.e. a hex string)
	// Note: no error handling - use wisely.
	public static void main(String[] args) throws Exception
	{
		int[] test = {0,0,0,1};
		Matrix tm = findMinimalEquation(test, 2, 1, 1);
		error("Sanity check: y = x_1x_2 (0 0 0 1)");
		error(tm.toString());
		error(buildLinearForm(0, tm, 0));

		// Check to see if LaTeX output is turned on...
		int mode = 0;
		if (args.length == 1)
		{
			mode = Integer.parseInt(args[0]);
			switch (mode)
			{
				case 0: disp("Normal output mode."); break;
				case 1: disp("LaTeX output mode."); break;
				case 2: disp("NonlinearOptimize.java output mode."); break;
			}
		}

		// Read in each hex string, one at a time, and find a minimal set of circuits for each one
		Scanner s = new Scanner(System.in);
		disp("Enter your GF(2^4) mapping as a hex string (GF(2^8) support to be added soon)");
		int n = 4;
		int n2n = (int)Math.pow(2, n);
		while (s.hasNext())
		{
			String tts = s.next();
			Matrix tt = new Matrix(n); // TODO: extend to high prime powers in the future
			for (int i = 0; i < tts.length(); i++)
			{
				int[] row = new int[n]; // four bits in a single entry

				int mrow = Integer.parseInt(tts.substring(i, i+1), n2n);
				for (int m = 0; m < n; m++)
				{
					if (((1 << (n - m - 1)) & mrow) > 0) row[m] = 1;
				}
				if (i == 0) 
				{
					tt.putRow(row, 0);
				}
				else
				{
					tt.appendRow(row);
				}
			}

			// Now, for each column in the truth table matrix (the coordinate functions), find a circuit!
			ArrayList<String> forms = new ArrayList<String>();
			for (int m = 0; m < tt.getLength(); m++)
			{
				int[] signal = tt.getColumn(m);
				disp("Searching for minimal circuit for the given signal: " + signal.length);
				disp(signal);
				Matrix eqn = findMinimalEquation(signal, 4, 1, 1);
				String f = buildLinearForm(m, eqn, mode);
				disp(f);
				forms.add(f);
				disp(eqn.toString());
			}
			disp("Final system of equations:");
			for (int m = 0; m < tt.getLength(); m++)
			{
				disp(tt.getColumn(m));
			}
			for (String ss : forms)
			{
				disp(ss);
			}
		}
	}
}
