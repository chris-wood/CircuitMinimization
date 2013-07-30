import java.util.Random;
import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Scanner;
import java.util.Random;

public class NonlinearOptimize
{
	static void error(String s) { System.err.println(s); }
	static void error(int[] v) {
		System.err.print("[");
		for (int i = 0; i < v.length - 1; i++) System.err.print(v[i] + " ");
		System.err.println(v[v.length - 1] + "]");
	}
	static void disp(String s) { System.out.println(s); }
	static void disp(ArrayList<String> strings) { for (String s : strings) System.out.println(s); }
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

	public enum RoundReturn { SIGNAL_NOT_FOUND, SIGNAL_MATCH, AND_EXCEEDED, SIGNALS_EXCEEDED }

	public static int evaluate(Matrix f, int[] x) throws Exception
	{
		int sum = 0;

		// Evaluate each term and then XOR the results together...
		for (int t = 0; t < f.getDimension(); t++)
		{
			int[] ft = f.getRow(t);
			if (x.length != ft.length) 
				throw new Exception("Evaluate input doesn't match function domain");
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

	public static int[] buildSignal(Matrix f) throws Exception
	{
		ArrayList<Signal> knowns = initializeKnown(f.getLength());
		ArrayList<int[]> inputs = new ArrayList<int[]>();
		for (int i = 0; i < knowns.size(); i++)
		{
			inputs.add(knowns.get(i).signal);	
		}
		Matrix known = new Matrix(inputs);
		int[] signal = new int[inputs.get(0).length];
		for (int kf = 0; kf < known.getLength(); kf++) 
		{
			int[] x = known.getColumn(kf);
			signal[kf] = evaluate(f, x);
		}
		return signal;
	}

	public static ArrayList<Signal> initializeKnown(int nvars) throws Exception
	{
		ArrayList<Signal> inputs = new ArrayList<Signal>();
		int run = 1;
		int order = 1 << nvars;
		for (int i = 0; i < nvars; i++)
		{
			int[] signal = new int[order];
			int j = 0;
			while (j < order)
			{
				for (int r = 0; r < run; r++) // skip the run
				{
					j++;
				}
				for (int r = 0; r < run && j < order; r++)
				{
					signal[j] = 1;
					j++;
				}
			}
			run *= 2;
			// disp(signal);
			Signal s = new Signal(signal, i);
			inputs.add(s);
		}
		return inputs;
	}

	// and/xor counts here
	static int xorCount = 0;
	static int andCount = 0;

	public static RoundReturn xorRound(ArrayList<Signal> used, 
		ArrayList<Signal> known, int[] signal, int numSignals, 
		double xp, int andCap) throws Exception
	{
		Random r = new Random(System.currentTimeMillis());
		for (int i = 0; i < numSignals; i++)
		{
			int u = r.nextInt(known.size() + used.size());
			int v = r.nextInt(known.size() + used.size());
			while (u == v) // enforce distinct signals
			{
				v = r.nextInt(known.size() + used.size());
			}
			if (r.nextDouble() < xp)
			{
				Signal su, sv;

				if (u < known.size()) su = known.get(u);
				else su = used.get(u - known.size());
				if (v < known.size()) sv = known.get(v);
				else sv = used.get(v - known.size());

				Signal x = Signal.xor(su, sv);
				if (equal(x.signal, signal)) 
				{
					known.add(x);
					return RoundReturn.SIGNAL_MATCH;
				}
				if (x.ac >= andCap) return RoundReturn.AND_EXCEEDED;
				known.add(x); // add the new signal
				if (known.size() == numSignals - 1) return RoundReturn.SIGNALS_EXCEEDED;
			}
		}
		return RoundReturn.SIGNAL_NOT_FOUND;
	}

	public static RoundReturn andRound(ArrayList<Signal> used, 
		ArrayList<Signal> known, int[] signal, int numSignals, 
		double ap, int andCap) throws Exception
	{
		Random r = new Random(System.currentTimeMillis());
		for (int i = 0; i < numSignals; i++)
		{
			int u = r.nextInt(known.size() + used.size());
			int v = r.nextInt(known.size() + used.size());
			while (u == v) // enforce distinct signals
			{
				v = r.nextInt(known.size() + used.size());
			}
			if (r.nextDouble() < ap)
			{
				Signal su, sv;

				if (u < known.size()) su = known.get(u);
				else su = used.get(u - known.size());
				if (v < known.size()) sv = known.get(v);
				else sv = used.get(v - known.size());

				Signal a = Signal.and(su, sv);

				if (equal(a.signal, signal)) 
				{
					known.add(a);
					return RoundReturn.SIGNAL_MATCH;
				}
				if (a.ac >= andCap) return RoundReturn.AND_EXCEEDED;
				known.add(a); // add the new signal
				if (known.size() == numSignals - 1) return RoundReturn.SIGNALS_EXCEEDED;
			}
		}
		return RoundReturn.SIGNAL_NOT_FOUND;
	}

	public static Signal optimizeEquation(Matrix f, int[] signal, 
		ArrayList<Signal> used, int numXorSignals, int numAndSignals, 
		int andCap, double xp, double ap, long rounds, long tries) 
		throws Exception
	{
		if (isZero(signal)) // a zero signal gets a zero input.
		{
			return null;
		}

		ArrayList<Signal> known = null;
		for (long t = 0; t < tries; t++)
		{
			long rc = 0;
			known = initializeKnown(f.getLength()); 

			xorCount = 0;
			andCount = 0;

			RoundReturn terminate = RoundReturn.SIGNAL_NOT_FOUND;
			while (rc < rounds && terminate == RoundReturn.SIGNAL_NOT_FOUND)
			{
				terminate = xorRound(used, known, signal, numXorSignals, xp, andCap);
				if (terminate == RoundReturn.SIGNAL_MATCH || terminate == RoundReturn.AND_EXCEEDED 
					|| terminate == RoundReturn.SIGNALS_EXCEEDED) break;
				terminate = andRound(used, known, signal, numAndSignals, ap, andCap);
				if (terminate == RoundReturn.SIGNAL_MATCH || terminate == RoundReturn.AND_EXCEEDED 
					|| terminate == RoundReturn.SIGNALS_EXCEEDED) break;
				rc++; // next set of rounds...
			}
			if (terminate == RoundReturn.SIGNAL_MATCH) 
			{
				// error("Found a matching circuit.");
				error(known.get(known.size() - 1).signal);
				error(known.get(known.size() - 1).rep);
				return known.get(known.size() - 1);
			}
		}

		return null; // the last function
	}

	public static ArrayList<ArrayList<Integer>> permute(ArrayList<Integer> set) throws Exception
	{
		if (set.size() == 1)
		{	
			ArrayList<ArrayList<Integer>> sets = new ArrayList<ArrayList<Integer>>();
			sets.add(set);
			return sets;
		}
		else
		{
			ArrayList<ArrayList<Integer>> sets = new ArrayList<ArrayList<Integer>>();
			for (int i = 0; i < set.size(); i++)
			{
				ArrayList<Integer> sc = new ArrayList<Integer>();
				for (int j = 0; j < set.size(); j++) sc.add(set.get(j));
				ArrayList<Integer> s = new ArrayList<Integer>();
				s.add(set.get(i));
				sc.remove(set.get(i));
				ArrayList<ArrayList<Integer>> smallerPermutations = permute(sc);
				for (ArrayList<Integer> sp : smallerPermutations)
				{
					ArrayList<Integer> ssc = new ArrayList<Integer>(s);
					ssc.addAll(sp);
					sets.add(ssc);
				}
			}
			return sets;
		}
	}

	public static ArrayList<Signal> peraltaSystemOptimization(PeraltaCircuit c) throws Exception
	{
		ArrayList<Matrix> forms = c.forms;
		int numXorSignals = c.numXorSignals;
		int numAndSignals = c.numAndSignals;
		int andCap = c.andCap;
		long rounds = c.rounds; 
		long tries = c.tries;
		double xp = c.xorProb;
		double ap = c.andProb;
		ArrayList<String> opt = new ArrayList<String>();

		Random r = new Random(System.currentTimeMillis());

		// Build permutations of the linear forms
		ArrayList<Integer> indices = new ArrayList<Integer>();
		for (int i = 0; i < forms.size(); i++) 
		{
			indices.add(i);
		}
		ArrayList<ArrayList<Integer>> permutations = permute(indices);
		while (true) // continue until we find the circuit.
		{
			xp = ap = 0.0;
			while (xp < 0.5)
			{
				xp = r.nextDouble();
			}
			while (ap < 0.5)
			{
				ap = r.nextDouble();
			}

			error("New run");
			error("xp = " + xp);
			error("ap = " + ap);
			error("original XOR count = " + c.origXorSignals);
			error("original AND count = " + c.origAndSignals);

			ArrayList<Signal> solution = new ArrayList<Signal>();
			ArrayList<Signal> used = new ArrayList<Signal>();
			for (int i = 0; i < permutations.get(0).size(); i++)
			{
				solution.add(null);
			}
			for (ArrayList<Integer> index : permutations)
			{
				// ArrayList<Signal> solution = new ArrayList<Signal>();
				boolean pass = true;
				error("Trying permutation: " + index.toString());
				for (int i = 0; i < index.size(); i++)
				{
					int ii = index.get(i);
					if (solution.get(ii) == null)
					{
						Signal s = optimizeEquation(forms.get(ii), buildSignal(forms.get(ii)), 
							used, numXorSignals, numAndSignals, andCap, xp, ap, rounds, tries);
						if (s == null)
						{
							pass = false;
							error("Couldn't find linear circuit for form: " + ii);
							// break;
						}
						else
						{
							error("Signal found for form: " + ii);
							// usedSignals.addAll(s.decompose());

							// permanent storage
							solution.set(ii, s);
							ArrayList<Signal> decomposition = s.decompose();
							for (Signal sig : decomposition)
							{
								if (!(used.contains(sig)))
								{
									used.add(sig); // only add circuits we haven't build before
								}
							}
						}
					}
				}

				// Check to see if all have been found...
				boolean sat = true;
				int xorCount = 0;
				int andCount = 0;
				int totalXorCount = 0; // for the input circuits
				int totalAndCount = 0; // for the input circuits
				for (int i = 0; i < index.size(); i++)
				{
					if (solution.get(i) == null)
					{
						sat = false;
						xorCount += c.origXorSignals.get(i);
						andCount += c.origAndSignals.get(i);
					}
					else
					{
						xorCount += solution.get(i).xc;
						andCount += solution.get(i).ac;
					}

					// Tally up the total gate count from the original input circuit
					totalXorCount += c.origXorSignals.get(i);
					totalAndCount += c.origAndSignals.get(i);
				}

				error("Original AND: " + totalAndCount);
				error("Original XOR: " + totalXorCount);
				error("New AND:      " + andCount);
				error("New XOR:      " + xorCount);

				if (sat && (andCount < totalAndCount)) // we found an optimized circuit for ALL input equations
				{
					return solution;
				} 
				else if (andCount < totalAndCount) // we still found a little better solution...
				{
					return solution;
				}
			}
		}
	}

	public static ArrayList<PeraltaCircuit> buildMatrices(String file) throws Exception
	{
		ArrayList<PeraltaCircuit> circuits = new ArrayList<PeraltaCircuit>();
		BufferedReader in = new BufferedReader(new FileReader(file)); 
		String line = "";
		while ((line = in.readLine()) != null && !line.isEmpty())
		{
			BufferedReader fin = new BufferedReader(new FileReader(line));
			String matRow = "";
			int p = 0;
			int q = 0;
			ArrayList<Matrix> forms = new ArrayList<Matrix>();

			// Read in the parameters
			String[] data = fin.readLine().split(" ");
			int xc = Integer.parseInt(data[0]);
			int ac = Integer.parseInt(data[1]);
			int cap = Integer.parseInt(data[2]);
			double xp = Double.parseDouble(data[3]);
			double ap = Double.parseDouble(data[4]);
			long rounds = Long.parseLong(data[5]);
			long tries = Long.parseLong(data[6]);
			int nvars = Integer.parseInt(data[7]);
			int nforms = Integer.parseInt(data[8]);

			ArrayList<Integer> signalAnds = new ArrayList<Integer>();
			ArrayList<Integer> signalXors = new ArrayList<Integer>();

			// Read the matrix forms...
			for (int i = 0; i < nforms; i++)
			{
				int ands = 0;
				int xors = 0;
				String form = fin.readLine();
				disp(form);
				ArrayList<int[]> rows = new ArrayList<int[]>();
				Scanner s = new Scanner(form);
				int[] row = null;
				int pc = 0;

				while (s.hasNext())
				{
					row = new int[nvars];
					pc = 0;
					while (s.hasNext())
					{
						String tok = s.next();
						if (tok.equals("+")) // or
						{
							xors++;
							rows.add(row);
							break;
						}
						else if (!(tok.equals("^"))) // and
						{
							int index = Integer.parseInt(tok); 
							if (index < 0 || index >= nvars) 
								throw new Exception("Invalid variable index");
							row[index] = 1; 
							pc++;
						}

						// special case
						if (tok.equals("^")) ands++;
					}
				}

				if (pc > 0)
				{
					rows.add(row);
				}

				Matrix m = new Matrix(rows);
				disp(m.toString());
				forms.add(m);

				// Add the number of signals for this particular form
				signalAnds.add(ands);
				signalXors.add(xors);
			}

			PeraltaCircuit circuit = new PeraltaCircuit(forms, xc, ac, xp, ap, cap, rounds, tries);
			circuit.origAndSignals = signalAnds;
			circuit.origXorSignals = signalXors;
			circuits.add(circuit);
		} 
		return circuits;
	}

	public static String matrixToFormString(Matrix m) throws Exception
	{
		String form = "";

		for (int r = 0; r < m.getDimension(); r++)
		{
			String tmp = "(";
			for (int c = 0; c < m.getLength(); c++)
			{
				if (m.getMatrix()[r][c] == 1)
				{
					tmp = tmp + "(" + c + ")";
				}
			}
			tmp = tmp + ")";
			if (r == (m.getDimension() - 1))
			{
				form = form + tmp;
			}
			else
			{
				form = form + tmp + "+";
			}
		}

		return form;
	}

	public static String matrixToSLP(Matrix m) throws Exception
	{
		String slp = "";

		for (int r = 0; r < m.getDimension(); r++)
		{
			String tmp = "(";
			int[] row = m.getRow(r);
			if (!(isZero(row)))
			{
				for (int c = 0; c < m.getLength(); c++)
				{
					if (row[c] == 1)
					{
						tmp = tmp + Signal.GLOBAL_VAR_STRING + c + " AND ";
					}
				}
				tmp = tmp.substring(0, tmp.length() - 4) + ")"; // strip last AND
				if (r == (m.getDimension() - 1))
				{
					slp = slp + tmp;
				}
				else
				{
					slp = slp + tmp + " + XOR + ";
				}
			}
		} 

		return slp;
	}

	public static void main(String[] args) throws Exception
	{
		if (args.length != 1)
		{
			System.err.println("usage: java NonlinearOptimize fileListFile");
			System.exit(-1);
		}

		// Read in parameters...
		ArrayList<PeraltaCircuit> circuits = buildMatrices(args[0]);
		disp("Rebuilt the following peralta circuits:");
		for (PeraltaCircuit c : circuits)
		{
			ArrayList<Signal> forms = peraltaSystemOptimization(c);
			if (forms != null)
			{
				System.out.println("Found a satisfying collection.");
				int xc = 0;
				int ac = 0;
				for (int i = 0; i < forms.size(); i++)
				{
					Signal s = forms.get(i);
					if (s != null)
					{
						xc += s.xc;
						ac += s.ac;
						disp(s.rep);
					}
					else
					{
						xc += c.origXorSignals.get(i);
						ac += c.origAndSignals.get(i);
						disp(matrixToFormString(c.forms.get(i)));
					}
				}
				System.out.println("XOR/AND count: " + xc + " " + ac);
			}
		}
	}
}
