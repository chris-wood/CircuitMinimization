import java.util.Random;
import java.util.ArrayList;
import java.util.HashMap;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Scanner;

public class SignalGraph
{
	public ArrayList<Signal> vertices;
	public ArrayList<Signal> initial;
	public ArrayList<Integer> leaves;
	public ArrayList<Integer> locked;
	public ArrayList<Pair> edges; // DAG, this is an ordered pair.
	public ArrayList<int[]> targets;
	public int numMatches;
	public int[][] distance;

	public int solAndCount;
	public int solXorCount;

	public enum RoundReturn { SIGNAL_NOT_FOUND, SIGNAL_MATCH, AND_EXCEEDED, XOR_EXCEEDED, SIGNALS_EXCEEDED }

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

	static ArrayList<Signal> initializeKnown(int nvars) throws Exception
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

	public SignalGraph(ArrayList<Signal> initialSignals, ArrayList<int[]> targets) throws Exception
	{
		numMatches = solAndCount = solXorCount = 0;
		initial = initialSignals;
		leaves = new ArrayList<Integer>();
		vertices = new ArrayList<Signal>();
		locked = new ArrayList<Integer>();
		int id = 0;
		for (Signal s : initialSignals) addSignal(s, id++);
		edges = new ArrayList<Pair>();

		this.targets = new ArrayList<int[]>();
		for (int t = 0; t < targets.size(); t++)
		{
			int len = targets.get(t).length;
			this.targets.add(new int[len]);
			for (int i = 0; i < len; i++) this.targets.get(t)[i] = targets.get(t)[i];
		}
	}

	public int fanout(int vi)
	{
		int fo = 0;
		for (Pair p : edges)
		{
			if (p.i == vi) fo++;
		}
		return fo;
	}

	public void addSignal(Signal x, int id) throws Exception
	{
		x.id = id;
		vertices.add(x);
		leaves.add(id);
		int newIndex = vertices.size() - 1;
		// disp("Adding initial signal: " + x.id);
	}

	// Add a signal and update the collection of leaves as well...
	public void addAndConnectSignal(Signal x, int op1_index, int op2_index) throws Exception
	{
		vertices.add(x);
		int newIndex = vertices.size() - 1;
		x.id = newIndex;
		leaves.add(newIndex);
		edges.add(new Pair(op1_index, newIndex));
		leaves.set(op1_index, 0);
		edges.add(new Pair(op2_index, newIndex));
		leaves.set(op2_index, 0);
	}

	public void updateDistances() 
	{
		int vc = vertices.size();
		distance = new int[vc][vc];
		int[][] edgeTo = new int[vc][vc];
		for (int i = 0; i < vc; i++)
		{
			for (int j = 0; j < vc; j++)
			{
				distance[i][j] = vertices.size();
			}
		}
		for (Pair e : edges) 
		{
			edgeTo[e.i][e.j] = 1;
			distance[e.i][e.j] = 1; // all edges have unit weight
		}

		// Floyd Warshall algorithm for min. distance
		for (int i = 0; i < vc; i++)
		{
			for (int v = 0; v < vc; v++)
			{
				if (edgeTo[v][i] == 0) continue;
				for (int w = 0; w < vc; w++)
				{
					if (distance[v][w] > distance[v][i] + distance[i][w])
					{
						distance[v][w] = distance[v][i] + distance[i][w];
					}
				}
			}
		}
	}

	public ArrayList<Integer> signalPath(Signal s)
	{
		ArrayList<Integer> path = new ArrayList<Integer>();
		ArrayList<Signal> queue = new ArrayList<Signal>();
		queue.add(s.siga);
		queue.add(s.sigb);

		while (queue.size() > 0)
		{
			Signal t = queue.get(queue.size() - 1);
			queue.remove(queue.size() - 1);
			if (t != null)
			{
				path.add(t.id);
				queue.add(t.siga);
				queue.add(t.sigb);
			}
		}

		return path;
	}

	public RoundReturn xorRound(int numSignals, int xorCap) throws Exception
	{
		Random r = new Random(System.currentTimeMillis());
		for (int i = 0; i < numSignals; i++)
		{
			Signal x = null;
			int u = 0;
			int v = 0;
			boolean repeat = false;
			int c = 0;
			do
			{
				if (repeat)
				{	
					c++;
					// disp("Trying to find a different XOR combination: " + c);
				}
				u = r.nextInt(vertices.size());
				v = r.nextInt(vertices.size());
				while (u == v) // enforce distinct signals
				{
					v = r.nextInt(vertices.size());
				}
				// disp("" + vertices.get(u).id);
				// disp("" + vertices.get(v).id);
				x = Signal.xor(vertices.get(u), vertices.get(v));
				// disp("vertices.contains = " + vertices.contains(x));
				repeat = true;
			} 
			while (vertices.contains(x));

			// disp("Selecting signals " + u + " and " + v + " for XOR round");

			boolean match = false;
			for (int ti = 0; ti < targets.size(); ti++)
			{
				int[] target = targets.get(ti);
				if (equal(x.signal, target)) 
				{
					disp("FOUND ONE SIGNAL.");
					disp(target);
					// System.exit(-1);
					addAndConnectSignal(x, u, v);

					// lock the signals along this path... 
					ArrayList<Integer> sigPath = signalPath(x);
					for (Integer pv : sigPath)
					{
						if (!(locked.contains(pv)))
						{
							locked.add(pv);
						}
					}
					match = true;
					targets.remove(ti);
					numMatches++;
				}
			}
			if (match && targets.size() == 0) 
			{
				// Return signal match...
				return RoundReturn.SIGNAL_MATCH;
			}

			// disp("Match? " + match);
			
			// Check to make sure we didn't go beyond the threshold
			if (x.xc > xorCap) return RoundReturn.XOR_EXCEEDED;

			// Update the graph accordingly to match the signal tree
			addAndConnectSignal(x, u, v);

			// Check to make sure we haven't gone over the cap on the number of maximum signals allowed
			// if (vertices.size() == numSignals - 1) return RoundReturn.SIGNALS_EXCEEDED;
		}
		return RoundReturn.SIGNAL_NOT_FOUND;
	}

	// TODO: update this function so that it matches the xorRound above
	public RoundReturn andRound(int numSignals, int andCap) throws Exception
	{
		Random r = new Random(System.currentTimeMillis());
		for (int i = 0; i < numSignals; i++)
		{
			Signal x = null;
			int u = 0;
			int v = 0;
			boolean repeat = false;
			int c = 0;
			do
			{
				if (repeat)
				{
					c++;
					// disp("Trying to find a different AND combination: " + c);
				}
				u = r.nextInt(vertices.size());
				v = r.nextInt(vertices.size());
				while (u == v) // enforce distinct signals
				{
					v = r.nextInt(vertices.size());
				}
				// disp("" + vertices.get(u).id);
				// disp("" + vertices.get(v).id);
				x = Signal.and(vertices.get(u), vertices.get(v));
				// disp("vertices.contains = " + vertices.contains(x));
				repeat = true;
			} 
			while (vertices.contains(x));

			boolean match = false;
			for (int ti = 0; ti < targets.size(); ti++)
			{
				int[] target = targets.get(ti);
				if (equal(x.signal, target)) 
				{
					disp("FOUND ONE SIGNAL.");
					disp(target);
					// System.exit(-1);
					addAndConnectSignal(x, u, v);

					// lock the signals along this path... 
					ArrayList<Integer> sigPath = signalPath(x);
					for (Integer pv : sigPath)
					{
						if (!(locked.contains(pv)))
						{
							locked.add(pv);
						}
					}
					match = true;
					targets.remove(ti);
					numMatches++;
				}
			}
			if (match && targets.size() == 0) 
			{
				return RoundReturn.SIGNAL_MATCH;
			}
			
			// Check to make sure we didn't go beyond the threshold
			if (x.ac > andCap) return RoundReturn.AND_EXCEEDED;

			// Update the graph accordingly to match the signal tree
			addAndConnectSignal(x, u, v);

			// Check to make sure we haven't gone over the cap on the number of maximum signals allowed
			// if (vertices.size() == numSignals - 1) return RoundReturn.SIGNALS_EXCEEDED;
		}
		return RoundReturn.SIGNAL_NOT_FOUND;
	}

	public void dropSignals(ArrayList<Integer> drops) throws Exception
	{
		// TODO: vertices = vertices - (drops - locked)
		for (Integer l : locked)
		{
			if (drops.contains(l))
			{
				drops.remove(l);
			}
		}
		for (Integer d : drops)
		{
			vertices.remove(d); 

			// sanity check with the edges
			for (int ei = 0; ei < edges.size(); ei++)
			{
				if (edges.get(ei).i == d)
				{
					// System.err.println("dropped vertex " + d + " is the tail of an edge to " + edges.get(ei).j);
					edges.remove(ei);
					// System.exit(-1);
				}
				if (edges.get(ei).j == d)
				{
					// System.err.println("dropped vertex " + d + " is the head of an edge from " + edges.get(ei).i);
					edges.remove(ei);
					// System.exit(-1);
				}
			}
		}
	}

	public HashMap<Integer, Integer> pathFans() throws Exception
	{
		HashMap<Integer, Integer> hp = new HashMap<Integer, Integer>();

		// Compute the weight of each signal path (degree of )
		for (Integer l : leaves)
		{
			if (l != 0)
			{
				int fo = 0; // fanout == out degree
				ArrayList<Signal> queue = new ArrayList<Signal>();
				queue.add(vertices.get(l));
				while (queue.size() > 0)
				{
					Signal t = queue.get(queue.size() - 1);
					queue.remove(queue.size() - 1);
					if (t != null)
					{
						queue.add(t.siga);
						// disp("adding id: " + t.siga.id);
						queue.add(t.sigb);
						// disp("Adding id: " + t.sigb.id);
						fo += fanout(t.id);
						// disp("Fanout for " + t.id + " = " + fanout(t.id));
					}
				}
				// disp("Leaf " + l + " fanout = " + fo);
				hp.put(l, fo);
			}
		}

		return hp;
	}

	public SLP optimize(long tries, long rounds, long seed, int pruneStrategy, int andCap, int xorCap, int maxDepth, int pruneFreq, double p_xor, int xc, double p_and, int ac) throws Exception
	{
		SLP slp = null;

		// TODO: how to handle this trivial case?
		// Trivial target - everything is zero.
		// if (isZero(target)) 
		// {
		// 	return null;
		// }

		ArrayList<int[]> cTargets = new ArrayList<int[]>();
		for (int[] t : targets)
		{
			cTargets.add(new int[t.length]);
			int ii = cTargets.size() - 1;
			for (int i = 0; i < t.length; i++)
			{
				cTargets.get(ii)[i] = t[i];
			}
		} 

		// Try to find the target t times
		for (long t = 0; t < tries; t++)
		{
			Random rng = new Random(seed);
			// disp("Try #" + t);

			// Reset! 
			int id = 0;
			vertices.clear();
			edges.clear();
			leaves.clear();
			locked.clear();
			ArrayList<Signal> initialSignals = initializeKnown(initial.size());
			for (Signal s : initialSignals) addSignal(s, id++);
			targets.clear();
			for (int[] tt : cTargets)
			{
				targets.add(new int[tt.length]);
				int ii = targets.size() - 1;
				for (int i = 0; i < tt.length; i++)
				{
					targets.get(ii)[i] = tt[i];
				}
			} 

			boolean failed = false;
			for (long r = 0; r < rounds && !failed; r++)
			{
				for (int i = 0; i < pruneFreq && !failed; i++)
				{
					// disp("XOR round: " + ((r * pruneFreq) + i));
					double p1 = rng.nextDouble();
					if (p1 < p_xor) 
					{
						RoundReturn res = xorRound(xc, xorCap);
						if (res == RoundReturn.SIGNAL_MATCH)
						{
							disp("BUILDING SLP");

							// TODO: build the SLP

							return new SLP(null, 0, 0);
						}
						else if (res == RoundReturn.XOR_EXCEEDED)
						{
							// disp("Failure: XORs exceeded");
							failed = true;
						}
					}
					double p2 = rng.nextDouble();
					if (p2 < p_and) 
					{
						RoundReturn res = andRound(ac, andCap);
						if (res == RoundReturn.SIGNAL_MATCH)
						{
							disp("BUILDING SLP");

							// TODO: build the SLP

							return new SLP(null, 0, 0);
						}
						else if (res == RoundReturn.AND_EXCEEDED)
						{
							// disp("Failure: ANDs exceeded");
							failed = true;
						}
					}
				}

				if (!failed)
				{
					// Prune algorithm:
					// 1. update distances
					// 2. compute heavy paths (and save them)
					// 3. prune the others (so long as they do not contain locked signals)
					// disp("Updating distances");

					switch (pruneStrategy)
					{
						case 0: // prune maximum weighted leaf trail
							updateDistances();
							HashMap<Integer, Integer> fatPaths = pathFans();
							// TODO: experiment with techniques - see what works better than the rest?
							// 1. above average?
							// 2. below average? 
							// 3. random droppage of half?
							int max = 0; // try #1: drop the vertices that have weight BELOW the max
							for (Integer v : fatPaths.values())
							{
								max = v > max ? v : max;
							}
							ArrayList<Integer> drops = new ArrayList<Integer>();
							for (Integer k : fatPaths.keySet())
							{
								if (fatPaths.get(k) == max) 
								{
									drops.add(k); // k == leaf ID
								}
							}

							// Drop the signals...
							dropSignals(drops);
							break;
						case 1:
						case 2:
							break;
					}
				}
			}
		}
		return slp;
	}

	public static void test() throws Exception
	{
		int[] signal = {0,1,0,1,1,1,1,1,1,0,0,0,1,0,0,1};
		int[] sig2   = {0,1,1,0,0,0,1,1,0,0,0,1,1,1,0,0};
		ArrayList<int[]> target = new ArrayList<int[]>();
		target.add(signal);
		target.add(sig2);
		ArrayList<Signal> known = initializeKnown(4);
		SignalGraph graph = new SignalGraph(known, target);
		// ArrayList<Signal> initialSignals, ArrayList<int[]> targets
		// SLP slp = graph.optimize(100000, 10, 123, 8, 20, 10, 1, 1.0, 2, 1.0, 1);
		SLP slp = graph.optimize(100000, 20, 123, 4, 10, 50, 50, 0.75, 2, 0.75, 1);
		// SLP optimize(long tries, long rounds, long seed, int andCap, int xorCap, int maxDepth, int pruneFreq, double p_xor, int xc, double p_and, int ac);
	}

	public static ArrayList<int[]> buildSignals(String file) throws Exception
	{
		ArrayList<int[]> signals = new ArrayList<int[]>();
		BufferedReader in = new BufferedReader(new FileReader(file)); 
		String line = "";
		while ((line = in.readLine()) != null && !line.isEmpty())
		{
			ArrayList<Integer> row = new ArrayList<Integer>();
			Scanner s = new Scanner(line);
			while (s.hasNext())
			{
				row.add(Integer.parseInt(s.next()));
			}
			int[] r = new int[row.size()];
			for (int i = 0; i < row.size(); i++) r[i] = row.get(i);
			signals.add(r);
		} 
		return signals;
	}

	public static void main(String[] args) throws Exception
	{
		if (args.length != 14)
		{
			test();
			System.err.println("usage: java SignalGraph nvars signalFile tries rounds seed pruneStrategy andCap xorCap maxDepth pruneFreq p_xor xc p_and ac");
			System.exit(-1);
		}

		// Parse the cmd line arguments
		String signalFile = args[1];
		long tries = Long.parseLong(args[2]);
		long rounds = Long.parseLong(args[3]);
		long seed = Long.parseLong(args[4]);
		int pruneStrategy = Integer.parseInt(args[5]);
		int andCap = Integer.parseInt(args[6]);
		int xorCap = Integer.parseInt(args[7]);
		int maxDepth = Integer.parseInt(args[8]);
		int pruneFreq = Integer.parseInt(args[9]);
		double p_xor = Double.parseDouble(args[10]);
		int xc = Integer.parseInt(args[11]);
		double p_and = Double.parseDouble(args[12]);
		int ac = Integer.parseInt(args[13]);

		// Open the file, which contains a list of signals, 
		// turn it into an arraylist of signals
		// run it, with the param list, through the optimizer
		// output the results (if successful)
		ArrayList<int[]> targets = buildSignals(signalFile);
		ArrayList<Signal> known = initializeKnown(Integer.parseInt(args[0]));
		SignalGraph graph = new SignalGraph(known, targets);

		disp("Trying to optimize:");
		for (int[] t : targets) disp(t);

		SLP slp = graph.optimize(tries, rounds, seed, pruneStrategy, andCap, xorCap, maxDepth, pruneFreq, p_xor, xc, p_and, ac);
		if (slp == null)
		{
			disp("Failed.");
		}
		else
		{
			disp("Passed.");
		}
	}
}