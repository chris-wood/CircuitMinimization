package caw.circuits.optimize.linear;

import java.util.ArrayList;

public class SolutionCircuit
{
	public int d = Integer.MAX_VALUE;
	public ArrayList<int[]> optimalBases = new ArrayList<int[]>();
	public ArrayList<int[]> optimalDistances = new ArrayList<int[]>();
	public ArrayList<Pair> pairs = new ArrayList<Pair>();

	public void reset()
	{
		d = Integer.MAX_VALUE;
		optimalBases = new ArrayList<int[]>();
		optimalDistances = new ArrayList<int[]>();
		pairs = new ArrayList<Pair>();
	}

	public synchronized void record(int t_d, ArrayList<int[]> t_optimalBases, ArrayList<int[]> t_optimalDistances, ArrayList<Pair> t_pairs)
	{
		if (t_d < d)
		{
			this.d = t_d;
			this.optimalBases = t_optimalBases;
			this.optimalDistances = t_optimalDistances;
			this.pairs = t_pairs;
		}
		else if (t_d == d) // add to the existing solution set...
		{
			for (int[] b : t_optimalBases)
			{
				this.optimalBases.add(b);
			}
			for (int[] d : t_optimalDistances)
			{
				this.optimalDistances.add(d);
			}
			for (Pair p : t_pairs)
			{
				this.pairs.add(p);
			}
		}
	}
}