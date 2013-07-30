import java.util.ArrayList;

public class Signal
{
	public enum SignalOp { AND, XOR };

	public static String GLOBAL_VAR_STRING = "X";
	public static int GLOBAL_VAR_ID = 0;

	public SignalOp op;
	public int id;
	public int ac;
	public int xc;
	public int[] signal;
	public Signal siga;
	public Signal sigb;
	public String rep;
	public int depth;

	public Signal(int[] signal, int ac, int xc, Signal siga, Signal sigb, SignalOp op, String rep, int depth)
	{
		GLOBAL_VAR_ID++;
		this.id = GLOBAL_VAR_ID;
		this.depth = depth;
		this.rep = rep;
		this.op = op;
		this.siga = copy(siga);
		this.sigb = copy(sigb);
		this.ac = ac;
		this.xc = xc;
		this.signal = new int[signal.length];
		for (int i = 0; i < signal.length; i++) this.signal[i] = signal[i];
	}

	public Signal(int[] signal, int id)
	{
		GLOBAL_VAR_ID++;
		this.id = GLOBAL_VAR_ID;
		this.depth = 0;
		this.id = id;
		this.rep = "" + id;
		this.siga = null;
		this.sigb = null;
		this.ac = 0;
		this.xc = 0;
		this.signal = new int[signal.length];
		for (int i = 0; i < signal.length; i++) this.signal[i] = signal[i];	
	}
	
	public static Signal xor(Signal u, Signal v) throws Exception
	{
		int[] result = signalXOR(u.signal, v.signal);
		int d = u.depth > v.depth ? u.depth : v.depth;
		Signal s = new Signal(result, u.ac + v.ac, u.xc + v.xc + 1, u, v, SignalOp.XOR, "(" + u.rep + ") + (" + v.rep + ")", d + 1);
		GLOBAL_VAR_ID++;
		s.id = GLOBAL_VAR_ID;
		return s;
	}

	public static Signal and(Signal u, Signal v) throws Exception
	{
		int[] result = signalAND(u.signal, v.signal);
		int d = u.depth > v.depth ? u.depth : v.depth;
		Signal s = new Signal(result, u.ac + v.ac + 1, u.xc + v.xc, u, v, SignalOp.AND, "(" + u.rep + ")(" + v.rep + ")", d + 1);
		GLOBAL_VAR_ID++;
		s.id = GLOBAL_VAR_ID;
		return s;
	}

	public Signal copy(Signal s) // unoptimized deep copy
	{
		if (s == null) return null;
		int[] sc = new int[s.signal.length];
		for (int i = 0; i < sc.length; i++) sc[i] = s.signal[i];
		Signal newSig = new Signal(sc, s.ac, s.xc, copy(s.siga), copy(s.sigb), s.op, s.rep, s.depth);
		// newSig.id = s.id;
		GLOBAL_VAR_ID++;
		newSig.id = GLOBAL_VAR_ID;
		return newSig;
	}

	public static int[] signalXOR(int[] u, int[] v) throws Exception
	{
		if (u.length != v.length) throw new Exception("Vector dimensions don't match.");
		int[] x = new int[u.length];
		for (int i = 0; i < u.length; i++) x[i] = u[i] ^ v[i];
		return x;
	}

	public static int[] signalAND(int[] u, int[] v) throws Exception
	{
		if (u.length != v.length) throw new Exception("Vector dimensions don't match.");
		int[] x = new int[u.length];
		for (int i = 0; i < u.length; i++) x[i] = u[i] * v[i];
		return x;
	}

	public ArrayList<String> toSLP()
	{
		ArrayList<String> slp = new ArrayList<String>();

		if (siga == null && sigb == null)
		{	
			slp.add(0, GLOBAL_VAR_STRING + id);
		}
		else if (siga == null)
		{
			ArrayList<String> slpA = siga.toSLP();
			for (int i = slpA.size() - 1; i >= 0; i--)
			{
				slp.add(0, "Y" + siga.id + " <= " + slpA.get(i));
			}
		}
		else if (sigb == null)
		{
			ArrayList<String> slpB = sigb.toSLP();
			for (int i = slpB.size() - 1; i >= 0; i--)
			{
				slp.add(0, "Y" + sigb.id + " <= " + slpB.get(i));
			}
		}
		else
		{
			ArrayList<String> slpA = siga.toSLP();
			for (int i = slpA.size() - 1; i >= 0; i--)
			{
				slp.add(0, "Y" + siga.id + " <= " + slpA.get(i));
			}
			ArrayList<String> slpB = sigb.toSLP();
			for (int i = slpB.size() - 1; i >= 0; i--)
			{
				slp.add(0, "Y" + sigb.id + " <= " + slpB.get(i));
			}
		}

		return slp;
	}

	public ArrayList<Signal> decompose()
	{
		ArrayList<Signal> sigs = new ArrayList<Signal>();
		sigs.add(this);
		if (siga == null && sigb == null)
		{	
			return sigs;
		}
		else if (siga == null)
		{
			sigs.addAll(sigb.decompose());
			return sigs;
		}
		else if (sigb == null)
		{
			sigs.addAll(siga.decompose());
			return sigs;
		}
		else
		{
			sigs.addAll(siga.decompose());
			sigs.addAll(sigb.decompose());
			return sigs;
		}
	}

	@Override
	public boolean equals(Object other)
	{
		if (other == null) return false;
		if (other == this) return true;
		if (!(other instanceof Signal))return false;
		Signal op = (Signal)other;
		for (int i = 0; i < signal.length; i++)
		{
			if (signal[i] != op.signal[i]) return false;
		}
		return true;
	}

	@Override
	public int hashCode()
	{
		int s = 1;
		int sum = 0;
		for (int i = 0; i < signal.length; i++)
		{
			s *= 2;
			if (signal[i] == 1) sum += s;
		}
		return sum;
	}
}

