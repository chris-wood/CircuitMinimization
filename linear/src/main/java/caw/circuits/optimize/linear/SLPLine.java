package caw.circuits.optimize.linear;

public class SLPLine
{
	public static enum SLPOperator { SLP_ADD, SLP_MULT };
	public static char[] SLPOperatorSymbol = {'+', 'x'};

	// Vars for the line
	public String out;
	public String o1;
	public String o2;
	public SLPOperator op;

	public SLPLine(String out, String o1, String o2, SLPOperator op)
	{
		this.out = new String(out);
		this.o1 = new String(o1);
		this.o2 = new String(o2);
		this.op = op;
	}

	@Override
	public String toString()
	{
		// return out + " = " + o1 + " " + SLPOperatorSymbol[op] + " " + o2; 
		return "DERP";
	}
}