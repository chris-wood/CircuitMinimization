package caw.circuits.optimize.linear;

public class SLPLine
{
	// Vars for the line
	public SLPElement out;
	public SLPElement o1;
	public SLPElement o2;
	public SLP.SLPOperator op;

	public SLPLine(SLPElement out, SLPElement o1, SLPElement o2, SLP.SLPOperator op)
	{
		this.out = out;
		this.o1 = o1;
		this.o2 = o2;
		this.op = op;
	}

	@Override
	public String toString()
	{
		return out.name + " = " + o1.name + " " + SLP.SLPOperatorSymbol[op.ordinal()] + " " + o2.name; 
	}
}