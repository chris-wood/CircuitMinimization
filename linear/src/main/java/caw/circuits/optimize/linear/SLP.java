package caw.circuits.optimize.linear;

import java.util.ArrayList;
import java.util.Collections;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;

public class SLP
{
	public enum SLPOperator { SLP_ADD, SLP_MULT, SLP_XNOR, SLP_OP_INVALID };
	public static char[] SLPOperatorSymbol = {'+', 'x', '#'};

	// composed of SLP lines, which are SLP elements with an operand

	public ArrayList<SLPLine> lines;
	// public ArrayList<SLPElement> ins;
	// public ArrayList<SLPElement> vars;
	// public ArrayList<SLPElement> outs;

	public SLP(ArrayList<SLPLine> lines)
	{
		this.lines = lines;
		// this.lines = new ArrayList<SLPLine>();
		// for (SLPLine l : lines) this.lines.add(l);
		// this.ins = ins;
		// this.ins = new ArrayList<SLPElement>();
		// for (SLPElement e : ins) this.ins.add(e);
		// this.vars = vars;
		// this.vars = new ArrayList<SLPElement>();
		// for (SLPElement e : vars) this.vars.add(e);
		// this.outs = outs;
		// this.outs = new ArrayList<SLPElement>();
		// for (SLPElement e : outs) this.outs.add(e);
	}

	// SLP bitwise operations
	public int BIT_ADD(int v1, int v2)
	{
		return (v1 ^ v2);
	}
	public int BIT_MULT(int v1, int v2)
	{
		return (v1 & v2);
	}
	public int BIT_XNOR(int v1, int v2)
	{
		if (v1 == v2) return 1;
		else return 0;
	}
	public int BIT_NEGATE(int v1)
	{
		if (v1 == 1) return 0;
		else return 1;
	}

	public void load_inputs(ArrayList<SLPElement> inputs)
	{
		for (int i = 0; i < lines.size(); i++)
		{
			for (SLPElement e : inputs)
			{
				if (lines.get(i).o1.name.equals(e.name))
				{
					lines.get(i).o1.val = e.val;
				}
				if (lines.get(i).o2.name.equals(e.name))
				{
					lines.get(i).o2.val = e.val;
				}
			}
		}
	}

	public ArrayList<Integer> extract_outputs(ArrayList<String> outputs)
	{
		ArrayList<Integer> outputValues = new ArrayList<Integer>();
		for (int i = 0; i < lines.size(); i++)
		{
			for (String e : outputs)
			{
				if (lines.get(i).out.name.equals(e))
				{
					outputValues.add(lines.get(i).out.val);
				}
			}
		}
		return outputValues;
	}

	public void execute() throws Exception
	{
		// debug
		System.err.println("Executing the SLP of size: " + lines.size());

		try
		{
			for (int i = 0; i < lines.size(); i++)
			{
				SLPLine line = lines.get(i);

				// Check for negation special case first
				if (line.o1.name.equals("ONE"))
				{
					line.out.val = BIT_NEGATE(line.o2.val);
				}
				else if (line.o2.name.equals("ONE"))
				{
					line.out.val = BIT_NEGATE(line.o1.val);
				}
				else
				{
					switch (line.op)
					{
						case SLP_ADD:
							line.out.val = BIT_ADD(line.o1.val, line.o2.val);
							break;

						case SLP_MULT:
							line.out.val = BIT_MULT(line.o1.val, line.o2.val);
							break;

						case SLP_XNOR:
							line.out.val = BIT_XNOR(line.o1.val, line.o2.val);
							break;
						default:
							throw new Exception("Error: invalid SLP operator encountered: " + line.op);
					}
				}
			}
		} 
		catch (Exception e)
		{
			System.err.println(e.getMessage());
			e.printStackTrace();
			throw e;
		}
	}

	public static SLPOperator find_matching_operator(String sym)
	{
		// debug
		System.err.println("Finding operator match for: " + sym);

		if (sym.equals("" + SLPOperatorSymbol[0]))
		{
			return SLPOperator.SLP_ADD;
		}
		else if (sym.equals("" + SLPOperatorSymbol[1]))
		{
			return SLPOperator.SLP_MULT;
		}
		else if (sym.equals("" + SLPOperatorSymbol[2]))
		{
			return SLPOperator.SLP_XNOR;
		}
		else
		{
			return SLPOperator.SLP_OP_INVALID;
		}
	}

	@Override
	public boolean equals(Object o)
	{
		if (o instanceof SLP)
		{
			SLP other = (SLP)o;
			if (other.lines.size() != lines.size())
			{
				return false;
			}
			for (int i = 0; i < lines.size(); i++)
			{
				if (!(lines.get(i).equals(other.lines.get(i))))
				{
					return false;
				}
			}
			return true;
		}
		else
		{
			return false;
		}
	}

	@Override
	public int hashCode()
	{
		return lines.hashCode();
	}

	@Override
	public String toString()
	{
		return "TODO";
	}
}