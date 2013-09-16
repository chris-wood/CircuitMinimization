package caw.circuits.optimize.linear;

import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;

public class SLP
{
	public ArrayList<String> lines;
	public ArrayList<SLPLine> loadedLines;
	public int inputVars;
	public int outputVars;
	public int xc;
	public int ac;

	public SLP(ArrayList<String> lines, int xc, int ac)
	{
		this.lines = lines;
		this.xc = xc;
		this.ac = ac;
	}

	// Parse an SLP file, which is composed of an arbitrary number of lines
	// x chars denote input vector elements
	// t chars denote temp variables (bits)
	// y chars denote output vector elements
	public static void parseSlpFile(String slpFile) throws IOException, FileNotFoundException
	{
		ArrayList<Matrix> matrices = new ArrayList<Matrix>();
		BufferedReader in = new BufferedReader(new FileReader(slpFile)); 
		int numInputs = 0;
		int numOutputs = 0;
		String line = "";

		// Loop until we reach the end
		while ((line = in.readLine()) != null && !line.isEmpty())
		{
			// line: out = op1 OP op2
			// check out, op1, op2 to determine what type of variable

			// TODO: want to make this code a little better... try and represent SLPs as DAGs using composite pattern
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
		return lines.hashCode() + xc + ac;
	}

	@Override
	public String toString()
	{
		StringBuilder builder = new StringBuilder();
		for (String s : lines)
		{
			builder.append(s);
		}
		return builder.toString();
	}
}