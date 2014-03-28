package caw.circuits.optimize.linear;

import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;

public class SLPStrings
{
	public ArrayList<String> lines;
	public ArrayList<SLPLine> loadedLines;
	public int inputVars;
	public int outputVars;
	public int xc;
	public int ac;

	public SLPStrings(ArrayList<String> lines, int xc, int ac)
	{
		this.lines = lines;
		this.xc = xc;
		this.ac = ac;
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