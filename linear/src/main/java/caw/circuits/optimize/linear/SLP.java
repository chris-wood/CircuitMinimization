package caw.circuits.optimize.linear;

import java.util.ArrayList;

public class SLP
{
	public ArrayList<String> lines;
	public int xc;
	public int ac;

	public SLP(ArrayList<String> lines, int xc, int ac)
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