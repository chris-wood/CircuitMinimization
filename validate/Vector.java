package caw.circuits.optimize.linear;

import java.math.BigInteger;

public class Vector
{
	public int id;
	public int[] row;
	public Vector(int[] r, int id)
	{
		this.id = id;
		row = new int[r.length];
		for (int i = 0; i < r.length; i++) row[i] = r[i];
	}

	public Vector(int[] r)
	{
		this.row = new int[r.length];
		for (int i = 0; i < r.length; i++) row[i] = r[i];
	}

	public Vector(String bits)
	{
		BigInteger bi = new BigInteger(bits, 2);
		this.row = new int[bits.length()];
		for (int i = 0; i < this.row.length; i++) row[i] = bi.testBit(i) ? 1 : 0;
	}

	@Override
	public boolean equals(Object o)
	{
		if (o instanceof Vector)
		{
			Vector ov = (Vector)o;
			for (int i = 0; i < ov.row.length; i++)
			{
				if (row[i] != ov.row[i])
				{	
					return false;
				}
			}
		}
		else
		{
			return false;
		}
		return true;
	}
}
