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
