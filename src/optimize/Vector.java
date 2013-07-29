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
}