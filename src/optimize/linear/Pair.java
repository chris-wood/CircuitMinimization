public class Pair
{
	public int i;
	public int j;
	public Pair(int ii, int jj) {this.i = ii; this.j = jj;}

	@Override
	public boolean equals(Object other){
		if (other == null) return false;
		if (other == this) return true;
		if (!(other instanceof Pair))return false;
		Pair op = (Pair)other;
		if (i == op.i && j == op.j) return true;
		return false;
	}

	@Override
	public int hashCode()
	{
		return (1000 * i) + j;
	}
}