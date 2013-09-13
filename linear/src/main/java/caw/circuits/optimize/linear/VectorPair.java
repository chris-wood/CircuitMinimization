package caw.circuits.optimize.linear;

import java.util.Arrays;

public class VectorPair
{
	public int[] i;
	public int[] j;
	public VectorPair(int[] ii, int[] jj) {this.i = Arrays.copyOfRange(ii, 0, ii.length); this.j = Arrays.copyOfRange(jj, 0, jj.length);}

	@Override
	public boolean equals(Object other){
		if (other == null) return false;
		if (other == this) return true;
		if (!(other instanceof VectorPair))return false;
		VectorPair op = (VectorPair)other;
		if (Arrays.equals(i, op.i) && Arrays.equals(j, op.j)) return true;
		return false;
	}

	@Override
	public int hashCode()
	{
		return (1000 * i.hashCode()) + j.hashCode();
	}
}