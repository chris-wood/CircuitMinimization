public class BasePair
{
	public int[] base;
	public int[] dist;
	public Pair p;

	public BasePair(int[] b, int[] d, Pair p)
	{
		this.p = new Pair(p.i, p.j);
		this.base = new int[b.length];
		this.dist = new int[d.length];
		for (int i = 0; i < b.length; i++) this.base[i] = b[i];
		for (int i = 0; i < d.length; i++) this.dist[i] = d[i];
	}
}