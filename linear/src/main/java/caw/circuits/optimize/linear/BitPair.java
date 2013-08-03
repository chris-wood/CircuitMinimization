package caw.circuits.optimize.linear;

public class BitPair 
{
	int n;
	int k;
	MatrixOptimize.WALK_DIR dir;
	String ls;
	String pre;
	String post;
	public BitPair(int n, int k, String pre, String post, MatrixOptimize.WALK_DIR dir, String ls) 
	{
		this.n = n; this.k = k; this.pre = pre; this.post = post; this.dir = dir; this.ls = ls; 
	}
}