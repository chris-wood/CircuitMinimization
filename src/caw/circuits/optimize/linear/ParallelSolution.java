public class ParallelSolution 
{
	public int gOptimal;
	public Matrix currMatrix;
	public int oldi;
	public int oldj;
	public int n;
	public int lastcol;

	public ParallelSolution(int go, Matrix cm, int oi, int oj, int nn, int lc) throws Exception
	{
		this.gOptimal = go;
		this.currMatrix = cm.copy();
		this.oldi = oi;
		this.oldj = oj;
		this.n = nn;
		this.lastcol = lc;
	}
}