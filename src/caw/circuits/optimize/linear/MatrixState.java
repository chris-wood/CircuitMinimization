public class MatrixState
{
	public Matrix matrix; // just a 1 or 0 in each spot...
	public int i; // new variable share
	public int j; // new variable share
	public int var; // the new variable that was introduced
	public Pair pair;

	public MatrixState(Matrix m, int ii, int jj, int newVar)
	{
		matrix = m;
		i = ii;
		j = jj;
		var = newVar;
	}

	public MatrixState(Matrix m, Pair p)  
	{
		matrix = m;
		pair = p;
		this.i = p.i;
		this.j = p.j;
	}
}