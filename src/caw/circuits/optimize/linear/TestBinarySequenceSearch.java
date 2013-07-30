import java.util.*;

public class TestBinarySequenceSearch
{
	public static void main(String args[]) throws Exception
	{
		if (args.length != 2) 
		{
			System.err.println("usage: java TestBinarySequenceSearch n k");
			System.exit(-1);
		}

		int n = Integer.parseInt(args[0]);
		int k = Integer.parseInt(args[1]);

		ArrayList<int[]> rows = new ArrayList<int[]>();
		for (int i = 0; i < n; i++)
		{
			rows.add(new int[10]);
		}

		// Run both procedures
		int[] f = new int[n];
		long start1 = System.currentTimeMillis();
		ArrayList<Vector> indices = MatrixOptimize.walkWeightedSequences(n, k);
		for (Vector v : indices)
		{
			int[] test = MatrixOptimize.XOR(f, v.row);
		}
		long end1 = System.currentTimeMillis();

		optimizedPeraltaDistance
	}
}