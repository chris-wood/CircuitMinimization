package caw.circuits.optimize.linear;

import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.FileReader;

public class VerifyLinearSLP
{

	public static boolean all_zeros(int[] v)
	{
		for (int i = 0; i < v.length; i++) 
		{
			if (v[i] == 1) return false;
		}
		return true;
	}

	public static String str(int[] v)
	{
		String s = "[ ";
		for (int i = 0; i < v.length; i++) 
		{
			s = s + v[i] + " ";
		}
		s = s + "]";
		return s;
	}

	public static boolean are_equal(int[] v1, int[] v2) throws Exception
	{
		if (v1.length != v2.length) throw new Exception("Vector dimensions don't match");
		for (int i = 0; i < v1.length; i++)
		{
			if (v1[i] != v2[i]) return false;
		}
		return true;
	}

	public static void main(String[] args) throws Exception
	{
		// cmd line arguments:
		// two file names, one for the matrix and one for the SLP
		if (args.length != 2)
		{
			System.err.println("usage: java -jar VerifyLinearSLP matrixFile slpFile");
			System.exit(-1);
		}

		// Parse the matrix and SLP
		Matrix mat = Matrix.parseMatrixFile(args[0]);
		SLPParser parser = new SLPParser();
		parser.parse(args[1], false);

		// Check that #rows = #outputs and #columns = #inputs
		if (mat.getLength() != parser.inputs.size())
		{
			System.out.println("Input variables do not match.");
			System.exit(-1);
		}
		if (mat.getDimension() != parser.outputs.size())
		{
			System.out.println("Output variables do not match.");	
			System.exit(-1);
		}

		// Generate all 2^n combinations of inputs, execute SLP and perform matrix product, and compare result
		int[] inputs = new int[parser.inputs.size()]; // all 0s to start
		do
		{	
			// Compute the output from this input
			System.out.println("Computing: " + str(inputs));
			int[] slpOut = parser.execute(inputs);
			int[] matOut = mat.multiply(inputs);

			if (!(are_equal(slpOut, matOut)))
			{
				ArrayList<Integer> iv = new ArrayList<Integer>();
				for (int i = 0; i < inputs.length; i++) iv.add(inputs[i]);
				System.out.println("SLP and Matrix don't match on input: " + iv.toString());
				System.out.println("SLP: " + str(slpOut));
				System.out.println("MAT: " + str(matOut));
				System.exit(-1);
			}

			// Advance the input...
			for (int i = 0; i < inputs.length; i++)
			{
				if (inputs[i] == 0)
				{
					inputs[i] = 1;
					break;
				}
				else
				{
					inputs[i] = 0;
				}
			}
		} while(!(all_zeros(inputs)));

		// made it this far... they must compute the same result!
		System.out.println("All inputs match. SLP and Matrix compute the same value.");
	}
}