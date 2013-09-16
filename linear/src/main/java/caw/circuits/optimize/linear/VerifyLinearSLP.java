package caw.circuits.optimize.linear;

import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.FileReader;

public class VerifyLinearSLP
{
	Matrix mat;
	SLP slp;

	public static void main(String[] args)
	{
		// cmd line arguments:
		// two file names, one for the matrix and one for the SLP
		if (args.length != 2)
		{
			System.err.println("usage: java -jar VerifyLinearSLP matrixFile slpFile");
			System.exit(-1);
		}

		// parse the matrix file and slp file to populate the data
		// TODO: currently rewriting the SLP representation code...
		// TODO: i'll come back to that later I guess...
	}
}