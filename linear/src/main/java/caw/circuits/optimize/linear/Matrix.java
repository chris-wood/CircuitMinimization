package caw.circuits.optimize.linear;

import java.util.ArrayList;
import java.math.BigInteger;

public class Matrix 
{
	private int[][] m;
	private int n;
	private int r;
	private int c;
	private int g;
	public String tag;

	// Empty constructor
	public Matrix(int nvars) throws Exception
	{
		this.r = 1;
		this.c = nvars;
		this.m = new int[1][nvars];
	}

	// Generic constructor
	public Matrix(int[][] mat) throws Exception
	{
		this.r = mat.length;    // assume > 0
		this.c = mat[0].length; // assume > 0
		this.n = 0; 

		if (r < 0 || c < 0) throw new Exception("Negative array dimension. Oh dear.");

		this.m = new int[r][c];
		for (int i = 0; i < r; i++)
		{
			for (int j = 0; j < c; j++)
			{
				if (mat[i][j] == 0 || mat[i][j] == 1)
				{
					this.m[i][j] = mat[i][j];
				}
				else
				{
					throw new Exception("Only binary matrices allowed.");
				}
			}
		}

		updateGateCount();
	}

	public Matrix(int[][] mat, int n) throws Exception
	{
		this.r = mat.length;    // assume > 0
		this.c = mat[0].length; // assume > 0
		this.n = n; // n == c to start

		if (r < 0 || c < 0) throw new Exception("Negative array dimension. Oh dear.");

		this.m = new int[r][c];
		for (int i = 0; i < r; i++)
		{
			for (int j = 0; j < c; j++)
			{
				if (mat[i][j] == 0 || mat[i][j] == 1)
				{
					this.m[i][j] = mat[i][j];
				}
				else
				{
					throw new Exception("Only binary matrices allowed.");
				}
			}
		}

		updateGateCount();
	}

	public Matrix(Matrix m) throws Exception
	{
		int[][] mat = m.getMatrix();
		this.r = mat.length;    // assume > 0
		this.c = mat[0].length; // assume > 0
		this.n = m.getN(); // n == c to start

		if (r < 0 || c < 0) throw new Exception("Negative array dimension. Oh dear.");

		this.m = new int[r][c];
		for (int i = 0; i < r; i++)
		{
			for (int j = 0; j < c; j++)
			{
				if (mat[i][j] == 0 || mat[i][j] == 1)
				{
					this.m[i][j] = mat[i][j];
				}
				else
				{
					throw new Exception("Only binary matrices allowed.");
				}
			}
		}

		updateGateCount();
	}

	public Matrix(ArrayList<int[]> rows) throws Exception
	{
		if (rows.size() == 0) throw new Exception("Empty set of rows is not allowed.");

		this.r = rows.size();
		this.c = rows.get(0).length;
		this.n = 0; 
		if (r < 0 || c < 0) throw new Exception("Negative array dimension. Oh dear.");
		this.m = new int[r][c];
		for (int i = 0; i < r; i++)
		{
			int[] row = rows.get(i);
			for (int j = 0; j < c; j++)
			{	
				if (row[j] == 0 || row[j] == 1)
				{
					this.m[i][j] = row[j];
				} 
				else
				{
					throw new Exception("Only binary matrices allowed.");
				}
			}
		}

		updateGateCount();
	}

	public Matrix(ArrayList<int[]> rows, int n) throws Exception
	{
		if (rows.size() == 0) throw new Exception("Empty set of rows is not allowed.");

		this.r = rows.size();
		this.c = rows.get(0).length;
		this.n = n; // THIS IS NOT USED FOR LINEAR OPTIMIZATION
		if (r < 0 || c < 0) throw new Exception("Negative array dimension. Oh dear.");
		this.m = new int[r][c];
		for (int i = 0; i < r; i++)
		{
			int[] row = rows.get(i);
			for (int j = 0; j < c; j++)
			{	
				if (row[j] == 0 || row[j] == 1)
				{
					this.m[i][j] = row[j];
				} 
				else
				{
					throw new Exception("Only binary matrices allowed.");
				}
			}
		}

		updateGateCount();
	}

	public Matrix copy() throws Exception
	{
		int[][] matrix = new int[r][c];
		for (int i = 0; i < r; i++) 
		{
			for (int j = 0; j < c; j++)
			{
				matrix[i][j] = m[i][j];
			}
		}
		return new Matrix(matrix, n);
	}

	public boolean isZero()
	{
		for (int i = 0; i < r; i++)
		{
			for (int j = 0; j < c; j++)
			{
				if (m[i][j] != 0) return false;
			}
		}
		return true;
	}

	public int[][] getMatrix() { return m; }
	public int getGateCount() { return g; }
	public int getDimension() { return r; }
	public int getLength() { return c; }
	public int getN() { return n; }
	public int getRowCount() { return r; }
	public int getColCount() { return c; }
	public ArrayList<Vector> getVectors() throws Exception
	{
		ArrayList<Vector> vectors = new ArrayList<Vector>();
		for (int i = 0; i < r; i++)
		{
			vectors.add(new Vector(this.getRow(i), i));
		}
		return vectors;
	}

	public ArrayList<int[]> getRows() throws Exception
	{
		ArrayList<int[]> rows = new ArrayList<int[]>();
		for (int i = 0; i < r; i++)
		{
			rows.add(this.getRow(i));
		}
		return rows;
	}

	public Matrix subMatrix(int nvars) throws Exception
	{
		return subMatrix(0, r, 0, nvars);
	}

	public Matrix subMatrix(int rs, int re, int cs, int ce) throws Exception
	{
		ArrayList<int[]> rows = new ArrayList<int[]>();
		for (int i = rs; i < re; i++)
		{
			int[] row = new int[re - rs + 1];
			for (int j = cs; j < ce; j++)
			{
				row[i] = m[i][j];
			}
			rows.add(row);
		}
		return new Matrix(rows);
	}

	private void updateGateCount() 
	{
		g = 0;

		// Count the number of new variables
		g += c - n; 

		// Count the gates contributed by the original variables
		for (int i = 0; i < r; i++)
		{
			int ng = -1;
			for (int j = 0; j < c; j++)
			{
				if (m[i][j] == 1) ng++;
			}
			if (ng != -1) g += ng;
		}
	}

	public int[] xorRows(int[] rowIndices) throws Exception
	{
		int[] result = new int[c];
		for (int i = 0; i < rowIndices.length; i++)
		{
			if (rowIndices[i] != 0)
			{
				result = MatrixOptimize.XOR(result, this.getRow(i));
			}
		}
		return result;
	}

	public int[] xorRows(BigInteger indices) throws Exception
	{
		int[] result = new int[c];
		for (int i = 0; i < r; i++)
		{
			if (indices.testBit(i))
			{
				result = MatrixOptimize.XOR(result, this.getRow(i));
			}
		}
		return result;
	}

	public int[] xorRows(int indices) throws Exception
	{
		int[] result = new int[c];
		for (int i = 0; i < r; i++)
		{
			if ((indices & (1 << i)) > 0)
			{
				result = MatrixOptimize.XOR(result, this.getRow(i));
			}
		}
		return result;
	}

	public void setEntry(int i, int j, int v) throws Exception
	{
		if ((i < 0 || i >= r) || (j < 0 || j >= c))
		{
			throw new Exception("Invalid column or row index.");
		}
		m[i][j] = v;
		updateGateCount();
	}

	public int[] getColumn(int index) throws Exception
	{
		if (index < 0 || index >= c)
		{
			throw new Exception("Invalid column index: " + index);
		}

		int[] col = new int[r];
		for (int i = 0; i < r; i++) {
			col[i] = m[i][index];
		}
		return col;
	}

	public int[] getRow(int index) throws Exception
	{
		if (index < 0 || index >= r)
		{
			throw new Exception("Invalid row index: " + index);
		}

		int[] row = new int[c];
		for (int i = 0; i < c; i++) {
			row[i] = m[index][i];
		}
		return row;
	}

	public void putColumn(int[] col, int index) throws Exception
	{
		if (col.length != r)
		{
			throw new Exception("Invalid column dimension.");
		}
		if (index < 0 || index >= c) 
		{
			throw new Exception("Invalid column index.");
		}

		// Drop in the new column values (replace old values...)
		for (int i = 0; i < m.length; i++) 
		{
			m[i][index] = col[i];
		}
		updateGateCount();
	}

	public void putRow(int[] row, int index) throws Exception
	{
		if (row.length != c)
		{
			throw new Exception("Invalid row dimension.");
		}
		if (index < 0 || index >= r) 
		{
			throw new Exception("Invalid row index.");
		}

		// Drop in the new column values (replace old values...)
		for (int i = 0; i < m[0].length; i++) 
		{
			m[index][i] = row[i];
		}
		updateGateCount();
	}

	public void appendColumn(int[] col) throws Exception
	{
		if (col.length != r)
		{
			throw new Exception("Invalid column dimension.");
		}

		// Allocate memory for the new matrix and append the new last column
		c++;
		int[][] newM = new int[r][c];
		for (int i = 0; i < m.length; i++) 
		{
			for (int j = 0; j < m[i].length; j++) 
			{
				newM[i][j] = m[i][j];
			}
		}
		for (int i = 0; i < r; i++) {
			newM[i][c - 1] = col[i];
		}

		// Swap references, let the GC take away the old m
		this.m = newM;
		updateGateCount();
	}

	public void appendRow(int[] row) throws Exception
	{
		if (row.length != c)
		{
			throw new Exception("Invalid row dimension.");
		}

		// Allocate memory for the new matrix and append the new last column
		r++;
		int[][] newM = new int[r][c];
		for (int i = 0; i < m.length; i++) 
		{
			for (int j = 0; j < m[i].length; j++) 
			{
				newM[i][j] = m[i][j];
			}
		}
		for (int i = 0; i < c; i++) {
			newM[r - 1][i] = row[i];
		}

		// Swap references, let the GC take away the old m
		this.m = newM;
		updateGateCount();
	}

	public boolean containsRow(int[] row) throws Exception
	{
		if (row.length != c)
		{
			throw new Exception("Invalid row dimension.");
		}
		for (int i = 0; i < m.length; i++)
		{
			int[] r_i = m[i];
			boolean match = true;
			for (int j = 0; j < r_i.length; j++)
			{
				if (r_i[j] != row[j]) 
				{
					match = false;
					break;
				}
			}
			if (match) return true;
		}
		return false;
	}

	public String toString()
	{
		String result = "";
		for (int i = 0; i < r; i++) 
		{
			result += "[";
			for (int j = 0; j < c; j++) 
			{
				result += "" + m[i][j] + " ";
			}
			result = result.substring(0, result.length() - 1) + "]\n";
		}
		return result;
	}
}
