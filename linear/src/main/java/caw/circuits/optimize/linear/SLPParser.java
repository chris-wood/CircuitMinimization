package caw.circuits.optimize.linear;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;

public class SLPParser
{
	public ArrayList<SLPElement> inputs = new ArrayList<SLPElement>();
	public ArrayList<SLPElement> outputs = new ArrayList<SLPElement>();
	public ArrayList<SLPLine> lines = new ArrayList<SLPLine>();
	public ArrayList<SLPElement> vars = new ArrayList<SLPElement>();
	public int numGates;

	public class MatchResult
	{
		SLPElement element;
		boolean newVariable;

		public MatchResult(SLPElement e, boolean nv)
		{
			element = e;
			newVariable = nv;
		}
	}

	public MatchResult find_matching_element(String id)
	{
		SLPElement match = null;

		for (SLPElement out : outputs) // check against outputs
		{
			if (id.equals(out.name))
			{
				match = out;
				break;
			}
		}
		if (match == null)
		{
			for (SLPElement in : inputs) // check against inputs
			{
				if (id.equals(in.name))
				{
					match = in;
					break;
				}
			}

			if (match == null)
			{
				for (SLPElement var : vars) // check against variables
				{
					if (id.equals(var.name))
					{
						match = var;
						break;
					}	
				}
			}
		}

		// If null at this point, then we must have a new variable, so return it
		boolean newVar = false;
		if (match == null)
		{
			newVar = true;
			match = new SLPElement(id, 0);
		}

		return new MatchResult(match, newVar);
	}

	public int[] execute(int[] ins) throws Exception
	{
		// Drop in the inputs so that the circuit can be executed
		for (int i = 0; i < inputs.size(); i++)
		{
			inputs.get(i).val = ins[i];
		}

		//// LOAD INPUTS AND BEGIN EXECUTION
		SLP slp = new SLP(lines);
		slp.load_inputs(inputs);
		slp.execute();
		//// END EXECUTION

		// Pull the output
		ArrayList<String> outputNames = new ArrayList<String>();
		for (int i = 0; i < outputs.size(); i++)
		{
			outputNames.add(outputs.get(i).name);
		}
		ArrayList<Integer> outputVals = slp.extract_outputs(outputNames);
		int[] outs = new int[outputVals.size()];
		for (int i = 0; i < outs.length; i++)
		{
			outs[i] = outputVals.get(i);
		}

		return outs;
	}

	public void parse(String fileName, boolean execute) throws Exception
	{
		BufferedReader in = new BufferedReader(new FileReader(fileName));

		// debug
		System.err.println("Now parsing the SLP metadata...");

		// Read the SLP metadata
		String line = in.readLine();
		String[] data = line.split(" ");
		numGates = Integer.parseInt(data[0]);

		// Read SLP inputs
		line = in.readLine();
		data = line.split(" ");
		int numInputs = Integer.parseInt(data[0]);
		line = in.readLine();
		data = line.split(" ");
		for (String d : data)
		{
			inputs.add(new SLPElement(d, 0));
		}
		System.err.println("Input signals: " + line);

		// Read SLP outputs
		line = in.readLine();
		data = line.split(" ");
		int numOutputs = Integer.parseInt(data[0]);
		line = in.readLine();
		data = line.split(" ");
		for (String d : data)
		{
			outputs.add(new SLPElement(d, 0));
		}
		System.err.println("Output signals: " + line);

		// Skip "begin"
		line = in.readLine();

		// debug
		System.err.println("Now parsing the SLP circuit lines...");
		
		// Now parse the lines of the SLP
		for (int i = 0; i < numGates; i++)
		{
			line = in.readLine();
			Scanner dataIn = new Scanner(line);

			// debug
			System.err.println("Parsing: " + line);

			String output = dataIn.next();
			SLPElement slpOutput = null;
			dataIn.next();
			String o1 = dataIn.next();
			SLPElement slpO1 = null;
			String op = dataIn.next();
			SLP.SLPOperator slpOp = SLP.find_matching_operator(op);
			String o2 = dataIn.next();
			SLPElement slpO2 = null;

			// Find matches to build up the SLP circuit
			MatchResult outputResult = find_matching_element(output);
			if (outputResult.newVariable)
			{
				vars.add(outputResult.element);
			}
			slpOutput = outputResult.element;
			MatchResult op1Results = find_matching_element(o1);
			if (op1Results.newVariable)
			{
				vars.add(op1Results.element);
			}
			slpO1 = op1Results.element;
			MatchResult op2Results = find_matching_element(o2);
			if (op2Results.newVariable)
			{
				vars.add(op2Results.element);
			}
			slpO2 = op2Results.element;

			// Build the line and add it to the list
			SLPLine slpLine = new SLPLine(slpOutput, slpO1, slpO2, slpOp);
			lines.add(slpLine);
		}

		// Skip "end"
		line = in.readLine();

		//////////////////////////////////////////////////////
		// Any further reading beyond end would go here
		// we'll put some dummy debug input data for testing purposes
		//////////////////////////////////////////////////////
		line = in.readLine();

		if (execute)
		{
			System.out.println("\n\n");
			while (line != null && line.length() > 0)
			{
				// Drop in the inputs so that the circuit can be executed
				Scanner dataIn = new Scanner(line);
				for (int i = 0; i < inputs.size(); i++)
				{
					int v = Integer.parseInt(dataIn.next());
					inputs.get(i).val = v;
				}

				// Display inputs and corresponding outputs
				SLP slp = new SLP(lines);
				ArrayList<Integer> inputVals = new ArrayList<Integer>();
				for (int i = 0; i < inputs.size(); i++)
				{
					inputVals.add(inputs.get(i).val);
				}
				System.out.println("-----");
				System.out.println("Input:  " + inputVals.toString());

				//// LOAD INPUTS AND BEGIN EXECUTION
				slp.load_inputs(inputs);
				slp.execute();
				//// END EXECUTION

				ArrayList<String> outputNames = new ArrayList<String>();
				for (int i = 0; i < outputs.size(); i++)
				{
					outputNames.add(outputs.get(i).name);
				}
				ArrayList<Integer> outputVals = slp.extract_outputs(outputNames);
				System.out.println("Output: " + outputVals.toString());
				System.out.println("-----\n");

				// See if there's more input to parse...
				line = in.readLine();
			}
		}
	}

	// Parse an SLP file, which is composed of an arbitrary number of lines
	// Format is as follows:
	// 
	// 13 gates
	// 5 inputs
	// a1 a2 ONE a4 a5
	// 1 outputs
	// a79
	// begin
	// a67 = a4  +  a1
	// a68 = a67  +  ONE
	// a69 = a5  +  a4
	// a70 = a2  +  a69
	// a71 = a70  +  a1
	// a72 = a69  x  a71
	// a73 = a68  +  a72
	// a74 = a5  +  a2
	// a75 = a73  x  a74
	// a76 = a5  +  a75
	// a77 = ONE  +  a71
	// a78 = a77  #  a75
	// a79 = a76  x  a78
	// end
	// 
	public static void main(String[] args)
	{
		if (args.length != 1)
		{
			System.err.println("usage: java -jar SLPParser slpFile");
			System.exit(-1);
		}

		try
		{
			SLPParser parser = new SLPParser();
			parser.parse(args[0], true);
		}
		catch (Exception e)
		{
			System.err.println(e.getMessage());
		}
	}
}