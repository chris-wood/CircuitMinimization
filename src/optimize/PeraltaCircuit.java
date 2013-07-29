import java.util.ArrayList;

public class PeraltaCircuit
{
	public ArrayList<Matrix> forms;
	public int numXorSignals;
	public int numAndSignals;
	public ArrayList<Integer> origXorSignals;
	public ArrayList<Integer> origAndSignals;
	public double xorProb;
	public double andProb;
	public int andCap;
	public long rounds;
	public long tries;

	public PeraltaCircuit(ArrayList<Matrix> fs, int nxs, int nas, 
		double xp, double ap, int ac, long r, long t)
	{
		this.forms = fs;
		this.numXorSignals = nxs;
		this.numAndSignals = nas;
		this.xorProb = xp;
		this.andProb = ap;
		this.andCap = ac;
		this.rounds = r;
		this.tries = t;
	}
}

