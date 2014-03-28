package caw.circuits.optimize.linear;

public class SLPElement
{
	public String name;
	public int val;

	public SLPElement(String n, int v) {this.name = new String(n); this.val = v; }

	@Override
	public boolean equals(Object o)
	{
		if (o instanceof SLPElement)
		{
			SLPElement other = (SLPElement)o;
			if (other.name.equals(this.name) && this.val == other.val)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		else
		{
			return false;
		}
	}

	@Override
	public int hashCode()
	{
		return name.hashCode() + val;
	}
}