package seqsection;

/**
 * Sequence of an organism
 * @author Kathleen
 *
 */
public class Sequence {
	private String organism;
	private String seq;
	
	/**
	 * Constructor of Sequence
	 * @param organism	name of organism
	 * @param seq	nucleotide sequence
	 */
	public Sequence(String organism, String seq) {
		this.organism = organism;
		this.seq = "";
		for(int i = 0; i<seq.length(); i++) {
			if(seq.charAt(i) != ' ') {
				this.seq += seq.charAt(i);
			}
		}
	}
	
	/**
	 * Constructor of Sequence with trim from first to last nucleotide
	 * @param s	nucleotide sequence
	 * @param first	first nucleotide of the orginal sequence
	 * @param last	last nucleotide of the original sequence
	 */
	public Sequence(Sequence s, int first, int last) {
		this.organism = s.organism;
		this.seq = s.seq.substring(first-1, last);
	}
	
	/**
	 * extend sequence with extSeq
	 * @param extSeq sequence added to orginal sequence
	 */
	public void extendSeq(String extSeq) {
		for(int i = 0; i<extSeq.length(); i++) {
			if(extSeq.charAt(i) != ' ') {
				this.seq += extSeq.charAt(i);
			}
		}
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "Sequence [organism =" + organism + ",\n sequence =" + seq + "]";
	}
	
	/**
	 * get the name of the organism
	 * @return	name of the organism
	 */
	public String getOrganism() {
		return organism;
	}
	
	/**
	 * get the nucleotide sequence
	 * @return sequence
	 */
	public String getSequence() {
		return seq;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		// TODO Auto-generated method stub
		return organism.hashCode()+seq.hashCode();
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		Sequence that = (Sequence) obj;
		if(!that.seq.equals(this.seq)) {
			return false;
		}
		else if(!this.organism.equals(that.organism)) {
			return false;
		}
		else {
			return true;
		}
	}
	

}
