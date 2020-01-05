package seqsection;

import java.util.ArrayList;

/**
 * Collection of sequences
 * @author Kathleen
 *
 */
public class SequenceCollection extends ArrayList<Sequence>{

	public SequenceCollection() {
		super();
		// TODO Auto-generated constructor stub
	}

	public SequenceCollection(SequenceCollection seqCol, int first, int last) {
		for(Sequence s: seqCol) {
			Sequence cutSeq = new Sequence(s, first, last);
			this.add(cutSeq);
		}
	}

//	@Override
//	public String toString() {
//		return "SequenceCollection []";
//	}
	

}
