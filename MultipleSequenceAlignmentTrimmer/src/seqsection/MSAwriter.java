package seqsection;

import java.io.IOException;
import java.io.PrintWriter;


/**
 * MSA writer
 * @author Kathleen
 *
 */
public class MSAwriter {
	private PrintWriter thePrintWriter;

	/**
	 * Constructor of MSA writer
	 * @param thePrintWriter
	 */
	public MSAwriter(PrintWriter thePrintWriter) {
		this.thePrintWriter = thePrintWriter;
	}
	
	/**
	 * write a collection of sequence in fasta format
	 * @param rec collection of sequence
	 * @throws IOException
	 */
	public void writeRecord(SequenceCollection rec) throws IOException {
			for(Sequence s: rec) {
				thePrintWriter.println(">" + s.getOrganism());
				thePrintWriter.println(s.getSequence());
			}
	
	}

}
