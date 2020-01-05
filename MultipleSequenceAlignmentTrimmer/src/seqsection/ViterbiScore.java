package seqsection;

/**
 * Viterbi score of an organism
 * @author Kathleen
 *
 */
public class ViterbiScore {
	private String organism;
	private String score;
	
	/**
	 * Constructor of ViterbiScore
	 * @param organism	name of organism
	 * @param score		viterbi score of organism
	 */
	public ViterbiScore(String organism, String score) {
		this.organism = organism;
		this.score = score;
	}
	/**
	 * get name of organism
	 * @return name
	 */
	public String getOrganism() {
		return organism;
	}
	/**
	 * get score
	 * @return viterbi score
	 */
	public String getScore() {
		return score;
	}
	
	
	
}
