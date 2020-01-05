package seqsection;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * ViterbiScoreReader reads viterbi score from a file
 * @author Kathleen
 *
 */
public class ViterbiScoreReader {
	private BufferedReader theBufferedReader;

	/**
	 * Constructor of ViterbiScoreReader
	 * @param theBufferedReader
	 */
	public ViterbiScoreReader(BufferedReader theBufferedReader) {
		this.theBufferedReader = theBufferedReader;
	}
	
	/**
	 * read viterbiScore from a file
	 * @return ArrayList carrying ViterbiScore
	 * @throws IOException
	 */
	public ArrayList<ViterbiScore> readViterbiScore()throws IOException{
		ArrayList<ViterbiScore> arraylist = new ArrayList<>();
		boolean done = false;
		
		while(!done) {
			String firstLine = theBufferedReader.readLine();
			if(firstLine!=null) {
			int start= firstLine.indexOf(">");
			String organism = firstLine.substring(start+1);
			theBufferedReader.readLine();
			String thirdLine = theBufferedReader.readLine();
			start = thirdLine.indexOf("= ");
			String score = thirdLine.substring(start+2);
			ViterbiScore v = new ViterbiScore(organism, score);
			arraylist.add(v);
			}
			else {
				done = true;
			}
		}
		return arraylist;
		
	}
	
	/**
	 * main method to read viterbi scores from file
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("Starting");
		try {
			String[] fileName = {"training", "generate", "random", "HBA", "HBB", "HBD"};
			for(String name: fileName) {
			File result = new File("result/HBBHMM_"+name+".txt");
			if (!result.exists()) {
				System.out.println("Can't find input file " + result.getAbsolutePath());
				System.exit(1);
			}
			
			//Build chain of readers.
			FileReader fr = new FileReader(result);
			BufferedReader br = new BufferedReader(fr);
			ViterbiScoreReader vsr = new ViterbiScoreReader(br);
					
			ArrayList<ViterbiScore> resultArrayList= vsr.readViterbiScore();
			
			br.close();
			fr.close();
			
			
			
			for(ViterbiScore vs: resultArrayList) {
				System.out.println(vs.getOrganism());
			}
			System.out.println();
			for(ViterbiScore vs: resultArrayList) {
				System.out.println(vs.getScore());
				
			}
			}
		} catch (IOException x) {
			System.out.println(x.getMessage());
		} 
		System.out.println("Done");
	}

}
