package seqsection;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;




/**
 * Multiple Sequence Alignment Reader
 * @author Kathleen
 *
 */
public class MSAreader {
	private BufferedReader theBufferedReader;

	/**
	 * Constructor of MSA reader
	 * @param theBufferedReader
	 */
	public MSAreader(BufferedReader theBufferedReader) {
		this.theBufferedReader = theBufferedReader;
	}

	/**
	 * read MSA alignment record
	 * @return a collection of sequences
	 * @throws IOException
	 * @throws RecordFormatException
	 */
	public SequenceCollection readRecord() throws IOException, RecordFormatException {

		theBufferedReader.readLine();
		theBufferedReader.readLine();
		theBufferedReader.readLine();

		SequenceCollection col = new SequenceCollection();

		boolean section1Done = false;
		boolean sectionDone = false;

	
			// section one read
			while (!section1Done) {
				String line = theBufferedReader.readLine();
				int space = line.indexOf(" ");
				if (space != 0) {

					String organism = line.substring(0, space);

					while (line.charAt(space) == ' ') {
						space++;
					}
					int first = space;
					int last = line.indexOf("	", first);
					String seq = line.substring(first, last);

					col.add(new Sequence(organism, seq));
				} else {
					theBufferedReader.readLine();
					section1Done = true;
				}
			}
			// other section read
			int seqNo = 0;
			while (section1Done && !sectionDone ) {
				String line = theBufferedReader.readLine();
				int space = line.indexOf(" ");
				if (space != 0) {

					String organism = line.substring(0, space);

					while (line.charAt(space) == ' ') {
						space++;
					}
					int first = space;
					int last = line.indexOf("	", first);
					String seq = line.substring(first, last);
					Sequence s = col.get(seqNo);
					if (!organism.equals(s.getOrganism())) {
						throw new RecordFormatException("organism name does not match");
					}
					s.extendSeq(seq);
					seqNo++;
					if(seqNo>=col.size()) {
						seqNo = 0;
					}
				} else {
					if (theBufferedReader.readLine() == null) {
						sectionDone = true;
					}
				}
			}
		

		return col;
	}
	
	/**
	 * Main method to read MSA file for trimming and write in fasta files
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("Starting");
		try {
			File fastq = new File("data/HBB_training_MSA.txt");
			if (!fastq.exists()) {
				System.out.println("Can't find input file " + fastq.getAbsolutePath());
				System.exit(1);
			}
			File fasta = new File("data/HBB_training_MSA.fasta");
			
			//Build chain of readers.
			FileReader fr = new FileReader(fastq);
			BufferedReader br = new BufferedReader(fr);
			MSAreader msar = new MSAreader(br);
					
			// Build chain of writers.
			FileWriter fw = new FileWriter(fasta);
			PrintWriter pw = new PrintWriter(fw);
			MSAwriter faw = new MSAwriter(pw);
			
			SequenceCollection col = msar.readRecord();
			
			
			//trim sequence
			SequenceCollection trimed = new SequenceCollection(col, 20,144);
			faw.writeRecord(trimed);
			
			for (Sequence s: col) {
				System.out.println(s.getOrganism());
			}
			
			pw.close();
			fw.close();
			br.close();
			fr.close();
			
		} catch (IOException x) {
			System.out.println(x.getMessage());
		} catch (RecordFormatException y) {
			System.out.println(y.getMessage());
		}
		System.out.println("Done");
	}

}
