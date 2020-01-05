package pa;

/**
 * main method for alignment of two sequences
 *
 * @author Kathleen
 * @version 11/05/2018
 */
public class Alignment
{
    public static void main(String[] args)
    {
        //AlignmentGrid a = new AlignmentGrid("trump","duck");
        AlignmentGrid a = new AlignmentGrid("acgggt", "gcgag");
        int[][] score = a.getScoreArray();
        a.printArray(score);
        int[][] direction = a.getDirectionArray();
        a.printArray(direction);
        a.trace();
    
    }
}

