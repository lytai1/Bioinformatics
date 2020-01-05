package pa;


/**
 * Square in AlignmentGrid
 *
 * @author Kathleen
 * @version 11/05/2018
 */
public class Square
{
    // instance variables - replace the example below with your own
    private int score;
    private int direction;

    /**
     * Constructor for objects of class square
     */
    public Square()
    {
        // Initialize instance variables
        this.direction = direction;
        this.score = score;
    }

    
    /**
     * set score of the square
     * @param score alignment score
     */
    public void setScore(int score)
    {
        // put your code here
        this.score = score;

    }
    
    /**
     * get the alignment score of the square
     * @return score
     */
    public int getScore()
    {
        // put your code here
        return score;

    }
    
    /**
     * set direction where the score calculated from
     * @param direction
     */
    public void setDirection(int direction)
    {
        // put your code here
        this.direction = direction;

    }
    
    /**
     * get direction where the score calculated from
     * @return direction
     */
    public int getDirection()
    {
        // put your code here
        return direction;

    }
}
