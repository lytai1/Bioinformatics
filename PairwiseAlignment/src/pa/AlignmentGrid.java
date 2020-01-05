package pa;
/**
 * Alignment grid for pairwise alignment of two DNA sequencde
 *
 * @author Kathleen
 * @version 11/05/2018
 */
/**
 * @author Kathleen
 *
 */
/**
 * @author Kathleen
 *
 */
public class AlignmentGrid
{
    // instance variables
    private Square[][] grid;
    private String sequence;
    private String sequence2;

    //direction constant
    private static int W = 0;
    private static int N = 1;
    private static int NW = 2;

    
    /**
     * Constructor of Alignment Grid
     * @param x DNA sequence
     * @param y DNA sequence
     */
    public AlignmentGrid(String x, String y)
    {
        // Initialize instance variables
        grid = new Square[y.length()+1][x.length()+1];
        sequence = x;
        sequence2 = y;

        // fill array with new Square
        for (int i = 0; i < grid.length; i++)
        {
            for (int j = 0; j< grid[0].length; j++)
            {
                grid[i][j] = new Square();
            }
        }

        grid[0][0].setScore(0);
        //fill the first row
        for (int i = 1; i < grid[0].length;i++)
        {
            grid[0][i] = new Square();
            grid[0][i].setScore(grid[0][i-1].getScore() - 2);
            grid[0][i].setDirection(W);
        }

        //fill the first column
        for (int i = 1; i < grid.length;i++)
        {
            grid[i][0].setScore(grid[i-1][0].getScore() - 2);
            grid[i][0].setDirection(N);
        }

        //fill the middle squares
        for (int j = 1; j < grid.length; j++)
        {
            for (int k = 1; k < grid[0].length; k++)
            {
                int up;
                int side;
                int left;

                //get the score from nearby grid
                up = grid[j-1][k].getScore();
                side = grid[j-1][k-1].getScore();
                left = grid[j][k-1].getScore();

                //calculate score
                up = up-2;
                left = left-2;

                char chr1 = sequence.charAt(k-1);
                char chr2 = sequence2.charAt(j-1);

                if (chr1==chr2)
                {
                    side = side +1;
                }
                else
                {
                    side = side - 1;
                }

                // compare the three directions
                // get greatest score
                int[] scoreList = {left, up, side};
                int max=scoreList[0];
                int maxIndex=-1;
                for (int i=0;i<scoreList.length;i++) {
                    if(scoreList[i]>max)
                    { max=scoreList[i];
                        grid[j][k].setDirection(i);//set direction
                    }
                }
                grid[j][k].setScore(max);//set score
            }
        }
    }

    /**
     * get Scoring array
     * @return array containing score of each square
     */
    public int[][] getScoreArray()
    {
        int[][] score = new int[sequence2.length() +1][sequence.length()+1];

        for (int i = 0; i < score.length; i++)
        {
            for (int j = 0; j< score[0].length; j++)
            {
                score[i][j] = grid[i][j].getScore();
            }
        }
        return score;
    }

    /**
     * get direction array
     * @return array containing direction of each grid
     */
    public int[][] getDirectionArray()
    {
        int[][] direction = new int[sequence2.length() +1][sequence.length()+1];

        for (int i = 0; i < direction.length; i++)
        {
            for (int j = 0; j< direction[0].length; j++)
            {
                direction[i][j] = grid[i][j].getDirection();
            }
        }
        return direction;
    }

    /**
     * Print array in console
     * @param array
     */
    public void printArray(int[][] array)
    {
        for (int i = 0; i < array.length; i++)
        {
            System.out.print("[");
            for (int j = 0; j< array[0].length; j++)
            {
                System.out.print(array[i][j] + " ");
            }
            System.out.println("]");
        }
    }

    /**
     * trace back the paired sequence and print out the aligned sequence in console
     */
    public void trace()
    {
        String paired1 = "";
        String paired2 = "";

        int i = sequence2.length();
        int j = sequence.length();

        while (i!=0&&j!=0)
        {
            int direction = grid[i][j].getDirection();
            if(direction == N)
            {
                paired1 = paired1 + "-";
                paired2 = paired2 + sequence2.charAt(i-1);
                i--;
            }
            else if (direction == W)
            {
                paired2 = paired2 + "-";
                paired1 = paired1 + sequence.charAt(j-1);
                j--;
            }
            else
            {
                paired1 = paired1 + sequence.charAt(j-1);
                paired2 = paired2 + sequence2.charAt(i-1);
                i--;
                j--;
            }
        }

        //reverse sequences
        String reversed1 = "";
        String reversed2 = "";
        for(int p = paired1.length()-1 ; p>=0 ;p--)
        {
            reversed1 = reversed1 + paired1.charAt(p);
        
        }
        for(int p = paired2.length()-1 ; p>=0 ;p--)
        {
            reversed2 = reversed2 + paired2.charAt(p);
        
        }
        paired1 = reversed1;
        paired2 = reversed2;
        
        //print the paired sequences
        System.out.println(paired1);
        System.out.println(paired2);
    }
}
