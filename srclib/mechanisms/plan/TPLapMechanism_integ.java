package priv.plan.mechanisms;

import org.jblas.DoubleMatrix;

// use method get_cell_probability to get the prob of a single cell and then truncate 
// on the border when convergence is detected. (old and non efficient)

// Truncated planar Laplacian mechanism
// The mechanism matrix has the form 
// (double) mechanism[i_x][i_y][j_x][j_y] where the first 
// coordinates represent the row (mechanism input) and the other two 
// coordinates represent the column (mechanism output). 


public class TPLapMechanism_integ extends Mechanism {
	
      private double eps;                    // used to construct Laplace matrix 
      private DoubleMatrix mechanism_matrix;  
      // stores the entries of the mechanism since every entry takes long time to be computed
    
  	  final int infinity = 200; // maximum index used in numerical summations.
   
    // this constructor creates the TLap mechanism
    // specified by the size of the grid, the cell_side (in km), the type, and privacy parameter
    public TPLapMechanism_integ(int p_hsize, int p_vsize, double p_cell_side, double p_eps) {
    	super(p_hsize,p_vsize,p_cell_side);
    	mechanism_matrix = new DoubleMatrix(p_hsize*p_vsize,p_hsize*p_vsize);	
        eps = p_eps;
    	fill_matrix_lap();
    	
    	String info = "// Truncated Laplacian mechanism  \n"
                + "// epsilon =" + eps +" for 1km \n"
                + "// Size: hsize = "+ hsize +", vsize = "+ vsize; 
    	this.add_info(info);
    }
    
    // fill the mechanism matrix elements with the Laplace values
    private void fill_matrix_lap() {
    	mechanism_matrix.fill(0);
        for (int in_x=0; in_x<hsize;in_x++)
    		for (int in_y=0; in_y<vsize;in_y++) 
               	for (int out_x=0;out_x<hsize;out_x++)
               		for (int out_y=0;out_y<vsize;out_y++) 
               			mechanism_matrix.put(in_y*hsize+in_x,out_y*hsize+out_x, get_Lap_prob(in_x,in_y,out_x,out_y));
    }
    
 	// returns the entry of the Lap mechanism 
 	// for input x,y and output x,y
 	private double get_Lap_prob(int p_in_x,int p_in_y,int p_out_x,int p_out_y) {
 		
 		
 		// a bit of optimization: use the symmetry of the mechanism with 
 		// the symmetry of the grid
 		// P(in_x,in_y,out_x,out_y) = P(hsize-1-in_x,vsize-1-in_y,hsize-1-out_x,vsize-1-out_y)
 		double ret = mechanism_matrix.get( (vsize-1-p_in_y)*hsize+(hsize-1-(hsize-1-p_in_x)),
 				(vsize-1-p_out_y)*hsize+(hsize-1-(hsize-1-p_out_x)) );
 		if (ret>0) return ret;
 		
 		int out_x_start= p_out_x; // starting x for summation.
 		int out_x_end =  p_out_x; // end x for summation.
 		int out_y_start= p_out_y; // starting y for summation.
 		int out_y_end =  p_out_y; // ending y for summation.
 		
 		int xstep=1, ystep=1;  // specifies whether the summations goes forward or backward.
 		
 		// adjust the boundaries of the sum if out is on the borders
 		if (p_out_x==0) {
 			out_x_end = out_x_end-infinity; 
 			xstep=-1;
 			}; 
 		if (p_out_x==hsize-1) {
 			out_x_end = out_x_end+infinity; 
 			xstep =1; 
 		}
 		if (p_out_y==0) {
 			out_y_end = out_y_end-infinity; 
 			ystep = -1;
 		}
 		if (p_out_y==vsize-1) {
 			out_y_end =out_y_end+infinity; 
 			ystep = 1;
 		}

        // now make the summation to compute the required probability
 		double sum=0;

 		int out_x_curr, out_y_curr;
 		double xprob; 
 		double yprob;
 		boolean x_min_precision_reached = false; 
 		boolean y_min_precision_reached = false; 
 		
 		out_x_curr=out_x_start;  
 		while ( ((xstep>0 && out_x_curr <= out_x_end) || (xstep< 0 && out_x_curr >= out_x_end)) 
 				&& !x_min_precision_reached) {
        	xprob=0;
        	out_y_curr=out_y_start;
        	y_min_precision_reached = false;
        	while (((ystep>0 && out_y_curr <= out_y_end) || (ystep< 0 && out_y_curr >= out_y_end))
        			&& !y_min_precision_reached) {
        		yprob = get_cell_prob(p_in_x,p_in_y,out_x_curr,out_y_curr);
    		    xprob = xprob + yprob;
        		out_y_curr=out_y_curr + ystep;
            	if (yprob<Math.pow(0.1,Settings.prop_decimal_places))
            			y_min_precision_reached = true;
        	} 
        	out_x_curr = out_x_curr+xstep;
        	sum = sum + xprob;
        	if (xprob<Math.pow(0.1,Settings.prop_decimal_places))
        			x_min_precision_reached = true;
        } 
 		
 		return sum;
 		
 	}

    
    
    // returns the total probability (by integration) of one cell having the center (p_out)
    // given that the input is (p_in). 
 	private double get_cell_prob(int p_in_x,int p_in_y,int p_out_x,int p_out_y) {
 
 		
 		double epsd=eps*cell_side;  // compute the priv for adjacent cells. 
 		int divisions = 10;  // number of divisions for integration

 		double delta= (double)1/divisions; //(unit distance) delta for integration
 		double out_x_start= (double)p_out_x-0.5; // starting x for integration.
 		double out_x_end =  (double)p_out_x+0.5; // end x for integration.
 		double out_y_start= (double)p_out_y -0.5; // starting y for integration.
 		double out_y_end =  (double)p_out_y +0.5; // ending y for integration.
 		
        double dist; // distance between input and output positions
 		double sum=0;
 		double prob;

 		double out_x_curr=out_x_start;  
        while (out_x_curr<out_x_end-0.5*delta) {
            double out_y_curr=out_y_start;  
        	while (out_y_curr<out_y_end-0.5*delta) {
    		    dist=Math.sqrt(Math.pow(out_x_curr-(double)p_in_x, 2)+
    		    		       Math.pow(out_y_curr-(double)p_in_y, 2));
        		prob=(double)(epsd*epsd)/(2*Math.PI)*Math.exp(-1*epsd*dist);
    		    sum = sum + (delta*delta)*prob;
        		out_y_curr=out_y_curr + delta;
        	}
        	out_x_curr=out_x_curr+delta;
        }
 		return sum;
 	}
 	
    // returns the mechanism matrix: override the super class method to improve performance 
    public DoubleMatrix get_mech_matrix() {
 		return mechanism_matrix;
    }

    // returns a mechanism matrix row probabilities: override the super class method to improve performance
    // used For EM. p_in is the out location in linear order, i.e. out_y*hsize +out_x. 
    public double[] get_mech_row_probabilities(int p_in) {
 	   return mechanism_matrix.getRow(p_in).toArray();
    }

    // returns the mechanism matrix column: override the super class method to improve performance
    // p_out is the out location in linear order, i.e. out_y*hsize +out_x. 
    public DoubleMatrix get_mech_column(int p_out) {
    	return mechanism_matrix.getColumn(p_out);
    }
 	
 	public double getProb(int p_in_x, int p_in_y, int p_out_x, int p_out_y) {
 		return mechanism_matrix.get(p_in_y*hsize+p_in_x,p_out_y*hsize+p_out_x);
 	}

}
