package YARF;


import java.io.Serializable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;


/**
 * This class handles the parallelization of many Gibbs chains over many CPU cores
 * to create one BART regression model. It also handles all operations on the completed model.
 * 
 * @author Adam Kapelner
 */
public class YARF extends Classifier implements Serializable {
	private static final long serialVersionUID = -6984205353140981153L;
	
	/** the number of CPU cores to build many different trees in a YARF model */
	protected int num_cores = 1; //default
	/** the number of trees in this RF model on all Gibbs chains */
	protected int num_trees = 500; //default
	
	/** should we print select messages to the screen */
	protected boolean verbose = true;
	/** ? */
	protected boolean mem_cache_for_speed = true;

	private YARFTree[] yarf_trees;

	
	/** the default constructor sets the number of total iterations each Gibbs chain is charged with sampling */
	public YARF(){}
	
	//init
	//set num cores, seed, verbose, etc
	//add data
	//set data feature names
	//add "other" data
	//set "other" data feature names
	//set num trees
	//init trees
	//give bootstrap samples (indices)
	//load all custom functions
	//BUILD
	
	
	public void setNumCores(int num_cores){
		this.num_cores = num_cores;
	}
	
	public void setNumTrees(int num_trees){
		this.num_trees = num_trees;
	}
	
	public void initTrees(){
		yarf_trees = new YARFTree[num_trees];
		for (int t = 0; t < num_trees; t++){
			yarf_trees[t] = new YARFTree(this);
		}
	}
	
	public void setBootstrapIndices(int[] bootstrap_indices, int t){
		yarf_trees[t].setBoostrapIndices(bootstrap_indices);
	}

	/** This function builds the forest by building all the trees */
	public void Build() {
		//run a build on all threads
		long t0 = System.currentTimeMillis();
		if (verbose){
			System.out.println("building YARF " + (mem_cache_for_speed ? "with" : "without") + " mem-cache speedup...");
		}

		ExecutorService tree_grow_pool = Executors.newFixedThreadPool(num_cores);
		for (int t = 0; t < num_trees; t++){
			final int tf = t;
	    	tree_grow_pool.execute(new Runnable(){
				public void run() {
					yarf_trees[tf].Build();
				}
			});
		}
		tree_grow_pool.shutdown();
		try {	         
	         tree_grow_pool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS); //effectively infinity
	    } catch (InterruptedException ignored){}	
		
		if (verbose){
			System.out.println("done building YARF in " + ((System.currentTimeMillis() - t0) / 1000.0) + " sec \n");
		}
	}


	/**
	 * Return the number of times each of the attributes were used during the construction of the sum-of-trees
	 * by Gibbs sample.
	 * 
	 * @param type	Either "splits" or "trees" ("splits" means total number and "trees" means sum of binary values of whether or not it has appeared in the tree)
	 * @return		The counts for all Gibbs samples further indexed by the attribute 1, ..., p
	 */
//	public int[][] getCountsForAllAttribute(final String type) {
//		final int[][] variable_counts_all_gibbs = new int[num_gibbs_total_iterations - num_gibbs_burn_in][p];		
//		
//		for (int g = 0; g < num_gibbs_total_iterations - num_gibbs_burn_in; g++){
//			final bartMachineTreeNode[] trees = gibbs_samples_of_bart_trees_after_burn_in[g];
//			int[] variable_counts_one_gibbs = new int[p];
//			for (bartMachineTreeNode tree : trees){	
//				if (type.equals("splits")){
//					variable_counts_one_gibbs = Tools.add_arrays(variable_counts_one_gibbs, tree.attributeSplitCounts());
//				}
//				else if (type.equals("trees")){
//					variable_counts_one_gibbs = Tools.binary_add_arrays(variable_counts_one_gibbs, tree.attributeSplitCounts());
//				}				
//				
//			}
//			variable_counts_all_gibbs[g] = variable_counts_one_gibbs;
//		}		
//		return variable_counts_all_gibbs;
//	}
	
	/**
	 * Return the proportion of times each of the attributes were used (count over total number of splits) 
	 * during the construction of the sum-of-trees by Gibbs sample.
	 * 
	 * @param type	Either "splits" or "trees" ("splits" means total number and "trees" means sum of binary values of whether or not it has appeared in the tree)
	 * @return		The proportion of splits for all Gibbs samples further indexed by the attribute 1, ..., p
	 */
//	public double[] getAttributeProps(final String type) {
//		int[][] variable_counts_all_gibbs = getCountsForAllAttribute(type);
//		double[] attribute_counts = new double[p];
//		for (int g = 0; g < num_gibbs_total_iterations - num_gibbs_burn_in; g++){
//			attribute_counts = Tools.add_arrays(attribute_counts, variable_counts_all_gibbs[g]);
//		}
//		Tools.normalize_array(attribute_counts); //will turn it into proportions
//		return attribute_counts;
//	}
	
	/**
	 * For all Gibbs samples after burn in, calculate the set of interaction counts (consider a split on x_j 
	 * and a daughter node splits on x_k and that would be considered an "interaction")
	 * 
	 * @return	A matrix of size p x p where the row is top split and the column is a bottom split. It is recommended to triangularize the matrix after ignoring the order.
	 */
//	public int[][] getInteractionCounts(){
//		int[][] interaction_count_matrix = new int[p][p];
//		
//		for (int g = 0; g < gibbs_samples_of_bart_trees_after_burn_in.length; g++){
//			bartMachineTreeNode[] trees = gibbs_samples_of_bart_trees_after_burn_in[g];
//			
//			for (bartMachineTreeNode tree : trees){
//				//get the set of pairs of interactions
//				HashSet<UnorderedPair<Integer>> set_of_interaction_pairs = new HashSet<UnorderedPair<Integer>>(p * p);
//				//find all interactions
//				tree.findInteractions(set_of_interaction_pairs);
//				//now tabulate these interactions in our count matrix
//				for (UnorderedPair<Integer> pair : set_of_interaction_pairs){
//					interaction_count_matrix[pair.getFirst()][pair.getSecond()]++; 
//				}
//			}	
//		}
//		
//		return interaction_count_matrix;
//	}

	/** Flush all unnecessary data from the Gibbs chains to conserve RAM */
	protected void FlushData() {
		for (int t = 0; t < num_cores; t++){
			yarf_trees[t].FlushData();
		}
	}

	/**
	 * The default BART evaluation of a new observations is done via sample average of the 
	 * posterior predictions. Other functions can be used here such as median, mode, etc. 
	 * Default is to use one CPU core.
	 * 
	 * @param record				The observation to be evaluated / predicted
	 */
	public double Evaluate(double[] record) {	
		return 0;
	}	
	
	/**
	 * The default BART evaluation of a new observations is done via sample average of the 
	 * posterior predictions. Other functions can be used here such as median, mode, etc.
	 * 
	 * @param record				The observation to be evaluated / predicted
	 * @param num_cores_evaluate	The number of CPU cores to use during evaluation
	 */
	public double Evaluate(double[] record, int num_cores_evaluate) {		
		return 0;
	}		
	
	/**
	 * After burn in, find the depth (greatest generation of a terminal node) of each tree for each Gibbs sample
	 * 
	 * @param thread_num	which CPU core (which Gibbs chain) to return results for
	 * @return				for each Gibbs chain return a vector of depths for all <code>num_trees</code> chains
	 */
	public int[] getDepthsForTrees(int thread_num){
		return null;
	}	
	
	/**
	 * After burn in, return the number of total nodes (internal plus terminal) of each tree for each Gibbs sample
	 * 
	 * @param thread_num	which CPU core (which Gibbs chain) to return results for
	 * @return				for each Gibbs chain return a vector of number of nodes for all <code>num_trees</code> chains
	 */
	public int[][] getNumNodesAndLeavesForTrees(){
		return null;
	}
	

	public void setVerbose(boolean verbose){
		this.verbose = verbose;
	}

	public void setSeed(int seed){
		StatToolbox.setSeed(seed);
	}
	
	public void setMemCacheForSpeed(boolean mem_cache_for_speed){
		this.mem_cache_for_speed = mem_cache_for_speed;
	}
	
	/** Must be implemented, but does nothing */
	public void StopBuilding() {}

	@Override
	public Classifier clone() {
		// TODO Auto-generated method stub
		return null;
	}	
}
