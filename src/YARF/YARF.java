package YARF;


import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TIntHashSet;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.script.Compilable;
import javax.script.CompiledScript;
import javax.script.Invocable;
import javax.script.ScriptContext;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;

import org.apache.commons.math3.stat.StatUtils;


/**
 * Builds a YARF model in parallel
 * 
 * @author Adam Kapelner
 */
public class YARF extends Classifier implements Serializable {
	private static final long serialVersionUID = -6984205353140981153L;
	
	/** the number of CPU cores to build many different trees in a YARF model */
	protected int num_cores;
	/** the number of trees in this YARF model */
	protected int num_trees;
	/** is this a regression problem? (or a classification problem) */
	protected boolean is_a_regression;
	
	/** an array of the raw training data by COLUMN i.e. consisting of xj = [x1j, ..., xnj] with the last entry being [y1, ..., yn] */ 
	private transient TIntObjectHashMap<double[]> X_by_col;

	/** other data which may be useful for custom functions */
	protected transient ArrayList<double[]> Xother;
	/** feature names in other data */
	protected String[] other_data_names;

	private YARFTree[] yarf_trees;
	private transient int[][] bootstrap_indices;
	
	protected int mtry;
	protected int nodesize;
	

	//convenient pre-computed data to have around
	protected transient TIntObjectHashMap<int[]> all_attribute_sorts;
	private transient TIntHashSet indicies_one_to_n;
	protected transient ArrayList<Integer> indicies_one_to_p;
	
	//everything that has to do with scripts
	private transient ScriptEngine nashorn_js_engine;
    private transient Compilable compilingEngine;
	private transient String shared_funs;
	protected transient Invocable mtry_fun;
	protected transient Invocable nodesize_fun;
	protected transient Invocable cost_calc_fun;
	protected transient Invocable node_assign_fun;
	
	
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

	
	
	/**
	 * Adds an observation / record to the "other" data array. The
	 * observation is converted to doubles and the entries that are 
	 * unrecognized are converted to {@link #MISSING_VALUE}'s.
	 * 
	 * @param x_i	The observation / record to be added as a String array.
	 */
	public void addOtherDataRow(String[] x_i){
		//initialize data matrix if it hasn't been initialized already
		if (Xother == null){
			Xother = new ArrayList<double[]>();
		}
		
		//now add the new record
		final double[] record = new double[x_i.length];
		for (int i = 0; i < x_i.length; i++){
			try {
				record[i] = Double.parseDouble(x_i[i]);
			}
			catch (NumberFormatException e){
				record[i] = MISSING_VALUE;
//				System.out.println("missing value at record #" + X_y.size() + " attribute #" + i);
			}
		}				
		Xother.add(record);		
	}
	
	public void setOtherDataNames(String[] other_data_names){
		this.other_data_names = other_data_names;
	}
	
	public void setNumCores(int num_cores){
		this.num_cores = num_cores;
	}
	
	public void setNumTrees(int num_trees){
		this.num_trees = num_trees;
	}
	
	public void setPredType(String pred_type){
		if (pred_type == "regression"){
			is_a_regression = true;
		}
	}
	
	public void setMTry(int mtry){
		this.mtry = mtry;
	}
	
	public void setNodesize(int nodesize){
		this.nodesize = nodesize;
	}
	
	//all JS functions stuff --- shared must be set FIRST!!!!
	public void setSharedFunctions(String shared_funs) throws ScriptException{
		this.shared_funs = shared_funs;
	}
	
	private Invocable stringToInvokableCompiledFunction(String fun) throws ScriptException{
        //lazy load for this stuff
		if (nashorn_js_engine == null){
			nashorn_js_engine = new ScriptEngineManager().getEngineByName("nashorn");
		}
		if (compilingEngine == null){
			compilingEngine = (Compilable) nashorn_js_engine;
		}
		
		String fun_and_shared_libraries = fun + "\n\n" + shared_funs;
		CompiledScript cscript = compilingEngine.compile(fun_and_shared_libraries); 
        cscript.eval(nashorn_js_engine.getBindings(ScriptContext.ENGINE_SCOPE));
        return (Invocable)cscript.getEngine();
	}
	
	public void setMTryFunction(String mtry_fun) throws ScriptException{
        this.mtry_fun = stringToInvokableCompiledFunction(mtry_fun);
	}
	
	public void setNodesizeFunction(String nodesize_fun) throws ScriptException{
		this.nodesize_fun = stringToInvokableCompiledFunction(nodesize_fun);
	}
	
	public void setCostCalcFunction(String cost_calc_fun) throws ScriptException{
		this.cost_calc_fun = stringToInvokableCompiledFunction(cost_calc_fun);
	}
	
	public void setNodeAssignFunction(String node_assign_fun) throws ScriptException{
		this.node_assign_fun = stringToInvokableCompiledFunction(node_assign_fun);
	}
	
	public boolean customFunctionMtry(){
		return mtry_fun != null;
	}
	
	public boolean customFunctionNodesize(){
		return nodesize_fun != null;
	}
	
	public boolean customFunctionCostCalc(){
		return cost_calc_fun != null;
	}
	
	public boolean customFunctionNodeAssign(){
		return node_assign_fun != null;
	}
	
	public void initTrees(){
		yarf_trees = new YARFTree[num_trees];
		for (int t = 0; t < num_trees; t++){
			yarf_trees[t] = new YARFTree(this);
			setBootstrapAndOutOfBagIndices(t);
		}
	}
	
	public void setBootstrapAndOutOfBagIndices(int t){
		//make a copy
		yarf_trees[t].setTrainingIndices(bootstrap_indices[t]);
		//now get oob indices - it begins as the full thing then we subtract 
		//out the bootstrap indices of the tree
		TIntHashSet oob_indices = new TIntHashSet(indicies_one_to_n);
		oob_indices.removeAll(bootstrap_indices[t]);
		yarf_trees[t].setOutOfBagIndices(oob_indices);
	}
	
//	public int[] getSortedIndicesAtAttribute(int j, int sub_indices){
//		synchronized(all_attribute_sorts){
//			int[] all_indices = all_attribute_sorts.get(j);
//			if (all_indices == null){ //lazy create for this attribute
//				double[] xj = new double[n];
//				for (int i = 0; i < n; i++){
//					xj[i] = X.get(i)[j];
//				}
//			}			
//		}
//
//	}

	/** This function builds the forest by building all the trees */
	public void Build() {

		all_attribute_sorts = new TIntObjectHashMap<int[]>(p);
		initTrees();
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
		return defaultNodeAssignmentAggregation(allNodeAssignments(record));
	}	
	
	private double defaultNodeAssignmentAggregation(double[] y_preds){
		if (is_a_regression){
			return StatUtils.mean(y_preds); //the sample average
		}
		double[] modes = StatUtils.mode(y_preds); //there could be multiple modes
		return modes[StatToolbox.randInt(modes.length)]; //return one at random in the spirit of "random forests"		
	}
	
	public double[] allNodeAssignments(double[] record){
		double[] y_preds = new double[num_trees];
		for (int t = 0; t < num_trees; t++){
			y_preds[t] = yarf_trees[t].Evaluate(record);	
		}
		return y_preds;		
	}
	
	public double[] allNodeAssignments(double[] record, int num_cores_evaluate){
		//speedup for the dumb user
		if (num_cores_evaluate == 1){
			return allNodeAssignments(record);
		}
		
		final double[] y_preds = new double[num_trees];
		ExecutorService tree_eval_pool = Executors.newFixedThreadPool(num_cores_evaluate);
		for (int t = 0; t < num_trees; t++){
			final int tf = t;
			tree_eval_pool.execute(new Runnable(){
				public void run() {
					y_preds[tf] = yarf_trees[tf].Evaluate(record);
				}
			});
		}
		tree_eval_pool.shutdown();
		try {	         
			tree_eval_pool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS); //effectively infinity
	    } catch (InterruptedException ignored){}
		
		return y_preds;	
	}
	
	/**
	 * The default BART evaluation of a new observations is done via sample average of the 
	 * posterior predictions. Other functions can be used here such as median, mode, etc.
	 * 
	 * @param record				The observation to be evaluated / predicted
	 * @param num_cores_evaluate	The number of CPU cores to use during evaluation
	 */
	public double Evaluate(double[] record, int num_cores_evaluate) {
		return defaultNodeAssignmentAggregation(allNodeAssignments(record, num_cores_evaluate));
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
	
	public void addBootstrapIndices(int[] indices_t, int tree){
		bootstrap_indices[tree] = indices_t;
	}
	
	public void finalizeTrainingData(){
		super.finalizeTrainingData();
		//initialize other data that requires data to be finalized
		X_by_col = new TIntObjectHashMap<double[]>(p);
		bootstrap_indices = new int[num_trees][n];
		indicies_one_to_n = new TIntHashSet();
		for (int i = 0; i < n; i++){
			indicies_one_to_n.add(i);
		}
		indicies_one_to_p = new ArrayList<Integer>();
		for (int j = 0; j < p; j++){
			indicies_one_to_p.add(j);
		}
	}
	
	
	/**
	 * Given a training data set indexed by row, this produces a training
	 * data set indexed by column
	 * 
	 * @param j		The feature to get
	 * @return		The nx1 vector of that feature
	 */
	protected double[] getXj(int j) {
		double[] x_dot_j = X_by_col.get(j);
		if (x_dot_j == null){ //gotta build it
			synchronized(X_by_col){ //don't wanna build it twice so sync it
				x_dot_j = new double[n];
				for (int i = 0; i < n; i++){
					x_dot_j[i] = X.get(i)[j];
				}
				X_by_col.put(j, x_dot_j);
			}	
		}
		return x_dot_j;
	 }
	
	protected int[] sortedIndices(int j, int[] sub_indices){
		int[] indices_sorted_j = all_attribute_sorts.get(j);
		if (indices_sorted_j == null){ //we need to build it
			synchronized(all_attribute_sorts){ //don't want to do this twice so sync it
				indices_sorted_j = getSortedIndices(j);
				all_attribute_sorts.put(j, indices_sorted_j);				
			}
		}
		if (sub_indices != null){ //that means we want some of them only
			int n_sub = sub_indices.length;
			int[] sorted_sub_indices = new int[n_sub];
			for (int i_s = 0; i_s < n_sub; i_s++){
				sorted_sub_indices[i_s] = indices_sorted_j[sub_indices[i_s]];
			}			
		}
		return indices_sorted_j;
	}
	
	
	private class SortPair implements Comparable<SortPair>{
	  public int ind;
	  public double value;

	  public SortPair(int ind, double value){
		  this.ind = ind;
		  this.value = value;
	  }

	  @Override 
	  public int compareTo(SortPair o){
	    return Double.compare(value, o.value);
	  }
	}
	
	//Java is about annoying as they come... why isn't this implemented????
	private int[] getSortedIndices(int j) {
		double[] xj = getXj(j);
		ArrayList<SortPair> temp = new ArrayList<SortPair>(n);
		for (int i = 0; i < n; i++){
			temp.add(new SortPair(i, xj[i]));
		}
		Collections.sort(temp);
		int[] indices = new int[n];
		for (int i = 0; i < n; i++){
			indices[i] = temp.get(i).ind;
		}		
		return indices;
	}

	public void StopBuilding() {
		for (int t = 0; t < num_trees; t++){
			yarf_trees[t].StopBuilding();
		}
	}
}
