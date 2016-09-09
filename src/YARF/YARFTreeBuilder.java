package YARF;

import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;

public class YARFTreeBuilder {

	private YARFTree tree;
	private YARF yarf;

	public YARFTreeBuilder(YARFTree tree) {
		this.tree = tree;
		this.yarf = tree.yarf;
		//the root node has a overall cost which we need to compute
		computeNodeCost(tree.root);
		//now we start splitting recursively
		splitNode(tree.root);
	}

	private void splitNode(YARFNode node) {
		
		if (tree.stop){
			return; //user hit the brakes so ditch
		}
		//if this node is too small (or whatever other reason), ditch from splitting and make it a leaf
		if (makeNodeLeaf(node)){
			node.is_leaf = true;
			node.assignYHat();
			if (YARF.DEBUG){node.printNodeDebugInfo("");}
			return; //ditch... because we're done...
		}
		
		//greedy search... set up the horses
		double lowest_total_split_cost = node.cost; //if you can't beat this... don't bother!!
		int lowest_cost_split_attribute = Integer.MIN_VALUE; //bad flag!
		double lowest_cost_split_value = Double.NaN; //bad flag!
		boolean lowest_send_missing_data_right = false;
		
		//which features can we split on in this node?
		int[] features_to_split_on = selectAttributesToTry(node);
		
		
		//this will house the optimal split
		YARFNode lowest_left_node = null;
		YARFNode lowest_right_node = null;
		if (YARF.DEBUG){System.out.println("features_to_split_on: " + Tools.StringJoin(features_to_split_on));}
//		System.out.println("n_n: " + n_n);
		
		//two options about missingness (a) there is no missingness in this feature (b) there is
		
		//there are two missingness options - check them in random order
		double r = StatToolbox.rand();
		boolean[] trueFalseRandomOrder = {r > 0.5, r <= 0.5};
		
		
		for (int j : features_to_split_on){	
			//first get the indices for this note sorted on attribute j and get the missings as well
			TIntArrayList ordered_nonmissing_indices_j = new TIntArrayList();
			TIntHashSet missing_indices_j = new TIntHashSet();
			yarf.sortedIndices(j, node.indices, ordered_nonmissing_indices_j, missing_indices_j);
			if (YARF.DEBUG){System.out.println("node.indices: " + Tools.StringJoin(node.indices));}
			if (YARF.DEBUG){System.out.println("ordered_nonmissing_indices_j: " + Tools.StringJoin(ordered_nonmissing_indices_j));}
			if (YARF.DEBUG){System.out.println("missing_indices_j: " + Tools.StringJoin(missing_indices_j));}

			int num_max_split_points = ordered_nonmissing_indices_j.size();

			double[] xj = yarf.getXj(j);
			
			missing_search : for (boolean send_missing_data_right : trueFalseRandomOrder){
				//iterate over all the cut points!
				for (int i_cut = 0; i_cut < num_max_split_points; i_cut++){
					double split_value = xj[ordered_nonmissing_indices_j.get(i_cut)];
					if (i_cut < num_max_split_points - 1 && 
							split_value == xj[ordered_nonmissing_indices_j.get(i_cut + 1)]){ //this is not a split point since the next in the sequence is the same
						continue;
					}
					//set up zygotes
					YARFNode putative_left = new YARFNode(node);
					YARFNode putative_right = new YARFNode(node);
					putative_left.indices = (TIntArrayList)ordered_nonmissing_indices_j.subList(0, i_cut + 1);
					putative_right.indices = (TIntArrayList)ordered_nonmissing_indices_j.subList(i_cut + 1, num_max_split_points);		
					
					//handle the indices from missingness L/R now
					if (!missing_indices_j.isEmpty()){ 
						if (send_missing_data_right){
							putative_right.indices.addAll(missing_indices_j);
						}
						else {
							putative_left.indices.addAll(missing_indices_j);
						}
					}

					//we should ditch if these don't work out
					if (makeNodeLeaf(putative_left) || makeNodeLeaf(putative_right)){
//						System.out.println(" !! makeNodeLeaf(putative_left) " + makeNodeLeaf(putative_left)
//								+ " || makeNodeLeaf(putative_right) " + makeNodeLeaf(putative_right));
						continue;
					}
					
					//these are now viable splits, so we compute cost
					
					computeNodeCost(putative_left);
					computeNodeCost(putative_right);
					
					double total_split_cost = totalChildrenCost(putative_left, putative_right);
					//System.out.println("total_split_cost: " + total_split_cost);
					if (total_split_cost < lowest_total_split_cost){
						if (YARF.DEBUG){System.out.println("beat cost @ " + total_split_cost);}
						lowest_total_split_cost = total_split_cost;
						lowest_cost_split_attribute = j;
						lowest_cost_split_value = split_value;
						lowest_left_node = putative_left;
						lowest_right_node = putative_right;
						lowest_send_missing_data_right = send_missing_data_right;
					}
				}
				if (missing_indices_j.isEmpty()){break missing_search;} //no need to check the other because it will be the same 
				//(this enforces randomness of L/R missingness sending)
			}
		}
		
		

		//if we have been unsuccessful in finding an advantageous split, we should not split
		if (lowest_total_split_cost >= node.cost){
			if (YARF.DEBUG){System.out.println("greedy search unsuccessful... for node: " + node.stringLocation(true));}
			node.is_leaf = true;
			node.assignYHat();
			if (YARF.DEBUG){node.printNodeDebugInfo("");}
			return;
		}
		
		if (YARF.DEBUG){System.out.println("greedy search successful!! for node: " + node.stringLocation(true)
				+ "\n previous cost: " + node.cost 
				+ "\n new cost: " + lowest_total_split_cost
				+ "\n lowest_cost_split_attribute: " + lowest_cost_split_attribute
				+ "\n lowest_left_node: " + lowest_left_node
				+ "\n lowest_right_node: " + lowest_right_node
				+ "\n send_missing: " + (lowest_send_missing_data_right ? "R" : "L")
				);}
		
		//otherwise we use the optimal split
		//first indicate the splitting rule in this node converting the two lucky zygotes into a fetus
		node.split_attribute = lowest_cost_split_attribute;
//		System.out.println("  splitValue " + yarf.getXj(lowest_cost_split_attribute)[lowest_cost_split_index]);
		node.split_value = lowest_cost_split_value;
		node.send_missing_data_right = lowest_send_missing_data_right;
		
		//then officially give birth to two children from the two fetuses
		node.left = lowest_left_node;
		node.right = lowest_right_node;

		if (YARF.DEBUG){node.printNodeDebugInfo("");}
		//and now recurse and split on the new children just created
		splitNode(node.left);
		splitNode(node.right);
	}

	private double totalChildrenCost(YARFNode putative_left, YARFNode putative_right) {
		if (yarf.customFunctionBothChildrenCostCalc()){
			return yarf.runBothChildrenCost(putative_left, putative_right);
		}
		if (yarf.is_a_regression){ //sum of SSEs
			return putative_left.cost + putative_right.cost; 
		}
		else { //average entropy among the children
			int nL = putative_left.nodeSize();
			int nR = putative_right.nodeSize();
			return (nL * putative_left.cost + nR * putative_right.cost) / (nL + nR);			
		}
	}

	private void computeNodeCost(YARFNode node) {
		if (yarf.customFunctionSingleNodeCostCalc()){
			node.cost = yarf.runSingleNodeCost(node);
		}
		else {
			double[] ys = node.node_ys();
			if (yarf.is_a_regression){ //the cost is the SSE (across both the left and right)
				node.y_pred = StatToolbox.sample_average(ys); //set it here for convenience
				node.cost = StatToolbox.sample_sum_sq_err(ys, node.y_pred);
			}
			else { //it's a classification - the "cost" is the negative entropy which is the negative of the gain...
				node.cost = StatToolbox.natural_negative_entropy(ys);
			}
		}
		//System.out.println("computeNodeCost node " + node + " cost = " + node.cost + " pred = " + node.y_pred + " size = " + node.nodeSize());
	}

	private boolean makeNodeLeaf(YARFNode node) {
//		System.out.println("makeNodeLeaf nodesize = " + node.nodeSize() + " yarfnodesize = " + yarf.nodesize);
	
		if (yarf.customFunctionNodesize()){
			return yarf.runNodesizeLegal(node);
		} else if (node.nodeSize() < yarf.nodesize){
			return true;
		}
		return false;
	}

	private int[] selectAttributesToTry(YARFNode node) {
		if (yarf.customFunctionMtry()){
			return yarf.runMtry(node);
		}
//		System.out.println("(yarf.p / 3.0)" +  (yarf.p / 3.0));
		int mtry = Math.max(1, (int)Math.floor(yarf.is_a_regression ? (yarf.p / (double)3) : Math.sqrt(yarf.p))); //at least it's 1!!!
		return StatToolbox.pickNRandomElements(yarf.indices_one_to_p_min_1, mtry);
	}

}
