package org.husonlab.fmhdist.ncbi;

import com.google.common.graph.Graph;
import com.google.common.graph.Graphs;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class TaxonomyTree {
	private Graph<Taxon> tree;
	private Map<Integer, Taxon> taxa;

	public TaxonomyTree(Graph<Taxon> tree, Map<Integer, Taxon> taxa) {
		this.tree = tree;
		this.taxa = taxa;
	}

	public Graph<Taxon> getTree() {
		return Graphs.copyOf(this.tree);
	}

	public Map<Integer, Taxon> getTaxa() {
		return new HashMap<>(this.taxa);
	}

	public Taxon getRoot() {
		// This is the root node in the NCBI taxonomy
		return this.taxa.get(1);
	}

	private String getSubtreeNewick(Taxon subtreeRoot) {
		if (this.tree.successors(subtreeRoot).isEmpty()) {
			return String.valueOf(subtreeRoot.getTaxonId());
		}
		StringBuilder result = new StringBuilder();
		result.append("(");
		List<String> subtrees = new ArrayList<>();
		for (Taxon child : this.tree.successors(subtreeRoot)) {
			subtrees.add(this.getSubtreeNewick(child));
		}
		result.append(String.join(",", subtrees));
		result.append(")" + subtreeRoot.getTaxonId());
		return result.toString();
	}

	public String getNewick() {
		return this.getSubtreeNewick(this.getRoot());
	}
}
