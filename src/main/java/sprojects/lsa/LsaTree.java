/*
 *  LsaTree.java Copyright (C) 2024 Algorithms in Bioinformatics
 *
 *  (Some files contain contributions from other authors, who are then mentioned separately.)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package sprojects.lsa;

import jloda.fx.util.ArgsOptions;
import jloda.graph.Node;
import jloda.graph.NodeArray;
import jloda.phylo.PhyloTree;
import jloda.util.FileUtils;
import jloda.util.UsageException;
import splitstree6.algorithms.trees.trees2trees.LSATree;

import java.io.BufferedReader;
import java.io.IOException;

/**
 * reads rooted networks as input and writes their LSA trees as output, or checks some property
 * May 2024
 */
public class LsaTree {
	/**
	 * run the program
	 *
	 * @param args commandline arguments
	 * @throws UsageException incorrect usage
	 * @throws IOException    problem reading or writing file
	 */
	public static void main(String[] args) throws UsageException, IOException {
		var options = new ArgsOptions(args, LsaTree.class, "Compute the LSA tree and check property");
		var infile = options.getOption("-i", "input", "Input file (stdin, .gz ok)", "stdin");
		var outfile = options.getOption("-o", "output", "Output file (stdout, .gz ok", "stdout");
		var property = options.getOption("-p", "property", "Name of property to check", "");
		options.done();

		try (var r = new BufferedReader(FileUtils.getReaderPossiblyZIPorGZIP(infile));
			 var w = FileUtils.getOutputWriterPossiblyZIPorGZIP(outfile)) {
			while (r.ready()) {
				var line = r.readLine();
				if (line.startsWith("(")) {
					var network = new PhyloTree();
					network.parseBracketNotation(line, true);
					var tree = new PhyloTree();
					try (NodeArray<Node> srcTargetMap = network.newNodeArray()) { // srcTargetMap maps network nodes to LSA nodes
						LSATree.computeLSA(network, tree, srcTargetMap);

						if (property.isBlank()) {
							w.write(tree.toBracketString(true) + ";\n");
						} else {
							System.err.println("Check property '" + property + "'");
							w.write(property + "\n");
						}
					}
				}
			}
		}
	}
}
