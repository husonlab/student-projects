package org.husonlab.fmhdist.cmd;

import javafx.geometry.Point2D;
import jloda.graph.Edge;
import jloda.graph.Node;
import jloda.graph.NodeArray;
import jloda.phylo.PhyloSplitsGraph;
import jloda.util.FileLineIterator;
import jloda.util.parse.NexusStreamParser;
import jloda.util.progress.ProgressSilent;
import org.jfree.svg.SVGGraphics2D;
import org.jfree.svg.SVGHints;
import splitstree6.algorithms.distances.distances2splits.NeighborNet;
import splitstree6.data.DistancesBlock;
import splitstree6.data.SplitsBlock;
import splitstree6.data.TaxaBlock;
import splitstree6.io.readers.NexusImporter;
import splitstree6.layout.splits.algorithms.PhylogeneticOutline;

import java.awt.*;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

/**
 * Class to create an SVG of a phylogenetic outline.
 */
public class OutlineVisualizer {
	/**
	 * Calculates the phylogenetic outline for the given distances and stores it
	 * as an SVG.
	 *
	 * @param input   Path to the input distances in Nexus format
	 * @param output  Path to the output SVG
	 * @param width   width of the generated image
	 * @param height  height of the generated image
	 * @param scale   scaling factor, i.e. how large is the outline
	 * @param xOffset offset of the outline in x direction
	 * @param yOffset offset of the outline in y direction
	 * @param labels  optional path to a TSV file specifying the display names
	 *                for taxa
	 */
	public void run(String input, String output, int width, int height, int scale, int xOffset, int yOffset, String labels) {
		Logger logger = Logger.getLogger(OutlineVisualizer.class.getName());
		try {
			logger.info("Reading distances...");
			FileReader reader = new FileReader(input);
			NexusStreamParser parser = new NexusStreamParser(reader);
			TaxaBlock taxa = new TaxaBlock();
			DistancesBlock distances = new DistancesBlock();
			NexusImporter.parse(parser, taxa, distances);

			logger.info("Calculating splits...");
			SplitsBlock splits = new SplitsBlock();
			NeighborNet nn = new NeighborNet();
			nn.compute(new ProgressSilent(), taxa, distances, splits);

			logger.info("Calculating outline...");
			PhyloSplitsGraph graph = new PhyloSplitsGraph();
			NodeArray<Point2D> nodes = new NodeArray<>(graph);
			BitSet usedSplits = new BitSet();
			ArrayList<ArrayList<Node>> loops = new ArrayList<>();
			PhylogeneticOutline.apply(new ProgressSilent(), true, taxa, splits, graph, nodes, usedSplits, loops, 0, 0);

			Map<String, String> labelMap = null;
			if (labels == null || labels.equals("")) {
				logger.info("No labels provided, using taxa ids...");
			} else {
				logger.info("Reading label file...");
				labelMap = new HashMap<>();
				try (FileLineIterator it = new FileLineIterator(labels)) {
					while (it.hasNext()) {
						String line = it.next();
						String[] elements = line.split("\t");
						String value = null;
						if (elements.length > 1) {
							value = elements[1];
						}
						labelMap.put(elements[0], value);
					}
				}
			}

			logger.info("Creating output...");
			SVGGraphics2D graphics = new SVGGraphics2D(width, height);

			graphics.setPaint(Color.black);
			graphics.setRenderingHint(SVGHints.KEY_BEGIN_GROUP, "edges");
			for (Edge e : graph.edges()) {
				Node s = e.getSource();
				Node t = e.getTarget();

				Point2D p1 = nodes.get(s);
				Point2D p2 = nodes.get(t);

				int x1 = (int) Math.floor(p1.getX() * scale) + xOffset;
				int y1 = (int) Math.floor(p1.getY() * scale) + yOffset;
				int x2 = (int) Math.floor(p2.getX() * scale) + xOffset;
				int y2 = (int) Math.floor(p2.getY() * scale) + yOffset;

				graphics.drawLine(x1, y1, x2, y2);
			}
			graphics.setRenderingHint(SVGHints.KEY_END_GROUP, "edges");
			graphics.setRenderingHint(SVGHints.KEY_BEGIN_GROUP, "leaves");
			for (Node n : graph.leaves()) {
				String label = n.getLabel();
				if (labelMap != null) {
					if (labelMap.containsKey(label)) {
						label = labelMap.get(label);
						// Listing only the key in the map without the label
						// indicates that the label should be skipped
						if (label == null) {
							continue;
						}
					}
				}
				Point2D p = nodes.get(n);
				int x = (int) Math.floor(p.getX() * scale) + xOffset;
				int y = (int) Math.floor(p.getY() * scale) + yOffset;
				graphics.drawString(label, x, y);
			}
			graphics.setRenderingHint(SVGHints.KEY_END_GROUP, "leaves");

			File imageFile = new File(output);
			if (imageFile.exists()) {
				imageFile.delete();
			}
			FileWriter writer = new FileWriter(imageFile);
			writer.write(graphics.getSVGDocument());
			writer.close();
		} catch (Exception e) {
			System.out.println("well, f****");
			e.printStackTrace();
		}
	}
}
