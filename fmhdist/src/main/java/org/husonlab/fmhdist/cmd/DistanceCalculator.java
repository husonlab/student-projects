package org.husonlab.fmhdist.cmd;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.husonlab.fmhdist.sketch.Distance;
import org.husonlab.fmhdist.sketch.FracMinHashSketch;

import jloda.thirdparty.HexUtils;
import jloda.util.FileLineIterator;
import splitstree6.data.DistancesBlock;
import splitstree6.data.TaxaBlock;
import splitstree6.io.writers.distances.NexusWriter;

/**
 * Class to calculate the evolutionary distances based on input sketches.
 */
public class DistanceCalculator {
	private boolean LOG_EMPTY_INTERSECTIONS = false;

	private FracMinHashSketch lineToSketch(String line) {
		try {
			String[] comp = line.replaceAll("\\s+", "").split(",");
			byte[] content = Files.readAllBytes(Paths.get(comp[0]));
			FracMinHashSketch result = FracMinHashSketch.parse(HexUtils.decodeHexString(new String(content)));
			if (comp.length > 1) {
				result.setName(comp[1]);
			} else {
				result.setName(comp[0]);
			}
			return result;
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Calculates three sets evolutionary distances: Mash Distance, FracMinHash
	 * containment distance and FracMinHash distance.
	 *
	 * @param input  A CSV of sketches that should be included in distance
	 *               calculations, each line must specify the path to the sketch file
	 * @param output Path to the output FracMinHash distance in Nexus format.
	 *               The other two distances are stored using ".mash" and ".containment"
	 *               suffixes.
	 */
	public void run(
			String input,
			String output
	) {
		Logger logger = Logger.getLogger(DistanceCalculator.class.getName());
		try {
			logger.info("Reading queries list...");
			FileLineIterator it = new FileLineIterator(input);
			List<FracMinHashSketch> sketches = it
					.stream()
					.map(this::lineToSketch)
					.collect(Collectors.toList());
			it.close();

			DistancesBlock distances_jaccard = new DistancesBlock();
			distances_jaccard.setNtax(sketches.size());

			DistancesBlock distances_mash = new DistancesBlock();
			distances_mash.setNtax(sketches.size());

			DistancesBlock distances_containment = new DistancesBlock();
			distances_containment.setNtax(sketches.size());

			TaxaBlock taxa = new TaxaBlock();
			int sParameter = 0;
			int kParameter = 0;
			int randomSeed = 0;
			long hashedMagicNumber = 0;

			List<String> emptyIntersectionLog = null;

			if (LOG_EMPTY_INTERSECTIONS) {
				emptyIntersectionLog = new ArrayList<>();
			}

			for (int i = 0; i < sketches.size(); i++) {
				if (sParameter != 0 || kParameter != 0 || randomSeed != 0 || hashedMagicNumber != 0) {
					if (sParameter != sketches.get(i).getSParam() || kParameter != sketches.get(i).getKSize() || randomSeed != sketches.get(i).getSeed() || hashedMagicNumber != sketches.get(i).getHashedMagicNumber()) {
						logger.severe("sketches have incompatible sketching parameters");
						return;
					}
				} else {
					sParameter = sketches.get(i).getSParam();
					kParameter = sketches.get(i).getKSize();
					randomSeed = sketches.get(i).getSeed();
					hashedMagicNumber = sketches.get(i).getHashedMagicNumber();
				}

				taxa.addTaxonByName(sketches.get(i).getName());
				for (int j = i; j < sketches.size(); j++) {
					double jaccard = Distance.calculateJaccardIndex(
							sketches.get(i).getValues(),
							sketches.get(j).getValues(),
							sParameter
					);
					// Containment is not symmetrical
					double containment_i = Distance.calculateContainmentIndex(
							sketches.get(i).getValues(),
							sketches.get(j).getValues(),
							sParameter
					);
					double containment_j = Distance.calculateContainmentIndex(
							sketches.get(j).getValues(),
							sketches.get(i).getValues(),
							sParameter
					);

					if (jaccard == 0 && emptyIntersectionLog != null) {
						emptyIntersectionLog.add(String.format("%s vs %s is empty\n", sketches.get(i).getName(), sketches.get(j).getName()));
					}

					distances_jaccard.setBoth(i + 1, j + 1, Distance.jaccardToDistance(jaccard, kParameter)); // for some reason, the method is 1-based
					distances_containment.set(i + 1, j + 1, Distance.containmentToDistance(containment_i, kParameter));
					distances_containment.set(j + 1, i + 1, Distance.containmentToDistance(containment_j, kParameter));
					distances_mash.setBoth(i + 1, j + 1, Distance.jaccardToMashDistance(jaccard, kParameter));
				}
			}

			logger.info("Exporting...");
			FileWriter outFile = new FileWriter(output, false);
			outFile.write("#nexus\n");
			NexusWriter writer = new NexusWriter();
			writer.write(outFile, taxa, distances_jaccard);
			outFile.close();

			outFile = new FileWriter(output + ".containment", false);
			outFile.write("#nexus\n");
			writer = new NexusWriter();
			writer.write(outFile, taxa, distances_containment);
			outFile.close();

			outFile = new FileWriter(output + ".mash", false);
			outFile.write("#nexus\n");
			writer = new NexusWriter();
			writer.write(outFile, taxa, distances_mash);
			outFile.close();

			if (emptyIntersectionLog != null) {
				outFile = new FileWriter(output + ".log", false);
				for (String line : emptyIntersectionLog) {
					outFile.write(line);
				}
				outFile.close();
			}

		} catch (Exception e) {
			System.out.println("well, f****");
			e.printStackTrace();
		}
	}
}
