package org.husonlab.fmhdist.cmd;

import jloda.thirdparty.HexUtils;
import jloda.util.FileLineIterator;
import org.husonlab.fmhdist.db.ReferenceDatabase;
import org.husonlab.fmhdist.sketch.Distance;
import org.husonlab.fmhdist.sketch.FracMinHashSketch;
import org.husonlab.fmhdist.sketch.GenomeSketch;
import org.husonlab.fmhdist.sketch.IncompatibleParameterException;
import org.husonlab.fmhdist.util.HashFunctionParser;
import org.sqlite.SQLiteException;
import splitstree6.data.DistancesBlock;
import splitstree6.data.TaxaBlock;
import splitstree6.io.writers.distances.NexusWriter;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.*;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 * Calculates the distances to the closest genomes in a reference database
 */
public class ReferenceDistanceCalculator {
	private int sParameter;
	private int kParameter;
	private int randomSeed;
	private long hashedMagicNumber;

	private FracMinHashSketch lineToSketch(String line) {
		try {
			String[] splitLine = line.replaceAll("\\s+", "").split(",");
			byte[] content = Files.readAllBytes(Paths.get(splitLine[0]));
			FracMinHashSketch result = FracMinHashSketch.parse(HexUtils.decodeHexString(new String(content)));
			if (splitLine.length > 1) {
				result.setName(splitLine[1]);
			} else {
				result.setName(splitLine[0]);
			}
			return result;
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	private List<FracMinHashSketch> prepareGenomesFromDatabase(String database)
			throws SQLException, IncompatibleParameterException, IOException {
		ReferenceDatabase db = ReferenceDatabase.open(database);
		Map<String, Integer> info = db.getNumericalInfo();
		String hashFunctionName = db.getUsedHashFunction();
		Collection<GenomeSketch> refSketches = db.getSketches();
		db.close();

		if (!info.containsKey("sketch_k") ||
			!info.containsKey("sketch_s") ||
			!info.containsKey("sketch_seed") ||
			hashFunctionName.equals("")) {
			throw new IncompatibleParameterException(
					"reference db does not provide all sketching parameters (s, k, seed, hash function)");
		}

		if (!HashFunctionParser.getSupportedFunctions().contains(hashFunctionName)) {
			throw new IncompatibleParameterException(
					String.format("hash function '%s' used in database is not supported", hashFunctionName));
		}
		this.kParameter = info.get("sketch_k");
		this.sParameter = info.get("sketch_s");
		this.randomSeed = info.get("sketch_seed");
		this.hashedMagicNumber = FracMinHashSketch
				.getHashedMagicNumber(HashFunctionParser.createHashFunction(hashFunctionName, randomSeed));

		List<FracMinHashSketch> result = new ArrayList<>();
		for (GenomeSketch g : refSketches) {
			result.add(g.getSketch());
		}
		return result;
	}

	private List<FracMinHashSketch> prepareGenomesFromSketchList(String database) throws IOException {
		FileLineIterator it = new FileLineIterator(database);
		List<FracMinHashSketch> sketches = it
				.stream()
				.map(this::lineToSketch)
				.collect(Collectors.toList());
		it.close();

		if (sketches.size() == 0) {
			throw new IOException("reference db is empty");
		}

		FracMinHashSketch first = sketches.get(0);
		this.kParameter = first.getKSize();
		this.sParameter = first.getSParam();
		this.randomSeed = first.getSeed();
		this.hashedMagicNumber = first.getHashedMagicNumber();

		return sketches;
	}

	/**
	 * Calculates the distances from the input query to the closest genomes in
	 * the database. Calculates three sets evolutionary distances: Mash
	 * Distance, FracMinHash containment distance and FracMinHash distance.
	 *
	 * @param input       Path to a CSV containing paths to sketches, one per line.
	 *                    This will be the query sequences.
	 * @param output      Path to the output FracMinHash distance in Nexus format.
	 *                    The other two distances are stored using ".mash" and ".containment"
	 *                    suffixes.
	 * @param database    Path to the database. This can be either an SQLite
	 *                    database or a csv file listing paths to sketches.
	 * @param maxDistance The maximum distance of a reference sequence to a
	 *                    query sequence to be included in the output.
	 */
	public void run(
			String input,
			String output,
			String database,
			double maxDistance) {
		Logger logger = Logger.getLogger(DistanceCalculator.class.getName());
		try {
			logger.info("Loading reference DB...");
			logger.fine("Try to parse as SQLite...");
			List<FracMinHashSketch> refSketches;
			try {
				refSketches = prepareGenomesFromDatabase(database);
			} catch (SQLiteException e) {
				logger.fine("Failed to read SQLite!");
				logger.fine("Try to parse as CSV...");
				refSketches = prepareGenomesFromSketchList(database);
			}

			logger.info("Reading queries list...");
			FileLineIterator it = new FileLineIterator(input);
			List<FracMinHashSketch> sketches = it
					.stream()
					.map(this::lineToSketch)
					.collect(Collectors.toList());
			it.close();

			logger.info("Finding closest reference genomes...");
			Set<FracMinHashSketch> resultSketchSet = new HashSet<>();
			for (FracMinHashSketch querySketch : sketches) {
				if (this.sParameter != querySketch.getSParam() ||
					this.kParameter != querySketch.getKSize() ||
					this.randomSeed != querySketch.getSeed() ||
					this.hashedMagicNumber != querySketch.getHashedMagicNumber()) {
					logger.severe("sketches have incompatible sketching parameters");
					return;
				}

				for (FracMinHashSketch refSketch : refSketches) {
					double jaccard = Distance.calculateJaccardIndex(querySketch.getValues(),
							refSketch.getValues(), sParameter);
					double distance = Distance.jaccardToDistance(jaccard, kParameter);
					if (distance <= maxDistance) {
						resultSketchSet.add(refSketch);
					}
				}
				resultSketchSet.add(querySketch);
			}

			logger.info("Calculating pairwise distances...");
			List<FracMinHashSketch> resultSketchesList = new ArrayList<>(resultSketchSet);

			DistancesBlock distances_jaccard = new DistancesBlock();
			distances_jaccard.setNtax(resultSketchSet.size());

			DistancesBlock distances_mash = new DistancesBlock();
			distances_mash.setNtax(resultSketchSet.size());

			DistancesBlock distances_containment = new DistancesBlock();
			distances_containment.setNtax(resultSketchSet.size());

			TaxaBlock taxa = new TaxaBlock();
			for (int i = 0; i < resultSketchSet.size(); i++) {
				taxa.addTaxonByName(resultSketchesList.get(i).getName());
				for (int j = i; j < resultSketchSet.size(); j++) {
					double jaccard = Distance.calculateJaccardIndex(
							resultSketchesList.get(i).getValues(),
							resultSketchesList.get(j).getValues(),
							sParameter);
					// Containment is not symmetrical
					double containment_i = Distance.calculateContainmentIndex(
							resultSketchesList.get(i).getValues(),
							resultSketchesList.get(j).getValues(),
							sParameter);
					double containment_j = Distance.calculateContainmentIndex(
							resultSketchesList.get(j).getValues(),
							resultSketchesList.get(i).getValues(),
							sParameter);

					// for some reason, the method is 1-based
					distances_jaccard.setBoth(i + 1, j + 1, Distance.jaccardToDistance(jaccard, kParameter));
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
			writer.write(outFile, taxa, distances_containment);
			outFile.close();

		} catch (Exception e) {
			System.out.println("well, f****");
			e.printStackTrace();
		}
	}
}
