package org.husonlab.fmhdist.cmd;

import jloda.util.FileLineIterator;
import jloda.util.ProgramExecutorService;
import jloda.util.Single;
import net.openhft.hashing.LongHashFunction;
import org.husonlab.fmhdist.db.ReferenceDatabase;
import org.husonlab.fmhdist.ncbi.Genome;
import org.husonlab.fmhdist.ncbi.NcbiApi;
import org.husonlab.fmhdist.ncbi.TaxonomyTree;
import org.husonlab.fmhdist.sketch.GenomeSketch;

import java.io.File;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 * Class to create a reference database given a set of NCBI accession codes.
 */
public class DatabaseCreator {
	private static NcbiApi api = new NcbiApi("https://api.ncbi.nlm.nih.gov/datasets/v2alpha");

	/**
	 * Creates a new reference database.
	 *
	 * @param input            path to a CSV, each line stating an NCBI accession code
	 * @param output           path to the output database
	 * @param kParameter       k-mer size for FracMinHash
	 * @param sParameter       scaling parameters for FracMinHash
	 * @param hashFunction     the hash function that should be used to create the FracMinHash sketch
	 * @param hashFunctionName the name of that hash function, it will be written to the database
	 * @param randomSeed       the random seed that was applied to create the hash function
	 */
	public void run(
			String input,
			String output,
			int kParameter,
			int sParameter,
			LongHashFunction hashFunction,
			String hashFunctionName,
			int randomSeed) {
		Logger logger = Logger.getLogger(DatabaseCreator.class.getName());

		if ((new File(output)).exists()) {
			(new File(output)).delete();
		}

		try {
			FileLineIterator it = new FileLineIterator(input);
			List<String> accessionCodes = it.stream().map(line -> line.replaceAll("\\s+", ""))
					.collect(Collectors.toList());
			it.close();

			logger.info("Preparing genomes...");
			List<Genome> genomes = api.getGenomes(accessionCodes);

			logger.info("Preparing taxonomy tree...");
			TaxonomyTree tree = api.getTaxonomyTreeForGenomes(genomes);

			logger.info("Sketching sequences...");
			Queue<GenomeSketch> sketches = new ConcurrentLinkedQueue<>();
			final Single<Throwable> exception = new Single<>();
			final ExecutorService executor = Executors
					.newFixedThreadPool(ProgramExecutorService.getNumberOfCoresToUse());

			try {
				genomes.forEach(genome -> executor.submit(() -> {
					if (exception.isNull()) {
						int retries = 5;
						// Sometimes, the connection to NCBI breaks - this is a quick workaround
						while (retries-- > 0) {
							try {
								GenomeSketch sketch = GenomeSketch.sketch(genome, kParameter, sParameter, hashFunction, randomSeed, false);
								sketches.add(sketch);
								break;
							} catch (Exception ex) {
								logger.warning(ex.getMessage());
							} catch (Throwable e) {
								// Somethings wrong here - no way to recover.
								logger.severe(e.getMessage());
								exception.setIfCurrentValueIsNull(e);
								break;
							}
						}
					}
				}));
			} finally {
				executor.shutdown();
				executor.awaitTermination(1000, TimeUnit.DAYS);
			}

			logger.info("Exporting database...");
			ReferenceDatabase db = ReferenceDatabase.create(output);
			db.insertTaxonomy(tree);
			db.insertSketches(sketches);
			db.insertFullInfo(kParameter, sParameter, randomSeed, hashFunctionName);
			db.close();
		} catch (Exception e) {
			System.out.println("well, f****");
			e.printStackTrace();
		}
	}
}
