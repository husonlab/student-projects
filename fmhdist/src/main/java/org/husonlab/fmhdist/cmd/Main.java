package org.husonlab.fmhdist.cmd;

import org.husonlab.fmhdist.util.HashFunctionParser;

import jloda.fx.util.ArgsOptions;
import jloda.util.ProgramExecutorService;
import jloda.util.UsageException;
import net.openhft.hashing.LongHashFunction;

public class Main {
	private final static String CREATE_DB_COMMAND = "db";
	private final static String COMPARE_SKETCH_COMMAND = "dist";
	private final static String COMPARE_REF_SKETCH_COMMAND = "ref_dist";
	private final static String SKETCH_COMMAND = "sketch";
	private final static String OUTLINE_COMMAND = "outline";

	public static void main(String[] args) throws UsageException {
		final ArgsOptions options = new ArgsOptions(args, Main.class,
				"Tool to calculate and compare FracMinHash sketches for sequence files");

		final String command = options.getCommand(
				new ArgsOptions.Command(
						CREATE_DB_COMMAND,
						"Create a new reference database from the given NCBI accession codes."
				),
				new ArgsOptions.Command(
						COMPARE_SKETCH_COMMAND,
						"Calculate the pairwise distances for all query sketches."
				),
				new ArgsOptions.Command(
						COMPARE_REF_SKETCH_COMMAND,
						"Find the closest database entries (in terms of distance) to given query sequences\n"
						+ "and calculate all pairwise distances for those + the query sequences."
				),
				new ArgsOptions.Command(
						SKETCH_COMMAND,
						"Calculate the sketch for all given sequences and store them on the file system."),
				new ArgsOptions.Command(
						OUTLINE_COMMAND,
						"Calculates the phylogenetic outline based on the given distances. The output is stored in SVG format.")

		);

		options.comment("Input/Output options");
		final String input = options.getOptionMandatory(
				"-i",
				"input",
				String.format(
						"CSV of" +
						"(a - for %s command) NCBI accession codes, a code per line\n" +
						"(b - for %s command) sequence file paths or URLs to fasta files (gzip ok), " +
						"a path per line. Optional name of the resulting sketch can be separated by comma\n" +
						"(c -  for %s and %s command) sketch file paths, a path per line. Optional name of " +
						"taxon in the distance matrix can be separated by comma. For the %s command, this is expected to be a distance mattrix in Nexus format",
						CREATE_DB_COMMAND, SKETCH_COMMAND, COMPARE_REF_SKETCH_COMMAND, COMPARE_SKETCH_COMMAND, OUTLINE_COMMAND),
				""
		);

		final String database = options.getOption(
				"-db",
				"database",
				String.format(
						"Path to reference database file (input for %s command). Can be either a " +
						"SQLite DB or a CSV (see --input param for more details)",
						COMPARE_REF_SKETCH_COMMAND
				),
				"database.db");

		final String output = options.getOption(
				"-o",
				"output",
				"Path to the output file (distances, sketches, coordinates, database)",
				""
		);

		final boolean saveCoordinates = options.getOption(
				"-c",
				"saveCoordinates",
				String.format(
						"When running %s, control if the coordinates of the k-mers that " +
						"are part of the sketch should be saved using <output>.coordinates.",
						SKETCH_COMMAND
				),
				false
		);

		options.comment("Algorithm parameters");
		final int kParameter = options.getOption(
				"-k",
				"kmerSize",
				"Word size k",
				21
		);
		final int sParameter = options.getOption(
				"-s",
				"scalingFactor",
				"Scaling factor s. Hash values h are only part of the sketch if h <= H/s",
				2000
		);
		final int randomSeed = options.getOption(
				"-rs",
				"randomSeed",
				"Seed for the selected hash function",
				42
		);
		final String hashFunctionName = options.getOption(
				"-hf",
				"hashFunction",
				"The hash function to use to calculate the 64-bit hash for each k-mer",
				HashFunctionParser.getSupportedFunctions(),
				HashFunctionParser.FARM_HASH_NAME
		);
		final double maxDistance = options.getOption(
				"-md",
				"maxDistance",
				"The maximum distance that a query sequence should have to have to a " +
				"reference sequence for the reference sequence to be included in the output",
				0.4
		);

		options.comment("Visualization options");
		final int width = options.getOption("-iw", "imageWidth", "The width of the generated outline image", 1000);
		final int height = options.getOption("-ih", "imageHeight", "The height of the generated outline image", 1000);
		final int scale = options.getOption("-is", "imageScale", "The scaling factor of the generated outline image", 1000);
		final int xOffset = options.getOption("-ix", "imageXOffset", "The x offset of the generated outline image", 500);
		final int yOffset = options.getOption("-iy", "imageYOffset", "The y offset of the generated outline image", 500);
		final String labels = options.getOption(
				"-il",
				"imageLabels",
				"Path to a TSV listing alternative labels for the generated outline. If not set, taxa name from distance " +
				"file will be used. If set but taxa are missing, their taxon name from the distance file will be used. List " +
				"taxa without value to hide the label in the output.",
				""
		);

		options.comment("Performance options");
		ProgramExecutorService
				.setNumberOfCoresToUse(options.getOption("-t", "threads", "Number of threads", 1));

		options.done();

		LongHashFunction hashFunction = HashFunctionParser.createHashFunction(hashFunctionName, randomSeed);
		switch (command) {
			case CREATE_DB_COMMAND:
				DatabaseCreator dbCreator = new DatabaseCreator();
				dbCreator.run(input, output, kParameter, sParameter, hashFunction, hashFunctionName, randomSeed);
				break;
			case COMPARE_SKETCH_COMMAND:
				DistanceCalculator distanceCalculator = new DistanceCalculator();
				distanceCalculator.run(input, output);
				break;
			case COMPARE_REF_SKETCH_COMMAND:
				ReferenceDistanceCalculator refDistanceCalculator = new ReferenceDistanceCalculator();
				refDistanceCalculator.run(input, output, database, maxDistance);
				break;
			case SKETCH_COMMAND:
				SequenceSketcher sketcher = new SequenceSketcher();
				sketcher.run(input, output, kParameter, sParameter, hashFunction, randomSeed, saveCoordinates);
				break;
			case OUTLINE_COMMAND:
				OutlineVisualizer visualizer = new OutlineVisualizer();
				visualizer.run(input, output, width, height, scale, xOffset, yOffset, labels);
		}
	}
}
