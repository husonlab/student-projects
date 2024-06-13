package org.husonlab.fmhdist.db;

import java.io.Closeable;
import java.io.IOException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

import org.husonlab.fmhdist.ncbi.Genome;
import org.husonlab.fmhdist.ncbi.Taxon;
import org.husonlab.fmhdist.ncbi.TaxonomyTree;
import org.husonlab.fmhdist.sketch.FracMinHashSketch;
import org.husonlab.fmhdist.sketch.GenomeSketch;
import org.sqlite.SQLiteConfig;
import org.sqlite.SQLiteException;

import jloda.thirdparty.HexUtils;

public class ReferenceDatabase implements Closeable {
	private Connection connection;
	private Logger logger;

	public static ReferenceDatabase create(String path) throws SQLException {
		SQLiteConfig config = new SQLiteConfig();
		ReferenceDatabase result = new ReferenceDatabase();
		result.connection = config.createConnection("jdbc:sqlite:" + path);
		result.connection.createStatement().execute("CREATE TABLE bloom_filters (taxon_id INTEGER PRIMARY KEY, bloom_filter TEXT NOT NULL) WITHOUT ROWID;");
		result.connection.createStatement().execute("CREATE TABLE tree (key TEXT PRIMARY KEY, value TEXT NOT NULL) WITHOUT ROWID;");
		result.connection.createStatement().execute("CREATE TABLE taxa (taxon_id INTEGER PRIMARY KEY, taxon_name TEXT, taxon_display_name TEXT, parent_id INTEGER REFERENCES taxa(taxon_id)) WITHOUT ROWID;");
		result.connection.createStatement().execute("CREATE TABLE info (key TEXT PRIMARY KEY, value TEXT NOT NULL) WITHOUT ROWID;");
		result.connection.createStatement().execute("CREATE TABLE genomes (taxon_id INTEGER PRIMARY KEY, genome_accession TEXT NOT NULL, genome_size INTEGER, fasta_url TEXT) WITHOUT ROWID;");
		result.connection.createStatement().execute("CREATE TABLE frac_min_hash_sketches (taxon_id INTEGER PRIMARY KEY, frac_min_hash_sketch TEXT NOT NULL) WITHOUT ROWID;");
		return result;
	}

	public static ReferenceDatabase open(String path) throws SQLException, SQLiteException {
		SQLiteConfig config = new SQLiteConfig();
		ReferenceDatabase result = new ReferenceDatabase();
		result.connection = config.createConnection("jdbc:sqlite:" + path);
		return result;
	}

	private ReferenceDatabase() {
		this.logger = Logger.getLogger(ReferenceDatabase.class.getName());
	}

	public Collection<GenomeSketch> getSketches() throws SQLException, IOException {
		this.logger.fine("Getting genomes and sketches from DB...");
		Collection<GenomeSketch> result = new ArrayList<>();
		PreparedStatement view_statement = this.connection.prepareStatement("SELECT b.genome_size as genome_size, a.frac_min_hash_sketch as sketch, a.taxon_id as taxon_id, b.genome_accession, b.fasta_url, c.taxon_name FROM frac_min_hash_sketches as a JOIN genomes as b ON a.taxon_id = b.taxon_id JOIN taxa as c ON a.taxon_id = c.taxon_id;");
		ResultSet rs = view_statement.executeQuery();
		while (rs.next()) {
			byte[] serializedSketch = HexUtils.decodeHexString(rs.getString("sketch"));
			int taxonId = rs.getInt("taxon_id");
			long genomeSize = rs.getLong("genome_size");
			String accession = rs.getString("genome_accession");
			String taxonName = rs.getString("taxon_name");
			String fastaUrl = rs.getString("fasta_url");
			Genome currentGenome = new Genome(taxonName, accession, taxonId, taxonName, fastaUrl, genomeSize);
			FracMinHashSketch sketch = FracMinHashSketch.parse(serializedSketch);
			sketch.setName(taxonName);
			result.add(new GenomeSketch(currentGenome, sketch));
		}

		this.logger.fine("Finished getting genomes and sketches from DB!");
		return result;
	}

	public void insertSketches(Collection<GenomeSketch> sketches) throws SQLException {
		this.logger.fine("Inserting genomes and sketches in DB...");
		PreparedStatement genome_statement = this.connection.prepareStatement("INSERT INTO genomes (taxon_id, genome_accession, genome_size, fasta_url) VALUES (?, ?, ?, ?);");
		PreparedStatement sketch_statement = this.connection.prepareStatement("INSERT INTO frac_min_hash_sketches (taxon_id, frac_min_hash_sketch) VALUES (?, ?);");
		PreparedStatement check_unique_statement = this.connection.prepareStatement("SELECT 1 FROM genomes WHERE taxon_id = ?");
		for (GenomeSketch g : sketches) {
			// Is Taxon already known?
			check_unique_statement.setInt(1, g.getGenome().getTaxonId());
			ResultSet rs = check_unique_statement.executeQuery();
			if (rs.next()) {
				this.logger.fine(String.format("multiple accession codes for taxon %d found. Skipping...", g.getGenome().getTaxonId()));
				continue;
			}

			genome_statement.setInt(1, g.getGenome().getTaxonId());
			genome_statement.setString(2, g.getGenome().getAccession());
			genome_statement.setLong(3, g.getGenome().getGenomeSize());
			genome_statement.setString(4, g.getGenome().getFastaUrl());
			genome_statement.executeUpdate();

			sketch_statement.setInt(1, g.getGenome().getTaxonId());
			sketch_statement.setString(2, HexUtils.encodeHexString(g.getSketch().getBytes()));
			sketch_statement.executeUpdate();
		}
		this.logger.fine("Finished inserting genome sketches in DB!");
	}

	public void insertTaxonomy(TaxonomyTree taxonomy) throws SQLException {
		this.logger.fine("Inserting taxa in DB...");
		PreparedStatement s = this.connection.prepareStatement("INSERT INTO taxa (taxon_id, taxon_name, taxon_display_name, parent_id) VALUES (?, ?, ?, ?)");
		for (Taxon t : taxonomy.getTaxa().values()) {
			s.setInt(1, t.getTaxonId());
			s.setString(2, t.getOrganismName());
			s.setString(3, t.getOrganismName());
			int parent = 0;
			if (!taxonomy.getRoot().equals(t)) {
				for (Taxon p : taxonomy.getTree().predecessors(t)) { // there is only one
					parent = p.getTaxonId();
				}
			}
			s.setInt(4, parent);
			s.executeUpdate();
		}
		this.logger.fine("Finished inserting taxa in DB!");
	}

	public void insertFullInfo(int kSize, int sParam, int seed, String hashFunction) throws SQLException {
		this.logger.fine("Inserting sketch creation info...");
		PreparedStatement s = this.connection.prepareStatement("INSERT INTO info (key, value) VALUES (?, ?);");
		int[] values = new int[]{kSize, sParam, seed};
		String[] keys = new String[]{"sketch_k", "sketch_s", "sketch_seed"};
		for (int i = 0; i < 3; i++) {
			s.setString(1, keys[i]);
			s.setInt(2, values[i]);
			s.executeUpdate();
		}
		s.setString(1, "hash_function");
		s.setString(2, hashFunction);
		s.executeUpdate();
	}

	public Map<String, Integer> getNumericalInfo() throws SQLException {
		this.logger.fine("Fetching sketch creation info from DB...");
		Map<String, Integer> result = new HashMap<>();
		PreparedStatement s = this.connection.prepareStatement("SELECT key, value FROM info;");
		ResultSet rs = s.executeQuery();
		while (rs.next()) {
			result.put(rs.getString("key"), rs.getInt("value"));
		}
		return result;
	}

	public String getUsedHashFunction() throws SQLException {
		this.logger.fine("Fetching used hash function from DB...");
		PreparedStatement s = this.connection.prepareStatement("SELECT value FROM info WHERE key=?");
		s.setString(1, "hash_function");
		ResultSet rs = s.executeQuery();
		while (rs.next()) {
			return rs.getString("value");
		}
		return "";
	}

	@Override
	public void close() throws IOException {
		try {
			this.connection.close();
		} catch (SQLException e) {
			throw new IOException(e);
		}
	}
}