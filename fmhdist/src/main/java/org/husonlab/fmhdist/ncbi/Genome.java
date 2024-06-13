package org.husonlab.fmhdist.ncbi;

import jloda.util.FileUtils;

public class Genome {
	private String organismName;
	private String accession;
	private int taxonId;
	private String assemblyName;
	private String downloadLink;
	private long genomeSize;
	private boolean isDraft;

	public Genome(String organismName, String accession, int taxonId, String assemblyName, String downloadLink, long genomeSize) {
		this.organismName = organismName;
		this.accession = accession;
		this.taxonId = taxonId;
		this.assemblyName = assemblyName;
		this.downloadLink = downloadLink;
		this.genomeSize = genomeSize;
		this.isDraft = false;
	}

	public Genome(String organismName, String path) {
		this.organismName = organismName;
		this.accession = organismName;
		this.assemblyName = organismName;
		this.downloadLink = path;

		// Estimate will include the headers as well - but this is a good enough
		// approximation for now.
		this.genomeSize = FileUtils.guessUncompressedSizeOfFile(path);
		this.isDraft = true;
	}

	public String getOrganismName() {
		return this.organismName;
	}

	public void setOrganismName(String organismName) {
		this.organismName = organismName;
	}

	public String getAccession() {
		return this.accession;
	}

	public void setAccession(String accession) {
		this.accession = accession;
	}

	public int getTaxonId() {
		return this.taxonId;
	}

	public void setTaxonId(int taxonId) {
		this.taxonId = taxonId;
	}

	public String getAssemblyName() {
		return this.assemblyName;
	}

	public void setAssemblyName(String assemblyName) {
		this.assemblyName = assemblyName;
	}

	public String getDownloadLink() {
		return this.downloadLink;
	}

	public void setDownloadLink(String downloadLink) {
		this.downloadLink = downloadLink;
	}

	public String getFastaUrl() {
		if (this.isDraft) {
			return this.downloadLink;
		}
		return this.downloadLink + "/" + this.accession + "_" + this.assemblyName.replaceAll("\\s+|\\/", "_") + "_genomic.fna.gz";
	}

	public void setGenomeSize(int genomeSize) {
		this.genomeSize = genomeSize;
	}

	public long getGenomeSize() {
		return this.genomeSize;
	}

	public boolean isDraft() {
		return this.isDraft;
	}

}
