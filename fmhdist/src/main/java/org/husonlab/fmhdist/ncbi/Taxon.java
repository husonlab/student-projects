package org.husonlab.fmhdist.ncbi;

public class Taxon {
	private String organismName;
	private int taxonId;

	public Taxon(String organismName, int taxonId) {
		this.organismName = organismName;
		this.taxonId = taxonId;
	}

	public String getOrganismName() {
		return this.organismName;
	}

	public void setOrganismName(String organismName) {
		this.organismName = organismName;
	}

	public int getTaxonId() {
		return this.taxonId;
	}

	public void setTaxonId(int taxonId) {
		this.taxonId = taxonId;
	}

	@Override
	public boolean equals(Object o) {
		if (o == this)
			return true;
		if (!(o instanceof Taxon)) {
			return false;
		}
		Taxon taxon = (Taxon) o;
		return this.taxonId == taxon.taxonId;
	}

	@Override
	public int hashCode() {
		return Integer.hashCode(this.taxonId);
	}

}
