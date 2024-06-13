package org.husonlab.fmhdist.db;

import org.junit.Test;

import java.io.File;
import java.sql.SQLException;

import static org.junit.Assert.assertTrue;

public class ReferenceDatabaseTests {
	@Test
	public void createReferenceDb() throws SQLException {
		ReferenceDatabase.create("target/test.db");
		assertTrue(new File("target/test.db").exists());
		(new File("target/test.db")).delete();
	}
}
