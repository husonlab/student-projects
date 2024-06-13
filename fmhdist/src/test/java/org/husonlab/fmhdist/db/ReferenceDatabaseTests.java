package org.husonlab.fmhdist.db;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.sql.SQLException;

import org.junit.Test;

public class ReferenceDatabaseTests {
	@Test
	public void createReferenceDb() throws SQLException {
		ReferenceDatabase.create("target/test.db");
		assertTrue(new File("target/test.db").exists());
		(new File("target/test.db")).delete();
	}
}
