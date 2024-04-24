module sprojects {
	requires transitive jloda_core;
	requires transitive jloda_fx;
	requires transitive splitstreesix;
	requires transitive javafx.controls;
	requires transitive javafx.graphics;
	requires transitive javafx.fxml;
	requires javafx.base;

	requires org.apache.commons.collections4;
	requires org.apache.commons.math4.legacy;
	requires org.apache.commons.math4.legacy.exception;

	requires countryboundaries;
	requires com.install4j.runtime;
	requires java.sql;
	requires java.sql.rowset;
	// requires org.xerial.sqlitejdbc;
	requires java.desktop;
	requires org.fxmisc.flowless;
	requires org.fxmisc.richtext;
	requires org.fxmisc.undo;
	requires com.google.zxing;
	requires org.xerial.sqlitejdbc;

	opens sprojects.mapview;
	opens sprojects.mapview.nodes;
	opens sprojects.mapview.mapbuilder;
	exports sprojects.mapview;

	opens sprojects.genetreeview;
	opens sprojects.genetreeview.layout;
	opens sprojects.genetreeview.io;
	opens sprojects.genetreeview.model;
	opens sprojects.genetreeview.util;
	exports sprojects.genetreeview;
}