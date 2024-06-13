module sprojects {
	requires transitive jloda_core;
	requires transitive jloda_fx;
	requires transitive splitstreesix;
	requires transitive javafx.controls;
	requires transitive javafx.graphics;
	requires transitive javafx.fxml;
	requires javafx.base;

	opens sprojects.genetreeview;
	opens sprojects.genetreeview.layout;
	opens sprojects.genetreeview.io;
	opens sprojects.genetreeview.model;
	opens sprojects.genetreeview.util;
	exports sprojects.genetreeview;
}