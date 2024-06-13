/*
 *  GeneTreeView.java Copyright (C) 2024 Algorithms in Bioinformatics
 *
 *  (Some files contain contributions from other authors, who are then mentioned separately.)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package sprojects.genetreeview;

import javafx.application.Application;
import javafx.fxml.FXMLLoader;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.stage.Stage;
import sprojects.genetreeview.model.Model;

import java.io.IOException;
import java.util.Objects;

public class GeneTreeView extends Application {
	private final Model model = new Model();
	private GeneTreeViewController controller;
	private sprojects.genetreeview.GeneTreeViewPresenter presenter;
	private Parent root;
	private Stage stage;

	public static void main(String[] args) {
		launch(args);
	}

	@Override
	public void start(Stage primaryStage) {
		this.stage = primaryStage;
		var fxmlLoader = new FXMLLoader();
		try (var ins = Objects.requireNonNull(GeneTreeViewController.class.getResource("GeneTreeView.fxml")).openStream()) {
			fxmlLoader.load(ins);
		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
		controller = fxmlLoader.getController();
		root = fxmlLoader.getRoot();

		presenter = new sprojects.genetreeview.GeneTreeViewPresenter(this);

		stage.setScene(new Scene(root));
		stage.getScene().getStylesheets().add("jloda/resources/css/white_pane.css");
		stage.sizeToScene();
		stage.setTitle("GeneTreeViewer");
		stage.show();

	}

	public GeneTreeViewController getController() {
		return controller;
	}

	public GeneTreeViewPresenter getPresenter() {
		return presenter;
	}

	public Parent getRoot() {
		return root;
	}

	public Stage getStage() {
		return stage;
	}

	public Model getModel() {
		return model;
	}
}
