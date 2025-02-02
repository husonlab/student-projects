/*
 *  MapViewController.java Copyright (C) 2024 Algorithms in Bioinformatics
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

package sprojects.mapview;

import javafx.fxml.FXML;
import javafx.scene.control.*;
import javafx.scene.layout.Pane;
import javafx.scene.layout.StackPane;

public class MapViewController {

	@FXML
	private MenuItem closeMenuItem;

	@FXML
	private MenuItem copyMenuItem;

	@FXML
	private Pane mainPane;

	@FXML
	private MenuBar menuBar;

	@FXML
	private MenuItem openMenuItem;

	@FXML
	private StackPane stackPane;

	@FXML
	private ProgressBar progressBar;
	@FXML
	private Label label;

	@FXML
	private Button redrawButton;

	@FXML
	private Slider chartSizeSlider;

	@FXML
	private Label infoLabel;

	@FXML
	private CheckBox checkBoxLegend;

	@FXML
	private ChoiceBox<String> choiceBoxColorScheme;

	@FXML
	private CheckBox showLabelsBox;

	@FXML
	private void initialize() {
		progressBar.setVisible(false);
	}

	public CheckBox getCheckBoxLegend() {
		return checkBoxLegend;
	}

	public CheckBox getShowLabelsBox() {
		return showLabelsBox;
	}

	public ChoiceBox<String> getChoiceBoxColorScheme() {
		return choiceBoxColorScheme;
	}

	public Label getInfoLabel() {
		return infoLabel;
	}

	public MenuItem getCloseMenuItem() {
		return closeMenuItem;
	}

	public Slider getChartSizeSlider() {
		return chartSizeSlider;
	}

	public MenuItem getCopyMenuItem() {
		return copyMenuItem;
	}

	public Pane getMainPane() {
		return mainPane;
	}

	public Label getLabel() {
		return label;
	}

	public MenuItem getOpenMenuItem() {
		return openMenuItem;
	}

	public StackPane getStackPane() {
		return stackPane;
	}

	public ProgressBar getProgressBar() {
		return progressBar;
	}

	public Button getRedrawButton() {
		return redrawButton;
	}
}
