import os
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
from PyQt5 import uic
from PyQt5.QtGui import QImage, QPixmap
from PyQt5.QtWidgets import QDialog, QGraphicsPixmapItem, QGraphicsScene
from PyQt5.QtCore import Qt

FORM_CLASS, _ = uic.loadUiType(os.path.join(os.path.dirname(__file__), 'qcrocoflow_Bathyview.ui'))

class BathyView(QDialog, FORM_CLASS):
    def __init__(self, parent=None):
        """Initialize the BathyView class."""

        super(BathyView, self).__init__(parent)
        self.setupUi(self)

        # Connect signals to slots
        self.AdjustBathydoubleSpinBox.valueChanged.connect(self.update_view)
        self.SavenewlvlPushButton.clicked.connect(self.save_new_level)
        self.greaterRadioButton.toggled.connect(self.set_radio_button_value)

        # Set alignment
        self.graphicsView.setAlignment(Qt.AlignCenter)

        # Initialize instance variables
        self.h = None
        self.nc = None
        self.previous_adjust_value = 0

    def set_radio_button_value(self):
        """Reverse h when the radio button state changes"""

        self.h = self.h * -1 if self.greaterRadioButton.isChecked() else self.original_h.copy()
        self.update_view()

    def set_h(self, h, nc):
        """Set h and nc, then update the view."""

        self.h = h
        self.original_h = h.copy()
        self.nc = nc
        self.update_view()

    def save_new_level(self):
        """Save the new level to the NetCDF file and close the dialog,
           or print a warning message if there is no data to save."""

        if self.h is not None and self.nc is not None:
            self.nc.variables['h'][:] = self.h
            print(self.h)
            self.close()
            return self.h
        else:
            print("\nNo data to save...")
            return None

    def update_view(self):
        """Update the view with the current h values and adjustment."""

        adjust_value = self.AdjustBathydoubleSpinBox.value()

        if self.h is not None:
            vmin = np.min(self.h)
            vmax = np.max(self.h)

            # Undo the previous adjustment
            self.h -= self.previous_adjust_value

            # Adjust values of h
            self.h += adjust_value

            # Store adjust_value for next time
            self.previous_adjust_value = adjust_value

            # Plot the data
            fig = plt.Figure()
            canvas = FigureCanvas(fig)
            ax = fig.add_subplot(111)

            # Set up color map and normalization
            colors = [(0, "lightblue"), (0.5, "blue"), (1, "black")]
            cmap = LinearSegmentedColormap.from_list("mycmap", colors)
            norm = BoundaryNorm([vmin, 0, vmax], cmap.N)

            # Display the image
            im = ax.imshow(self.h, origin='lower', cmap=cmap, norm=norm)
            fig.colorbar(im)


            # Render the image
            canvas.draw()
            width, height = map(int, fig.get_size_inches() * fig.get_dpi())  # Convertir en entiers
            buf = canvas.buffer_rgba()
            data = buf.tobytes()  # Convertir en bytes
            qimage = QImage(data, width, height, QImage.Format_RGBA8888)
            pixmap = QPixmap(qimage)


            # Display the image in a QGraphicsScene
            scene = QGraphicsScene()
            item = QGraphicsPixmapItem(pixmap)
            scene.addItem(item)
            self.graphicsView.setScene(scene)