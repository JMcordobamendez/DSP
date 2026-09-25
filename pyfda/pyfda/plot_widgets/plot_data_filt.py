# -*- coding: utf-8 -*-
#
# This file is part of the pyfda project hosted at https://github.com/chipmuenk/pyfda
#
# Copyright © pyfda Project Contributors
# Licensed under the terms of the MIT License
# (see file LICENSE in root directory for details)

"""
Widget for applying the current filter design to measured data loaded from a
file (csv, txt, wav, npy, ...) and plotting original and filtered data overlaid.
"""
import csv
import logging
import os
import re

import numpy as np
import scipy.signal as sig

from pyfda.plot_widgets.mpl_widget import MplWidget
from pyfda.pyfda_rc import params
from pyfda.filterbroker import fb_get
import pyfda.libs.pyfda_dirs as dirs
from pyfda.libs import pyfda_io_lib as io
from pyfda.libs.compat import (
    QCheckBox, QWidget, QFrame, QComboBox, QLabel, QPushButton, QHBoxLayout, pyqtSignal)

logger = logging.getLogger(__name__)

# Dict containing class name : display name
classes = {'PlotDataFilt': 'Data Filt'}

FILE_TYPES = ('csv', 'txt', 'wav', 'npy')


class PlotDataFilt(QWidget):
    """
    Load a column of data from a file, filter it with the current design and
    plot original and filtered data (and optionally their spectra) overlaid.
    """
    sig_rx = pyqtSignal(object)  # incoming

    def __init__(self):
        super().__init__()
        self.needs_calc = True   # flag whether plot needs to be recalculated
        self.tool_tip = self.tr("Filter data from a file with the current design")
        self.tab_label = self.tr("Data Filt")

        self.data = None        # 2D array (samples, columns) of loaded data
        self.col_names = []     # names of the data columns
        self.file_name = ""
        self.x = None           # selected data column
        self.y = None           # filtered data (None: not filtered yet)

        self._construct_ui()

    # ------------------------------------------------------------------------------
    def _construct_ui(self):
        """
        Initialize the widget, consisting of:
        - Matplotlib widget with NavigationToolbar
        - Frame with control elements
        """
        self.but_load = QPushButton(self.tr("Load data ..."), self)
        self.but_load.setToolTip(self.tr(
            "<span>Load data from a file (csv, txt, wav, npy). CSV import "
            "uses the settings of the CSV options in the coefficient tab.</span>"))
        self.lbl_file = QLabel(self.tr("No file loaded"), self)

        self.lbl_col = QLabel(self.tr("Column:"), self)
        self.cmb_col = QComboBox(self)
        self.cmb_col.setToolTip(self.tr("<span>Data column to be filtered.</span>"))
        self.cmb_col.setSizeAdjustPolicy(QComboBox.AdjustToContents)

        self.lbl_time = QLabel(self.tr("Time:"), self)
        self.cmb_time = QComboBox(self)
        self.cmb_time.setToolTip(self.tr(
            "<span>Column with the time axis, or '<i>n / f_S</i>' to calculate the "
            "time axis from the sampling frequency of the current design.</span>"))
        self.cmb_time.setSizeAdjustPolicy(QComboBox.AdjustToContents)

        self.but_filter = QPushButton(self.tr("Filter data"), self)
        self.but_filter.setToolTip(self.tr(
            "<span>Apply the current filter design to the selected data column.</span>"))
        self.but_filter.setEnabled(False)

        self.chk_zero_phase = QCheckBox(self.tr("Zero phase"), self)
        self.chk_zero_phase.setToolTip(self.tr(
            "<span>Filter forwards and backwards (<i>filtfilt</i>) for zero phase "
            "delay. The effective magnitude response is squared.</span>"))

        self.chk_spectrum = QCheckBox(self.tr("Spectrum"), self)
        self.chk_spectrum.setToolTip(self.tr(
            "<span>Show the magnitude spectra of original and filtered data "
            "in a second plot.</span>"))

        self.but_export = QPushButton(self.tr("Export ..."), self)
        self.but_export.setToolTip(self.tr(
            "<span>Save original and filtered data as a CSV file.</span>"))
        self.but_export.setEnabled(False)

        lay_h_controls = QHBoxLayout()
        lay_h_controls.addWidget(self.but_load)
        lay_h_controls.addWidget(self.lbl_file)
        lay_h_controls.addWidget(self.lbl_col)
        lay_h_controls.addWidget(self.cmb_col)
        lay_h_controls.addWidget(self.lbl_time)
        lay_h_controls.addWidget(self.cmb_time)
        lay_h_controls.addStretch(10)
        lay_h_controls.addWidget(self.chk_zero_phase)
        lay_h_controls.addWidget(self.but_filter)
        lay_h_controls.addWidget(self.chk_spectrum)
        lay_h_controls.addWidget(self.but_export)

        # This widget encompasses all control subwidgets:
        self.frm_controls = QFrame(self, objectName="frm_controls")
        self.frm_controls.setLayout(lay_h_controls)

        self.mplwidget = MplWidget(self)
        self.mplwidget.lay_v_main_mpl.addWidget(self.frm_controls)
        self.mplwidget.lay_v_main_mpl.setContentsMargins(*params['mpl_margins'])
        self.mplwidget.mpl_toolbar.a_he.setEnabled(False)
        self.mplwidget.mpl_toolbar.a_ui_num_levels = 2
        self.setLayout(self.mplwidget.lay_v_main_mpl)

        self.draw()

        # ----------------------------------------------------------------------
        # GLOBAL SIGNALS & SLOTs
        # ----------------------------------------------------------------------
        self.sig_rx.connect(self.process_sig_rx)
        # ----------------------------------------------------------------------
        # LOCAL SIGNALS & SLOTs
        # ----------------------------------------------------------------------
        self.mplwidget.mpl_toolbar.sig_tx.connect(self.process_sig_rx)
        self.but_load.clicked.connect(self.load_data)
        self.cmb_col.currentIndexChanged.connect(self.select_data)
        self.cmb_time.currentIndexChanged.connect(self.draw)
        self.but_filter.clicked.connect(self.filter_data)
        self.chk_zero_phase.clicked.connect(self._refilter)
        self.chk_spectrum.clicked.connect(self.draw)
        self.but_export.clicked.connect(self.export_data)

    # ------------------------------------------------------------------------------
    def process_sig_rx(self, dict_sig=None):
        """
        Process signals coming from the navigation toolbar and from sig_rx
        """
        if self.isVisible():
            if 'data_changed' in dict_sig or self.needs_calc:
                # filter design has changed, filtered data is outdated
                self._refilter()
                self.needs_calc = False
            elif 'view_changed' in dict_sig\
                    or ('mpl_toolbar' in dict_sig and dict_sig['mpl_toolbar'] == 'home'):
                self.draw()
            elif 'mpl_toolbar' in dict_sig and dict_sig['mpl_toolbar'] == 'ui_level':
                self.frm_controls.setVisible(self.mplwidget.mpl_toolbar.a_ui_level == 0)
        else:
            if 'data_changed' in dict_sig or 'view_changed' in dict_sig:
                self.needs_calc = True

    # ------------------------------------------------------------------------------
    def load_data(self):
        """
        Select a file via file dialog and load it using pyfda's import routines.
        """
        file_name, file_type = io.select_file(
            self, title=self.tr("Load data to be filtered"), mode="r",
            file_types=FILE_TYPES)
        if file_name is None:  # operation cancelled
            return
        self.load_file(file_name, file_type)

    # ------------------------------------------------------------------------------
    def load_file(self, file_name: str, file_type: str) -> bool:
        """
        Load data from `file_name`, convert it to a 2D float array with one
        column per channel and populate the column combo boxes.
        """
        data = io.file2array(file_name, file_type, as_str=True)
        if data is None:
            logger.error("Couldn't load data from '%s'.", file_name)
            return False
        data = np.atleast_1d(np.squeeze(data))
        if data.ndim == 1:
            data = data[:, np.newaxis]
        elif data.ndim != 2:
            logger.error("Unsuitable data shape %s.", np.shape(data))
            return False
        if data.dtype.kind in {'U', 'S', 'O'}:
            data = self._to_float(data)
            if data is None:
                return False
            # drop header row(s) not removed by the csv import (all cells non-numeric)
            while len(data) > 1 and np.all(np.isnan(data[0])):
                data = data[1:]

        self.data = data
        self.file_name = file_name
        self.col_names = self._read_col_names(file_name, file_type, data.shape[1])
        self.lbl_file.setText(os.path.basename(file_name))
        self.lbl_file.setToolTip(file_name)
        logger.info("Loaded %d samples x %d column(s) from '%s'.",
                    data.shape[0], data.shape[1], file_name)

        self.cmb_col.blockSignals(True)
        self.cmb_time.blockSignals(True)
        self.cmb_col.clear()
        self.cmb_col.addItems(self.col_names)
        self.cmb_time.clear()
        self.cmb_time.addItem("n / f_S")
        self.cmb_time.addItems(self.col_names)
        # When there are several columns and the first one is monotonically increasing,
        # assume it's a time column and preselect the second column for filtering
        if data.shape[1] > 1 and np.all(np.diff(data[:, 0].real) > 0):
            self.cmb_time.setCurrentIndex(1)
            self.cmb_col.setCurrentIndex(1)
        self.cmb_col.blockSignals(False)
        self.cmb_time.blockSignals(False)

        self.select_data()
        return True

    # ------------------------------------------------------------------------------
    @staticmethod
    def _to_float(data: np.ndarray) -> np.ndarray | None:
        """
        Convert array of str to float (or complex), treating empty or non-numeric
        cells as NaN. Decimal commas are accepted when no dot is present.
        """
        def conv(s):
            s = str(s).strip()
            if ',' in s and '.' not in s:
                s = s.replace(',', '.')
            try:
                return complex(s.replace('i', 'j')) if 'j' in s or 'i' in s else float(s)
            except ValueError:
                return np.nan
        try:
            out = np.vectorize(conv, otypes=[complex])(data)
        except (TypeError, ValueError) as e:
            logger.error("Couldn't convert data to numbers:\n%s", e)
            return None
        if np.all(np.isnan(out)):
            logger.error("Data doesn't contain any numeric values.")
            return None
        if not np.any(np.iscomplex(out)):
            out = out.real
        return out

    # ------------------------------------------------------------------------------
    @staticmethod
    def _read_col_names(file_name: str, file_type: str, n_cols: int) -> list:
        """
        Try to read column names from the header line of a csv / txt file,
        otherwise return generic names.
        """
        names = [f"Col {i + 1}" for i in range(n_cols)]
        if file_type not in {'csv', 'txt'}:
            return names
        try:
            with open(file_name, 'r', newline=None) as f:
                line = f.readline()
        except (IOError, UnicodeDecodeError):
            return names
        if not line.strip() or re.fullmatch(r'[eEjJ()0-9,;:|\.\+\-\s]+', line.strip()):
            return names  # empty or no header (only numeric characters)
        delim = params['CSV']['delimiter']
        if delim.lower() == 'auto':
            try:
                delim = csv.Sniffer().sniff(line, delimiters=['\t', ';', ',', '|', ' ']).delimiter
            except csv.Error:
                return names
        header = [h.strip() for h in next(csv.reader([line], delimiter=delim))]
        if len(header) == n_cols:
            return [h if h else names[i] for i, h in enumerate(header)]
        return names

    # ------------------------------------------------------------------------------
    def select_data(self):
        """
        Pick the selected column as data to be filtered, invalidate filtered data.
        """
        if self.data is None or self.cmb_col.currentIndex() < 0:
            return
        self.x = self.data[:, self.cmb_col.currentIndex()]
        self.y = None
        self.but_filter.setEnabled(True)
        self.but_export.setEnabled(False)
        self.draw()

    # ------------------------------------------------------------------------------
    def _refilter(self):
        """ Update filtered data when it has been calculated before """
        if self.y is not None:
            self.filter_data()
        else:
            self.draw()

    # ------------------------------------------------------------------------------
    def filter_data(self):
        """
        Filter the selected data column with the current filter design, using
        second-order sections when available for better numerical stability.
        """
        if self.x is None:
            return
        x = self.x
        nan_mask = np.isnan(x)
        if np.any(nan_mask):
            logger.warning("Replacing %d non-numeric value(s) by zero.", np.sum(nan_mask))
            x = np.where(nan_mask, 0, x)

        sos = np.asarray(fb_get('sos'))
        bb = np.asarray(fb_get('ba', 0))
        aa = np.asarray(fb_get('ba', 1))
        # coefficients are stored as complex, use real part when possible
        if not np.any(np.iscomplex(bb)) and not np.any(np.iscomplex(aa)):
            bb, aa = bb.real, aa.real
        # remove trailing zeros from the denominator
        aa = np.trim_zeros(aa, 'b') if np.any(aa) else np.array([1.])

        zero_phase = self.chk_zero_phase.isChecked()
        try:
            if sos.ndim == 2 and sos.shape[0] > 0 and sos.shape[1] == 6\
                    and not np.any(np.iscomplex(sos)):
                sos = sos.real
                if zero_phase:
                    self.y = sig.sosfiltfilt(sos, x)
                else:
                    self.y = sig.sosfilt(sos, x)
            else:
                if zero_phase:
                    self.y = sig.filtfilt(bb, aa, x)
                else:
                    self.y = sig.lfilter(bb, aa, x)
        except ValueError as e:
            logger.error("Filtering failed:\n%s", e)
            self.y = None
        self.but_export.setEnabled(self.y is not None)
        self.draw()

    # ------------------------------------------------------------------------------
    def _time_axis(self):
        """ Return time axis and its label """
        idx = self.cmb_time.currentIndex()
        if self.data is not None and idx > 0:
            return self.data[:, idx - 1].real, self.col_names[idx - 1]
        f_s = fb_get('f_s')
        n = np.arange(len(self.x))
        return n / f_s, fb_get('plt_t_label')

    # ------------------------------------------------------------------------------
    def draw(self):
        """ (Re-)draw the figure """
        self.mplwidget.fig.clf()
        if self.chk_spectrum.isChecked():
            self.ax_t = self.mplwidget.fig.add_subplot(2, 1, 1)
            self.ax_f = self.mplwidget.fig.add_subplot(2, 1, 2)
        else:
            self.ax_t = self.mplwidget.fig.add_subplot(1, 1, 1)
            self.ax_f = None

        if self.x is None:
            self.ax_t.text(0.5, 0.5, self.tr("Load data with 'Load data ...', then press "
                                             "'Filter data'"),
                           ha='center', va='center', transform=self.ax_t.transAxes)
            self.ax_t.set_xticks([])
            self.ax_t.set_yticks([])
            self.mplwidget.redraw()
            return

        t, t_label = self._time_axis()
        name = self.col_names[self.cmb_col.currentIndex()]
        cplx = np.iscomplexobj(self.x) or (self.y is not None and np.iscomplexobj(self.y))
        if cplx:
            logger.info("Complex data, only the real part is plotted.")
        self.ax_t.plot(t, self.x.real, label=self.tr("Original") + f" ({name})",
                       alpha=0.6 if self.y is not None else 1)
        if self.y is not None:
            self.ax_t.plot(t, self.y.real, label=self.tr("Filtered"))
        self.ax_t.set_xlabel(t_label)
        self.ax_t.set_ylabel(name)
        self.ax_t.legend(loc='best')
        self.ax_t.grid(True)

        if self.ax_f is not None:
            f_s = fb_get('f_s')
            N = len(self.x)
            f = np.fft.rfftfreq(N, d=1/f_s)
            for d, lbl in ((self.x, self.tr("Original")), (self.y, self.tr("Filtered"))):
                if d is None:
                    continue
                d = np.nan_to_num(d.real)
                X = np.abs(np.fft.rfft(d)) / N
                X[1:] *= 2  # single-sided spectrum
                self.ax_f.plot(f, 20 * np.log10(np.maximum(X, 1e-12)), label=lbl)
            self.ax_f.set_xlabel(fb_get('plt_f_label'))
            self.ax_f.set_ylabel(r"$|X(f)|$ in dB")
            self.ax_f.legend(loc='best')
            self.ax_f.grid(True)

        self.mplwidget.fig.tight_layout()
        self.mplwidget.redraw()

    # ------------------------------------------------------------------------------
    def export_data(self):
        """ Save time axis, original and filtered data as a CSV file """
        if self.y is None:
            return
        file_name, file_type = io.select_file(
            self, title=self.tr("Export filtered data"), mode="w", file_types=('csv',))
        if file_name is None:
            return
        t, t_label = self._time_axis()
        name = self.col_names[self.cmb_col.currentIndex()]
        try:
            with open(file_name, 'w', newline='') as f:
                wr = csv.writer(f, delimiter=',')
                wr.writerow(["t", name, name + "_filtered"])
                for row in zip(t, self.x.real, self.y.real):
                    wr.writerow([repr(float(v)) for v in row])
            logger.info("Exported filtered data to '%s'.", file_name)
        except IOError as e:
            logger.error("Couldn't write '%s':\n%s", file_name, e)


# ------------------------------------------------------------------------------
def main():
    import sys
    from pyfda.libs.compat import QApplication

    app = QApplication(sys.argv)
    mainw = PlotDataFilt()
    app.setActiveWindow(mainw)
    mainw.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
    # module test using python -m pyfda.plot_widgets.plot_data_filt
