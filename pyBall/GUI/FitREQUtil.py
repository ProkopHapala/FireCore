import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

EV_TO_KCAL = 23.060549

class BlitManager:
    """Manages blitting for fast Matplotlib updates in Qt/VisPy applications."""
    def __init__(self, canvas, animated_artists=()):
        self.canvas = canvas
        self._bg = None
        self._artists = []
        for a in animated_artists:
            self.add_artist(a)
        self.cid = canvas.mpl_connect("draw_event", self.on_draw)

    def on_draw(self, event):
        if event is not None and event.canvas != self.canvas:
            return
        self._bg = self.canvas.copy_from_bbox(self.canvas.figure.bbox)
        self._draw_animated()

    def add_artist(self, art):
        art.set_animated(True)
        self._artists.append(art)

    def _draw_animated(self):
        fig = self.canvas.figure
        for a in self._artists:
            fig.draw_artist(a)

    def update(self):
        if self._bg is None:
            self.on_draw(None)
        else:
            self.canvas.restore_region(self._bg)
            self._draw_animated()
            self.canvas.blit(self.canvas.figure.bbox)
        self.canvas.flush_events()

def plot_stacked_1d(ax, x, total, components, labels, colors, ref=None, title=None, xlabel=None, ylabel=None):
    """
    Plot stacked area energy decomposition.
    components: list of arrays, each of same length as x.
    """
    ax.clear()
    if ref is not None:
        ax.plot(x, ref, 'k--', label='Ref', lw=1.5)
    
    ax.plot(x, total, 'r-', label='Total', lw=2)
    
    current_bottom = np.zeros_like(x)
    for comp, label, color in zip(components, labels, colors):
        ax.fill_between(x, current_bottom, current_bottom + comp, color=color, alpha=0.3, label=label)
        current_bottom += comp
        
    if title: ax.set_title(title)
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    ax.legend(loc='best', fontsize='small')
    ax.grid(True, ls=':')

def extract_cuts(V, row, col):
    """Extract radial and angular cuts from 2D grid V."""
    ny, nx = V.shape
    row = max(0, min(row, ny - 1))
    col = max(0, min(col, nx - 1))
    r_cut = V[row, :]
    a_cut = V[:, col]
    return r_cut, a_cut

def reshape_to_grid_proper(vals, r, a, rows):
    """
    Pad rows of a scan to a rectangular grid.
    This version is compatible with split_scan_imshow_new.py logic.
    """
    ny = len(rows)
    nx = max(e - s for s, e in rows)
    V = np.full((ny, nx), np.nan)
    R = np.full((ny, nx), np.nan)
    A = np.full((ny,), np.nan)
    rv = np.full((nx,), np.nan)
    for iy, (s, e) in enumerate(rows):
        n = e - s
        V[iy, :n] = vals[s:e]
        R[iy, :n] = r[s:e]
        try:
            if n > 0:
                irel = int(np.nanargmax(r[s:e]))
                A[iy] = a[s + irel]
            else:
                A[iy] = np.nan
        except:
            A[iy] = np.nanmean(a[s:e])
            
    # Always extract rv from the first row BEFORE sorting
    rv = R[0, :].copy()
    
    # Sort by angle to make the grid continuous
    order = np.argsort(A)
    A = A[order]
    V = V[order, :]
    R = R[order, :]
    sorted_rows = [rows[i] for i in order]
    
    return V, R, A, rv, sorted_rows
