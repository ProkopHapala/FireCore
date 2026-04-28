"""Generic JSON-Lines results database helper.

Each line in the file is expected to be a *single* JSON object whose top-level
keys are unique solution identifiers (strings). The corresponding value is an
arbitrary JSON-serialisable dict (e.g. {"rmse": 1e-3, ...}).

The class focuses only on IO, uniqueness, and append-only behaviour.  Any
problem-specific logic (parsing the key, computing metrics, etc.) can be passed
in via call-backs or handled externally.
"""
from __future__ import annotations

import json, pathlib, io, typing, contextlib
from typing import Dict, Iterator, Callable, Any

JsonDict = Dict[str, Any]

class ResultsDB:
    """Append-only JSON-Lines database (one JSON object per line)."""

    def __init__(self, path: str | pathlib.Path):
        self.path = pathlib.Path(path)
        self._loaded: bool = False
        self._records: Dict[str, JsonDict] = {}

    # ---------------------------------------------------------------------
    # Internal helpers
    # ---------------------------------------------------------------------
    def _load(self):
        if self._loaded:
            return
        if not self.path.is_file():
            self._loaded = True
            return
        with self.path.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    obj = json.loads(line)
                except json.JSONDecodeError:
                    continue  # silently skip malformed line
                for k, v in obj.items():
                    self._records[k] = v
        self._loaded = True

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def __contains__(self, key: str) -> bool:  # membership test
        self._load()
        return key in self._records

    def get(self, key: str, default=None):
        self._load()
        return self._records.get(key, default)

    def keys(self):
        self._load()
        return self._records.keys()

    def items(self):
        self._load()
        return self._records.items()

    def append(self, key: str, data: JsonDict):
        """Append a new record *if* the key is not already present."""
        self._load()
        if key in self._records:
            return False  # duplicate
        # write atomically: append then flush
        with self.path.open("a", encoding="utf-8") as f:
            json.dump({key: data}, f, separators=(",", ":"))
            f.write("\n")
        self._records[key] = data
        return True

    # Iterable interface -------------------------------------------------
    def __iter__(self) -> Iterator[tuple[str, JsonDict]]:
        self._load()
        return iter(self._records.items())

    # Convenience --------------------------------------------------------
    def filter(self, predicate: Callable[[str, JsonDict], bool]) -> Dict[str, JsonDict]:
        """Return new dict with records matching *predicate*."""
        return {k: v for k, v in self if predicate(k, v)}

    # ------------------------------------------------------------------
    # Pareto front utility (2-D minimisation)
    # ------------------------------------------------------------------
    def pareto_front(
        self,
        metric_fn: Callable[[str, JsonDict], tuple[float, float]],
        *,
        unique_x: bool = False,
    ) -> list[str]:
        """Return list of keys on the 2-D Pareto front.

        Parameters
        ----------
        metric_fn : callable
            Maps (key, record) → (x, y) where both axes are to be minimised.
        unique_x : bool, optional
            If True, keep only the record with minimal *y* for each distinct *x*
            **before** computing the Pareto front (e.g. ensure one solution per
            nBasis bucket).
        """
        pts = []  # (key, x, y)
        for k, rec in self.items():
            try:
                x, y = metric_fn(k, rec)
            except Exception:
                continue
            pts.append((k, x, y))

        if unique_x:
            best_for_x = {}
            for k, x, y in pts:
                if (x not in best_for_x) or (y < best_for_x[x][2]):
                    best_for_x[x] = (k, x, y)
            pts = list(best_for_x.values())

        # sort by x ascending
        pts.sort(key=lambda t: t[1])

        pareto_keys: list[str] = []
        best_y = float("inf")
        for k, _x, y in pts:
            if y < best_y:
                pareto_keys.append(k)
                best_y = y
        return pareto_keys
