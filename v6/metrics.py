"""Simple in-memory counters for lightweight app metrics."""

from threading import Lock

_counters = {}
_lock = Lock()


def init_metrics() -> None:
    """Initialize required counters once at startup."""
    with _lock:
        _counters.setdefault("visits", 0)
        _counters.setdefault("gene_queries", 0)


def increment(key: str, amount: int = 1) -> None:
    """Increase a named counter by `amount`."""
    with _lock:
        _counters[key] = _counters.get(key, 0) + amount


def get_value(key: str) -> int:
    """Return current counter value (0 if missing)."""
    with _lock:
        return _counters.get(key, 0)
