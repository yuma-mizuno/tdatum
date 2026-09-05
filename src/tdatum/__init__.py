"""T-data and mutation loops for SageMath."""

try:
    import sage.all  # Initialize Sage before importing individual Sage modules.
except ModuleNotFoundError as exc:
    if exc.name == "sage":
        raise ImportError("tdatum requires SageMath; use sage -python or a SageMath kernel.") from exc
    raise

from .t_datum import MutationLoop, TDatum

__all__ = ["TDatum", "MutationLoop"]
__version__ = "0.1.0.dev0"
