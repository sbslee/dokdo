from .collapse import collapse
from .make_manifest import make_manifest
from .add_metadata import add_metadata
from .summarize import summarize
from .prepare_lefse import prepare_lefse

commands = {
    "collapse": collapse,
    "make-manifest": make_manifest,
    "add-metadata": add_metadata,
    "summarize": summarize,
    "prepare-lefse": prepare_lefse,
}
