from .collapse import collapse
from .make_manifest import make_manifest
from .add_metadata import add_metadata
from .merge_metadata import merge_metadata
from .summarize import summarize

commands = {
    "collapse": collapse,
    "make-manifest": make_manifest,
    "add-metadata": add_metadata,
    "merge-metadata": merge_metadata,
    "summarize": summarize,
}
