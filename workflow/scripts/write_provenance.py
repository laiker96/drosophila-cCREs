"""Persist the fully resolved workflow configuration in the result namespace."""

import json
import os
from pathlib import Path
import tempfile


output = Path(str(snakemake.output.config))
log = Path(str(snakemake.log[0]))
output.parent.mkdir(parents=True, exist_ok=True)
log.parent.mkdir(parents=True, exist_ok=True)
resolved = json.loads(str(snakemake.params.config))
content = json.dumps(resolved, indent=2, sort_keys=True) + "\n"
temporary_name = None
try:
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=f".{output.name}.",
        dir=output.parent,
        delete=False,
    ) as handle:
        temporary_name = handle.name
        handle.write(content)
    os.replace(temporary_name, output)
finally:
    if temporary_name and os.path.exists(temporary_name):
        os.unlink(temporary_name)
log.write_text("Wrote fully resolved workflow configuration\n", encoding="utf-8")
