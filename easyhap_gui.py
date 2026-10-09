"""Single Windows EXE: no arguments launches the GUI; arguments use the CLI."""

from __future__ import annotations

import multiprocessing
import sys


def main() -> None:
    multiprocessing.freeze_support()
    if len(sys.argv) > 1 and sys.argv[1] not in {"--gui", "gui"}:
        from easyhap.cli import main as cli_main

        raise SystemExit(cli_main())
    from easyhap.gui import main as gui_main

    gui_main()


if __name__ == "__main__":
    main()
