"""
utils.py

Subprocess helpers and logging setup.
"""

from __future__ import annotations

import logging
import shlex
import subprocess
from pathlib import Path
from typing import List, Optional

LOG = logging.getLogger("lai-pipeline")


def setup_logging(level: str) -> None:
    """Configure root logger format and level."""
    raise NotImplementedError


def shjoin(cmd: List[str]) -> str:
    """Return a shell-quoted string representation of a command list."""
    raise NotImplementedError


def run(
    cmd: List[str],
    *,
    cwd: Optional[Path] = None,
    check: bool = True,
    capture_stdout: bool = False,
    capture_stderr: bool = True,
    text: bool = True,
) -> subprocess.CompletedProcess:
    """Run a command, log it, and return the CompletedProcess."""
    raise NotImplementedError


def popen_lines(cmd: List[str], *, cwd: Optional[Path] = None) -> subprocess.Popen:
    """Open a subprocess and return the Popen object for line-by-line streaming."""
    raise NotImplementedError


def count_stream_lines(proc: subprocess.Popen, *, label: str) -> int:
    """Count lines from a Popen stdout, logging progress every 1M lines."""
    raise NotImplementedError