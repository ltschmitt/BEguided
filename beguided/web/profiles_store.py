"""In-memory store and validator for user-uploaded base editor profiles TSV.

The web app lets users upload a custom base editor profiles TSV that
replaces the operator-configured file for the lifetime of the server.
Uploads are validated against a fixed schema and rejected if they fail
any check. The store itself is a thread-safe singleton; the scheduler
is single-worker but uploads come from the async event loop.
"""

from __future__ import annotations

import base64
import io
import threading
from typing import List, Optional

import pandas as pd


REQUIRED_PROFILE_COLUMNS = (
    "name",
    "editing",
    "general_edit_start",
    "general_edit_end",
    "high_edit_start",
    "high_edit_end",
    "pam",
)
MAX_PROFILE_ROWS = 10
MAX_PROFILE_PAM_LEN = 6


class ProfilesUploadError(ValueError):
    """Raised when an uploaded profiles TSV fails validation."""

    def __init__(self, message: str, status_code: int = 400) -> None:
        super().__init__(message)
        self.status_code = status_code


def validate_uploaded_profiles(content: bytes, max_bytes: int) -> pd.DataFrame:
    """Validate an uploaded profiles TSV and return the parsed DataFrame.

    Checks (in order): file size, parse-ability, required columns, row
    count, PAM length. Extras columns are dropped; only the required
    columns are returned.

    Raises:
        ProfilesUploadError: with a user-facing message and HTTP status.
    """
    if len(content) > max_bytes:
        size_mb = len(content) / (1024 * 1024)
        max_mb = max_bytes / (1024 * 1024)
        raise ProfilesUploadError(
            f"File too large: {size_mb:.2f} MB (max {max_mb:g} MB)"
        )

    try:
        df = pd.read_csv(io.BytesIO(content), sep="\t")
    except Exception as exc:
        raise ProfilesUploadError(f"Failed to parse TSV: {exc}")

    missing = set(REQUIRED_PROFILE_COLUMNS) - set(df.columns)
    if missing:
        raise ProfilesUploadError(
            f"Missing required columns: {', '.join(sorted(missing))}. "
            f"Found: {', '.join(df.columns.tolist())}"
        )

    if len(df) > MAX_PROFILE_ROWS:
        raise ProfilesUploadError(
            f"Too many profiles: {len(df)} (max {MAX_PROFILE_ROWS})"
        )

    pam_lengths = df["pam"].astype(str).str.len()
    bad_mask = pam_lengths > MAX_PROFILE_PAM_LEN
    if bad_mask.any():
        offenders: List[str] = []
        for i in df.index[bad_mask]:
            offenders.append(
                f"'{df.at[i, 'name']}' (PAM '{df.at[i, 'pam']}', "
                f"{int(pam_lengths.at[i])} bases)"
            )
        raise ProfilesUploadError(
            f"PAM too long (max {MAX_PROFILE_PAM_LEN} bases): "
            f"{', '.join(offenders)}"
        )

    return df[list(REQUIRED_PROFILE_COLUMNS)].copy()


class CustomProfilesStore:
    """Singleton holding one uploaded base editor profiles TSV in memory."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._bytes: Optional[bytes] = None
        self._names: List[str] = []

    def set(self, content: bytes) -> None:
        with self._lock:
            self._bytes = content
            try:
                df = pd.read_csv(io.BytesIO(content), sep="\t")
                self._names = sorted(df["name"].astype(str).tolist())
            except Exception:
                self._names = []

    def clear(self) -> None:
        with self._lock:
            self._bytes = None
            self._names = []

    def is_loaded(self) -> bool:
        with self._lock:
            return self._bytes is not None

    def get_bytes(self) -> Optional[bytes]:
        with self._lock:
            return self._bytes

    def get_b64(self) -> Optional[str]:
        with self._lock:
            if self._bytes is None:
                return None
            return base64.b64encode(self._bytes).decode("ascii")

    def names(self) -> List[str]:
        with self._lock:
            return list(self._names)


_store = CustomProfilesStore()


def get_store() -> CustomProfilesStore:
    """Return the module-level custom profiles store singleton."""
    return _store
