"""In-memory store and validator for user-uploaded base editor profiles.

The web app lets users upload a custom base editor profiles spreadsheet
that replaces the operator-configured file for the lifetime of the
server. Uploads are validated against a fixed schema and rejected if
they fail any check. The store itself is a thread-safe singleton; the
scheduler is single-worker but uploads come from the async event loop.

Both XLSX (``application/vnd.openxmlformats-officedocument.spreadsheetml.sheet``)
and TSV are accepted. CSV is also accepted defensively. The downstream
pipeline only knows how to read TSV, so the store normalises whatever
format was uploaded back to TSV bytes; the route hands the parser a
DataFrame directly and stores its TSV serialisation. This keeps the
disk format invariant (always TSV) and the audit label in the
scheduler (``tsv_b64``) honest.
"""

from __future__ import annotations

import base64
import io
import threading
from pathlib import Path
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

# Extensions accepted by ``validate_uploaded_profiles``. ``.xls`` is listed
# for completeness; we cannot actually parse it without the optional
# ``xlrd`` dependency, and uploads in that format are rejected with a
# clear "re-save as .xlsx" hint at the route layer.
SUPPORTED_PROFILE_SUFFIXES = (".xlsx", ".xls", ".tsv", ".csv")

# Magic-byte prefixes used to reject uploads whose extension lies about
# their contents (mirrors ``_sniff_upload_kind`` in ``app.py`` for the
# targets route).
_XLSX_MAGIC = b"PK\x03\x04"
_XLS_MAGIC = b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1"


class ProfilesUploadError(ValueError):
    """Raised when an uploaded profiles spreadsheet fails validation."""

    def __init__(self, message: str, status_code: int = 400) -> None:
        super().__init__(message)
        self.status_code = status_code


def _sniff_profile_bytes(content: bytes, suffix: str) -> str:
    """Return the inferred profile kind (``"xlsx"``/``"xls"``/``"tsv"``/``"csv"``/``"binary"``).

    Anything that contains a NUL byte in the first 1 KiB is treated as
    binary unless the extension is an Office format whose magic matches.
    """
    head = content[:1024]
    if not head:
        return "binary"
    sfx = suffix.lower()
    if sfx in (".xlsx", ".xls"):
        if sfx == ".xlsx" and head.startswith(_XLSX_MAGIC):
            return "xlsx"
        if sfx == ".xls" and head.startswith(_XLS_MAGIC):
            return "xls"
        return "binary"
    try:
        head_text = head.decode("utf-8")
    except UnicodeDecodeError:
        return "binary"
    if "\x00" in head_text:
        return "binary"
    if sfx == ".csv":
        return "csv" if "," in head_text else "tsv"
    if "\t" in head_text:
        return "tsv"
    if "," in head_text:
        return "csv"
    return "tsv"


def _read_profile_dataframe(content: bytes, suffix: str) -> pd.DataFrame:
    """Parse a profile spreadsheet payload into a DataFrame.

    The reader is chosen from ``suffix`` first, then validated against
    the magic bytes via ``_sniff_profile_bytes``. ``.xls`` is rejected
    explicitly because parsing it requires the optional ``xlrd`` package
    and the user is much better served by re-saving the workbook as
    ``.xlsx`` (which is what the template is).
    """
    sfx = suffix.lower()
    if sfx not in SUPPORTED_PROFILE_SUFFIXES:
        raise ProfilesUploadError(
            f"Unsupported file type: '{sfx or '<no extension>'}'. "
            f"Allowed: {', '.join(SUPPORTED_PROFILE_SUFFIXES)}."
        )

    kind = _sniff_profile_bytes(content, sfx)
    if kind == "binary":
        raise ProfilesUploadError(
            f"Uploaded file contents do not match the {sfx} extension. "
            "Allowed: .xlsx, .tsv, .csv."
        )

    try:
        if kind == "xlsx":
            return pd.read_excel(io.BytesIO(content), engine="openpyxl")
        if kind == "xls":
            raise ProfilesUploadError(
                "Legacy .xls format is not supported. Please re-save the "
                "workbook as .xlsx and upload that file."
            )
        if kind == "csv":
            return pd.read_csv(io.BytesIO(content))
        return pd.read_csv(io.BytesIO(content), sep="\t")
    except ProfilesUploadError:
        raise
    except Exception as exc:
        raise ProfilesUploadError(f"Failed to parse file: {exc}")


def validate_uploaded_profiles(
    content: bytes,
    suffix: str,
    max_bytes: int,
    default_profiles_file: Optional[str] = None,
) -> pd.DataFrame:
    """Validate an uploaded profiles spreadsheet and return its DataFrame.

    Checks (in order): file size, parse-ability, required columns, row
    count, PAM length, name collisions with the default profiles file.
    Extra columns are dropped; only the required columns are returned.

    Args:
        content: Raw bytes from the upload (any of ``.xlsx``/``.tsv``/``.csv``).
        suffix: Lower-cased filename suffix (e.g. ``".xlsx"``) used to pick
            the parser and to flag mismatches between the extension and
            the magic bytes.
        max_bytes: Maximum file size in bytes.
        default_profiles_file: Optional path to the operator-configured
            default profiles file. When provided, any uploaded name that
            already exists in this file triggers a collision error.

    Raises:
        ProfilesUploadError: with a user-facing message and HTTP status.
    """
    if len(content) > max_bytes:
        size_mb = len(content) / (1024 * 1024)
        max_mb = max_bytes / (1024 * 1024)
        raise ProfilesUploadError(
            f"File too large: {size_mb:.2f} MB (max {max_mb:g} MB)"
        )

    df = _read_profile_dataframe(content, suffix)

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

    if default_profiles_file and Path(default_profiles_file).exists():
        try:
            default_df = pd.read_csv(default_profiles_file, sep="\t")
            default_names = set(default_df["name"].astype(str).tolist())
            uploaded_names = df["name"].astype(str).tolist()
            collisions = sorted(n for n in uploaded_names if n in default_names)
            if collisions:
                raise ProfilesUploadError(
                    f"Name collision with default profile(s): "
                    f"{', '.join(collisions)}. Custom profiles must use unique "
                    f"names not present in the default profiles file."
                )
        except ProfilesUploadError:
            raise
        except Exception:
            # If we can't read the default file for collision check, skip it
            # rather than failing the upload.
            pass

    return df[list(REQUIRED_PROFILE_COLUMNS)].copy()


def _dataframe_to_tsv_bytes(df: pd.DataFrame) -> bytes:
    """Serialise a validated profiles DataFrame back to TSV bytes.

    The store keeps the canonical TSV representation so that
    ``get_bytes()`` always returns the same format, regardless of what
    the user uploaded. This also means ``load_base_editor_profiles`` and
    the workdir writer in ``tasks.py`` keep reading TSV with no change.
    """
    return df.to_csv(sep="\t", index=False).encode("utf-8")


class CustomProfilesStore:
    """Singleton holding one uploaded base editor profiles spreadsheet in memory.

    Whatever format the user uploads (XLSX or TSV), the bytes returned by
    :meth:`get_bytes` are always TSV. Use :meth:`set_dataframe` from the
    upload route once ``validate_uploaded_profiles`` has returned a clean
    DataFrame.
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._bytes: Optional[bytes] = None
        self._names: List[str] = []

    def set_dataframe(self, df: pd.DataFrame) -> None:
        """Store the validated DataFrame as canonical TSV bytes."""
        tsv_bytes = _dataframe_to_tsv_bytes(df)
        with self._lock:
            self._bytes = tsv_bytes
            self._names = sorted(df["name"].astype(str).tolist())

    def set(self, content: bytes) -> None:
        """Backwards-compatible entry point used by tests.

        Treats ``content`` as TSV bytes (the format the store has always
        held). New callers should prefer :meth:`set_dataframe` after
        validating an arbitrary spreadsheet format.
        """
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
