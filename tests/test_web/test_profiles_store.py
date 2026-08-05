"""Unit tests for beguided.web.profiles_store."""

from __future__ import annotations

import base64

import pytest


def test_validate_uploaded_profiles_accepts_minimal_columns():
    from beguided.web.profiles_store import validate_uploaded_profiles

    body = (
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\n"
        "ABE8e\tABE\t2\t13\t3\t8\tNGG\n"
    )
    df = validate_uploaded_profiles(body.encode(), ".tsv", 1024 * 1024)
    assert list(df.columns) == [
        "name",
        "editing",
        "general_edit_start",
        "general_edit_end",
        "high_edit_start",
        "high_edit_end",
        "pam",
    ]
    assert len(df) == 1
    assert df.iloc[0]["name"] == "ABE8e"


def test_validate_uploaded_profiles_drops_extra_columns():
    from beguided.web.profiles_store import validate_uploaded_profiles

    body = (
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\tcas\tsource\textra\n"
        "ABE8e\tABE\t2\t13\t3\t8\tNGG\tSpCas9\tNature2022\tdropme\n"
    )
    df = validate_uploaded_profiles(body.encode(), ".tsv", 1024 * 1024)
    assert "cas" not in df.columns
    assert "source" not in df.columns
    assert "extra" not in df.columns


def test_validate_uploaded_profiles_size_limit():
    from beguided.web.profiles_store import (
        ProfilesUploadError,
        validate_uploaded_profiles,
    )

    body = b"x" * 2048
    with pytest.raises(ProfilesUploadError) as exc:
        validate_uploaded_profiles(body, ".tsv", 1024)
    assert "too large" in str(exc.value).lower()


def test_validate_uploaded_profiles_missing_columns():
    from beguided.web.profiles_store import (
        ProfilesUploadError,
        validate_uploaded_profiles,
    )

    body = "name\tediting\nABE8e\tABE\n"
    with pytest.raises(ProfilesUploadError) as exc:
        validate_uploaded_profiles(body.encode(), ".tsv", 1024 * 1024)
    assert "Missing required columns" in str(exc.value)


def test_validate_uploaded_profiles_row_limit():
    from beguided.web.profiles_store import (
        MAX_PROFILE_ROWS,
        ProfilesUploadError,
        validate_uploaded_profiles,
    )

    header = (
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\n"
    )
    body = header + "".join(
        f"BE{i}\tABE\t2\t13\t3\t8\tNGG\n" for i in range(MAX_PROFILE_ROWS + 1)
    )
    with pytest.raises(ProfilesUploadError) as exc:
        validate_uploaded_profiles(body.encode(), ".tsv", 1024 * 1024)
    assert "Too many profiles" in str(exc.value)


def test_validate_uploaded_profiles_pam_length_limit():
    from beguided.web.profiles_store import (
        MAX_PROFILE_PAM_LEN,
        ProfilesUploadError,
        validate_uploaded_profiles,
    )

    header = (
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\n"
    )
    too_long = header + f"ABE8e\tABE\t2\t13\t3\t8\t{'N' * (MAX_PROFILE_PAM_LEN + 1)}\n"
    with pytest.raises(ProfilesUploadError) as exc:
        validate_uploaded_profiles(too_long.encode(), ".tsv", 1024 * 1024)
    assert "PAM too long" in str(exc.value)

    exactly_max = header + f"ABE8e\tABE\t2\t13\t3\t8\t{'N' * MAX_PROFILE_PAM_LEN}\n"
    df = validate_uploaded_profiles(exactly_max.encode(), ".tsv", 1024 * 1024)
    assert len(df) == 1


def test_store_roundtrip():
    from beguided.web.profiles_store import CustomProfilesStore

    s = CustomProfilesStore()
    assert s.is_loaded() is False
    body = (
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\n"
        "ABE8e\tABE\t2\t13\t3\t8\tNGG\n"
        "CBE6d\tCBE\t2\t10\t6\t6\tNGG\n"
    )
    s.set(body.encode())
    assert s.is_loaded() is True
    assert s.names() == ["ABE8e", "CBE6d"]
    assert s.get_bytes() == body.encode()
    assert base64.b64decode(s.get_b64()) == body.encode()

    s.clear()
    assert s.is_loaded() is False
    assert s.names() == []
    assert s.get_bytes() is None
    assert s.get_b64() is None
