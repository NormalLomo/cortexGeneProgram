"""Process-stable BLAKE2b seeds for spatial null calculations."""
from __future__ import annotations

import hashlib


def stable_seed(*parts: object) -> int:
    payload = "\x1f".join(map(str, parts)).encode("utf-8")
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "big")
