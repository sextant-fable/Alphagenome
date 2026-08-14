#!/usr/bin/env python3
"""Flag claim-boundary, terminology and common templated-prose risks."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


FLAG_RULES = (
    ("style", r"\b(?:leverage|delve(?: into)?|underscore|unveil|pivotal|transformative|groundbreaking|revolutionary)\b", "Prefer a concrete technical verb if the word adds no technical meaning."),
    ("style", r"\b(?:comprehensive|powerful|robust|novel|notably)\b", "Check that a nearby result defines the claimed property."),
    ("style", r"\b(?:it is worth noting|first and foremost|in conclusion)\b", "Replace formulaic signposting with the actual logical relation or delete it."),
    ("claim", r"\bfirst(?: ever)?\s+(?:model|method|approach|system|tool|demonstration|study|report|framework)\b|\b(?:state[- ]of[- ]the[- ]art|sota|universal|always|never)\b", "Verify the scope against a frozen comparison or remove the universal claim."),
    ("claim", r"\bDPY-27\b.{0,80}\b(?:biological validation|dosage[- ]compensation recovery)\b", "P10 is a negative internal diagnostic and cannot support this framing."),
    ("term", r"(?<!size-matched )\bModel C\b", "Use `size-matched Model C` when parameter matching is relevant."),
    ("term", r"\b(?:fine[- ]tuning|full fine[- ]tuning)\b", "Use `LoRA adaptation` unless all model parameters were optimized."),
)

BLOCKER_RULES = (
    ("placeholder", r"\bEXTERNAL_[A-Z0-9_]+\b|\bEXT(?:_|\\_)[A-Z0-9_\\]+", "Freeze the corresponding external artifact before drafting final prose."),
    ("boundary", r"\b(?:DPY-27 response recovery|DPY-27 biological validation)\b", "P10 is a negative internal diagnostic and cannot support this framing."),
)


def find_issues(text: str) -> list[dict[str, object]]:
    issues: list[dict[str, object]] = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        for category, pattern, guidance in BLOCKER_RULES:
            for match in re.finditer(pattern, line, flags=re.IGNORECASE):
                issues.append(
                    {
                        "severity": "blocker",
                        "category": category,
                        "line": line_number,
                        "match": match.group(0),
                        "guidance": guidance,
                    }
                )
        for category, pattern, guidance in FLAG_RULES:
            for match in re.finditer(pattern, line, flags=re.IGNORECASE):
                issues.append(
                    {
                        "severity": "review",
                        "category": category,
                        "line": line_number,
                        "match": match.group(0),
                        "guidance": guidance,
                    }
                )
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path, help="UTF-8 manuscript text or Markdown file")
    parser.add_argument("--strict", action="store_true", help="Exit non-zero when blockers are found")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    text = args.path.read_text(encoding="utf-8")
    issues = find_issues(text)
    blockers = [issue for issue in issues if issue["severity"] == "blocker"]

    if args.json:
        print(json.dumps({"path": str(args.path), "issues": issues}, indent=2))
    elif issues:
        for issue in issues:
            print(
                f"{issue['severity'].upper()} line {issue['line']}: "
                f"{issue['category']} `{issue['match']}` - {issue['guidance']}"
            )
    else:
        print("PASS: no configured claim-boundary, terminology or style flags.")

    return 1 if args.strict and blockers else 0


if __name__ == "__main__":
    raise SystemExit(main())
