#!/usr/bin/env python3
"""Extract the embedded GFA subgraph from a Pangene/gfa-server HTML page."""

from __future__ import annotations

import argparse
from html import unescape
from html.parser import HTMLParser
from pathlib import Path


class GfaTextareaParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.in_gfa = False
        self.parts: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "textarea" and dict(attrs).get("id") == "gfa-text":
            self.in_gfa = True

    def handle_endtag(self, tag: str) -> None:
        if tag == "textarea":
            self.in_gfa = False

    def handle_data(self, data: str) -> None:
        if self.in_gfa:
            self.parts.append(data)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract the readonly GFA textarea from a Pangene viewer HTML file."
    )
    parser.add_argument("--input", required=True, help="Pangene/gfa-server HTML file")
    parser.add_argument("--output", required=True, help="Output GFA file")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    html = Path(args.input).read_text(encoding="utf-8")
    parser = GfaTextareaParser()
    parser.feed(html)

    gfa = unescape("".join(parser.parts)).strip()
    if not gfa:
        raise SystemExit(
            "No GFA text found. Expected a <textarea id=\"gfa-text\"> element."
        )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(gfa + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
