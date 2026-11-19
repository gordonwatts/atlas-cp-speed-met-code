"""Measure original-file basket requests during an otherwise ordinary C endpoint."""
from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
import json
import runpy
import sys
from unittest.mock import patch

# Set the wheel's C++ runtime before importing uproot's native dependencies.
import atlas_met_portable
import typer
import uproot

app = typer.Typer(add_completion=False)
SOURCE = uproot.source.file.MemmapSource
BRANCH = uproot.behaviors.TBranch.TBranch


def union_size(intervals) -> int:
    total = right = 0
    for first, last in sorted(intervals):
        if last > right:
            total += last - max(right, first)
            right = last
    return total


class BasketAudit:
    def __init__(self, input_file: Path):
        self.input_file = input_file
        self.requests = Counter()
        self.planned = defaultdict(set)
        self.branch_roles = defaultdict(set)
        self.leaves = {}
        self.embedded = defaultdict(set)
        self.events = 0

    def original(self, owner) -> bool:
        return Path(owner.file_path) == self.input_file

    def source_read(self, intervals) -> None:
        self.requests.update((int(first), int(last)) for first, last in intervals)

    def record_branch(self, branch, role) -> str:
        path = branch.object_path
        self.branch_roles[path].add(role)
        if path not in self.leaves:
            leaves = []
            for leaf in branch.member("fLeaves"):
                count = leaf.member("fLeafCount", none_if_missing=True)
                leaves.append({"name": leaf.member("fName"), "class": leaf.classname,
                               "count_leaf": None if count is None else count.member("fName")})
            count_branch = branch.count_branch
            self.leaves[path] = {
                "name": branch.name,
                "fLeaves": leaves,
                "typename": branch.typename,
                "declared_count_branch": None if count_branch is None else count_branch.object_path,
            }
        return path

    def planned_basket(self, branch, index, ranges_or_basket, role) -> None:
        path = self.record_branch(branch, role)
        if isinstance(ranges_or_basket, tuple):
            self.planned[tuple(map(int, ranges_or_basket))].add((path, int(index)))
        else:
            self.embedded[path].add(int(index))

    def capture(self, demo: Path, histogram: Path, nominal: bool) -> None:
        chunk = SOURCE.chunk
        chunks = SOURCE.chunks
        plan = BRANCH.entries_to_ranges_or_baskets
        basket = BRANCH.basket

        def one(source, start, stop):
            result = chunk(source, start, stop)
            if self.original(source):
                self.source_read(((start, stop),))
            return result

        def many(source, ranges, notifications):
            result = chunks(source, ranges, notifications)
            if self.original(source):
                self.source_read(ranges)
            return result

        def partition(branch, first, last):
            result = plan(branch, first, last)
            if self.original(branch.file):
                self.events = max(self.events, int(last))
                for index, entry in result:
                    self.planned_basket(branch, index, entry, "uproot_projection")
            return result

        def direct(branch, index):
            result = basket(branch, index)
            if self.original(branch.file):
                if index < branch._num_normal_baskets:
                    first = branch.member("fBasketSeek")[index]
                    entry = (first, first + branch.basket_compressed_bytes(index))
                else:
                    entry = result
                self.planned_basket(branch, index, entry, "direct_basket")
            return result

        args = [str(demo), str(self.input_file), str(histogram), "--stop", "20000"]
        if nominal:
            args.append("--nominal-only")
        previous = sys.argv
        try:
            sys.argv = args
            with (patch.object(SOURCE, "chunk", one), patch.object(SOURCE, "chunks", many),
                  patch.object(BRANCH, "entries_to_ranges_or_baskets", partition),
                  patch.object(BRANCH, "basket", direct)):
                try:
                    runpy.run_path(str(demo), run_name="__main__")
                except SystemExit as exit:
                    if exit.code:
                        raise
        finally:
            sys.argv = previous
        if self.events != 20000 or not self.requests or not histogram.is_file():
            raise RuntimeError("No complete original-file source audit")

    def summarize(self, output: Path, nominal: bool) -> None:
        data = self.requests.keys() & self.planned.keys()
        non_baskets = self.requests.keys() - data
        if not data or self.embedded:
            raise RuntimeError(f"Expected normal Athena baskets, embedded={dict(self.embedded)}")
        size = sum(last - first for first, last in data)
        if union_size(data) != size or any(len(self.planned[extent]) != 1 for extent in data):
            raise RuntimeError("Shared or overlapping branch basket extents")
        by_branch = defaultdict(Counter)
        header_bytes = 0
        payload_bytes = 0
        with uproot.open(self.input_file, handler=SOURCE) as original:
            for first, last in sorted(data):
                parent = original["CollectionTree"]
                header = original.file.source.chunk(first, first + uproot.reading._key_format_big.size)
                key = uproot.reading.ReadOnlyKey(
                    header, uproot.source.cursor.Cursor(first), {}, original.file, parent, read_strings=False
                )
                if key.fSeekKey != first or key.fNbytes != last - first:
                    raise RuntimeError(f"Invalid basket byte extent {(first, last)}")
                (path, index), = self.planned[first, last]
                counts = by_branch[path]
                counts["distinct_baskets"] += 1
                counts["requests"] += self.requests[first, last]
                counts["compressed_bytes_including_key"] += last - first
                counts["uncompressed_bytes_including_key"] += key.data_uncompressed_bytes + key.fKeylen
                counts["compressed_payload_bytes"] += last - first - key.fKeylen
                counts["uncompressed_payload_bytes"] += key.data_uncompressed_bytes
                header_bytes += key.fKeylen
                payload_bytes += key.data_uncompressed_bytes
        combined_requests = sum((last - first) * repeats for (first, last), repeats in self.requests.items())
        data_requests = sum((last - first) * self.requests[first, last] for first, last in data)
        non_basket_extent = union_size(non_baskets)
        non_basket_overlay = size + non_basket_extent - union_size(self.requests)
        branch_rows = {
            path: {**self.leaves[path], **counts,
                   "access": sorted(self.branch_roles[path]),
                   "declared_count_branch_also_fetched": self.leaves[path]["declared_count_branch"] in by_branch}
            for path, counts in sorted(by_branch.items())
        }
        receipt = {
            "scope": "nominal" if nominal else "full",
            "events": self.events,
            "original_file": str(self.input_file),
            "unit": "bytes (KB is decimal 1000 bytes)",
            "interpretation": {
                "instrumented": "Actual MemmapSource chunk/chunks intervals intersected with planned leaf/direct TBasket intervals.",
                "distinct_baskets": "Disjoint on-disk extents; compressed and decompressed TKey totals counted once.",
                "excluded": "Metadata and other nonbasket source requests; logical request duplicates; no kernel readahead/page-cache or physical IO claim.",
                "fLeaves": "ROOT branch fLeaves names and TLeaf::fLeafCount names from original metadata. Fetched count-branch test uses the same recorded branch census.",
                "uncompressed": "Uncompressed payload bytes from TKey fObjlen; key-inclusive adds fKeylen. Includes actual on-disk basket record overhead in the key-inclusive compressed total.",
            },
            "unique_branch_basket_extents": len(data),
            "unique_actual_branches": len(branch_rows),
            "unique_basket_compressed_bytes_including_key": size,
            "unique_basket_uncompressed_bytes_including_key": payload_bytes + header_bytes,
            "unique_basket_compressed_payload_bytes": size - header_bytes,
            "unique_basket_uncompressed_payload_bytes": payload_bytes,
            "basket_compressed_KB_per_event_including_key": size / self.events / 1000,
            "basket_uncompressed_KB_per_event_including_key": (payload_bytes + header_bytes) / self.events / 1000,
            "duplicate_basket_requests": sum(self.requests[row] - 1 for row in data),
            "logical_duplicate_basket_requested_bytes": data_requests - size,
            "planned_ranges_not_requested": len(self.planned.keys() - data),
            "non_basket_range_count": len(non_baskets),
            "non_basket_bytes_outside_baskets": non_basket_extent - non_basket_overlay,
            "non_basket_range_bytes_overlapping_baskets": non_basket_overlay,
            "all_logical_source_request_bytes": combined_requests,
            "all_unique_source_range_bytes": union_size(self.requests),
            "branches": branch_rows,
        }
        with output.open("x") as handle:
            json.dump(receipt, handle, indent=2)
            handle.write("\n")


@app.command()
def diagnose(
    input_file: Path, demo: Path, output: Path, nominal: bool = False,
) -> None:
    output.mkdir()
    audit = BasketAudit(input_file.absolute())
    audit.capture(demo.absolute(), output / "histogram.npz", nominal)
    audit.summarize(output / "logical-root-requests.json", nominal)


if __name__ == "__main__":
    app()
