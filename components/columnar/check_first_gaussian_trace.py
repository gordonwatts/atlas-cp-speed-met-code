#!/usr/bin/env python3
"""Verify actual first-Gaussian evaluations separately from stock logical draws."""
from collections import Counter
from itertools import groupby, zip_longest
import json
from pathlib import Path
import struct

import typer

ACTIONS = {"rootSeed", "setRandomSeed", "gaussian", "firstGaussian"}


def events(path, *, stop=None, nominal=False):
    def load_rows(handle):
        for line in handle:
            entry,state,action,*fields = line.rstrip("\n").split("\t")
            entry,state = int(entry),int(state)
            if stop is not None and entry >= stop:
                break
            if action in ACTIONS and (not nominal or state == 0):
                yield entry,state,action,fields

    with path.open() as handle:
        for entry,rows in groupby(load_rows(handle),key=lambda row: row[0]):
            yield entry,[(state,action,fields) for _,state,action,fields in rows]


def main(stock_status: Path, columnar_status: Path, output: Path,
         entries: int = typer.Option(20000,min=1), nominal_only: bool = False):
    count = Counter()
    first_failures = []

    def check(entry, field, expected, actual):
        if expected != actual:
            count["failures"] += 1
            if len(first_failures) < 8:
                first_failures.append({"entry":entry,"field":field,
                                       "stock_expected":expected,"candidate":actual})

    stock = events(stock_status,stop=entries,nominal=nominal_only)
    candidate = events(columnar_status)
    for entry,a,c in zip_longest(range(entries),stock,candidate):
        if a is None or c is None or entry is None or a[0] != entry or c[0] != entry:
            check(entry,"entry",entry,(None if a is None else a[0], None if c is None else c[0]))
            break

        setters = [(state,fields[0]) for state,action,fields in a[1] if action=="setRandomSeed"]
        found_setters = [(state,fields[0]) for state,action,fields in c[1] if action=="setRandomSeed"]
        check(entry,"setter",setters,found_setters)
        roots = {state:int(fields[0]) % 2**32 for state,action,fields in a[1] if action=="rootSeed"}
        stock_draws = [(state,fields) for state,action,fields in a[1] if action=="gaussian"]
        found_draws = [(state,fields) for state,action,fields in c[1] if action=="gaussian"]

        def bits(rows):
            return [(state,[struct.pack("=d",float.fromhex(value)).hex() for value in values])
                    for state,values in rows]

        check(entry,"logical_gaussian_bits",bits(stock_draws),bits(found_draws))
        active = {state for state,_ in stock_draws}
        expected = []
        previous = None
        for state,_ in setters:
            key = roots[state]
            use = "skipped"
            if state in active:
                use = "reused" if previous == key else "evaluated"
                previous = key
            expected.append((state,key,use))
        observed = [(state,int(fields[0]),fields[1]) for state,action,fields in c[1] if action=="firstGaussian"]
        check(entry,"first_gaussian_evaluation",expected,observed)
        count["events"] += 1
        count["logical_draws"] += len(found_draws)
        count["physical_evaluations"] += sum(use=="evaluated" for _,_,use in observed)
        count["cache_hits"] += sum(use=="reused" for _,_,use in observed)
        count["skipped_states"] += sum(use=="skipped" for _,_,use in observed)

    report = {"passed":count["events"]==entries and not count["failures"],
              "nominal_only":nominal_only,"counts":dict(count),"first_failures":first_failures}
    with output.open("x") as handle:
        json.dump(report,handle,indent=2)
        handle.write("\n")
    if not report["passed"]:
        raise typer.Exit(1)


if __name__ == "__main__":
    typer.run(main)
