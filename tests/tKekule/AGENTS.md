# Kekule Bond-Order Optimizer Tests

## Purpose

Tests for the C++ Kekule bond-order optimizer (`pyBall.Kekule`): bond-order relaxation on aromatic fragments and generation of Kekule structure variants with donor/acceptor groups.

## Ownership

- Bond-order relaxation test script
- Kekule variant generation with substituents (-OH, -NH2, =O, =N-)

## Local Contracts

- **Run from this directory** — scripts use relative paths to data.
- **Use `run.sh` / `make.sh` scripts** when available; never invoke `make` directly.

## Work Guidance

- `run.py` — test bond-order relaxation on anthracene fragments
- `generate_2.py` — generate Kekule structure variants with donor/acceptor groups

## Verification

- Run `python3 run.py` from this directory
- Bond orders should converge to chemically meaningful values (alternating single/double)

## Child DOX Index

- No child subtrees
