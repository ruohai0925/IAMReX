<!--
SPDX-FileCopyrightText: 2026 Yadong Zeng<zdsjtu@gmail.com>

SPDX-License-Identifier: BSD-3-Clause
-->

# Security Policy

IAMReX is a research CFD code. Most issues in it are correctness or
performance bugs rather than security vulnerabilities, but a few classes of
report are still security-relevant: memory-safety bugs reachable from
untrusted input files (inputs decks, restart/checkpoint files, mesh files),
or supply-chain issues in the build (e.g. a compromised dependency fetch in
CI).

## Supported Versions

Security fixes are made against the `development` branch and included in
the next tagged release. We do not backport fixes to older release tags.

## Reporting a Vulnerability

Please do **not** open a public GitHub issue for a suspected security
vulnerability. Instead, report it privately to ruohai372@gmail.com with:

- A description of the issue and its potential impact.
- Steps to reproduce it (a minimal inputs file or test case helps a lot).
- The IAMReX commit or release tag affected.

We will acknowledge your report, investigate, and let you know the planned
fix and disclosure timeline. We're a small academic team, so please allow a
few days for a first response.
