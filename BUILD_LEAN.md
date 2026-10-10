# Building the accompanying Lean proof

The proof uses Lean, Mathlib and Zeta23. Reproduce it in the pinned Zeta23 project with the supplied auxiliary definition file. The release-root `lean-toolchain` records the version; the build commands run inside the Zeta23 checkout, whose Lake configuration supplies the imported libraries.

| Dependency | Exact revision |
|---|---|
| Lean | `leanprover/lean4:v4.33.0-rc2`  |
| [Zeta23](https://github.com/anthropics/zeta-23-lean) | 
| [Mathlib](https://github.com/leanprover-community/mathlib4) | pinned by that Zeta23 checkout |

Keep the checkout's `lake-manifest.json` and `lean-toolchain`. Its `lakefile.toml` declares the `ChallengeDeps` and `Solution` libraries with `srcDir = "comparator"`, so the exact build target is `Solution.Candidate`.

From this release directory, the following commands use a POSIX shell, Git and an installed Lean toolchain manager (`elan`):

```sh
git clone https://github.com/anthropics/zeta-23-lean zeta23-build
git -C zeta23-build checkout --detach 3635e74826a4c1fcece7d1cd2b6fa75e43a00510
mkdir -p zeta23-build/comparator/ChallengeDeps zeta23-build/comparator/Solution
cp lean_support/CandidateSpec.lean zeta23-build/comparator/ChallengeDeps/CandidateSpec.lean
cp r_arun_chandru_proof.lean zeta23-build/comparator/Solution/Candidate.lean
cd zeta23-build
lake exe cache get
lake build Solution.Candidate
```

The cache command obtains the pinned dependencies and available compiled library artifacts. Lake compiles remaining imports as required. Retain sufficient memory and time for the large proof file. The supplied `lean_support/CandidateSpec.lean` is the exact existing 335-byte auxiliary file used with this proof. It contains the fixed comparison and candidate definitions imported by the proof. Preserve both supplied files byte-for-byte, including the proof's attribution and license notices:

The exact arithmetic checks in `reproducibility/` have their own Python/Node runner and can be run independently of the Lean build. License and inherited attribution texts are provided in `LICENSE`, `THIRD_PARTY_NOTICES.txt`, and the frozen proof source.
