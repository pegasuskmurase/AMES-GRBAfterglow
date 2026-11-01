# Developer notes and initial audit

## Scope and method

This is a static source audit of `include/GRBAfterglow.hpp`,
`src/GRBAfterglow.cpp`, and `src/module/main_bindings.cpp`.
The constructor and distinct public method names are inventoried below;
overloads share one row. “Calls” counts same-file method-name call expressions,
excluding definitions and comments. It is a review aid, not a whole-program
reachability analysis. No callers found does **not** mean safe to remove:
C++ clients outside this repository may call public methods.

The constructor `GRBafterglow(Source&)` is implemented and bound as
`AMES.GRBafterglow(source)`. The engine references its Source rather than owning it.

## Public method inventory

| C++ method | Definition in .cpp | Python name | Same-file calls |
| --- | --- | --- | --- |
| `help` | Yes | Not bound | 0 |
| `setGRBafterglowParam` | Yes | `setGRBafterglowParam` | 0 |
| `setGRBafterglowRSParam` | Yes | `setGRBafterglowRSParam` | 0 |
| `setExtendedEmissionParam` | Yes | Not bound | 0 |
| `setPlateauEmissionParam` | Yes | Not bound | 0 |
| `setOutputFolder` | Yes | `setOutputFolder` | 0 |
| `showInfo` | Yes | `showInfo` | 0 |
| `haveReverseShock` | Yes | `haveReverseShock` | 0 |
| `haveOneZone` | Yes | `haveOneZone` | 0 |
| `haveSSCSpec` | Yes | `haveSSCSpec` | 0 |
| `haveEICSpec` | Yes | Not bound | 0 |
| `haveEdgeEffect` | Yes | Not bound | 0 |
| `haveAttenuSSA` | Yes | `haveAttenuSSA` | 0 |
| `haveAttenuFF` | Yes | Not bound | 0 |
| `haveAttenuGGSource` | Yes | `haveAttenuGGSource` | 0 |
| `haveAttenuGGCosmic` | Yes | Not bound | 0 |
| `FluxOutput` | Yes | `FluxOutput` | 0 |
| `Flux1D` | Yes | `Flux1D` | 0 |
| `Flux` | Yes | `Flux` | 0 |
| `setParamsFromMCMC` | Yes | `set_params_mcmc` | 0 |
| `CalculateLightCurveMCMC` | Yes | `calculate_lightcurve` | 0 |
| `Spectrum` | Yes | `Spectrum` | 5 |
| `SpectrumOneZone` | Yes | Not bound | 1 |
| `ReverseShock` | Yes | Not bound | 1 |
| `InitJet` | Yes | Not bound | 5 |
| `EvolutionThickShell` | Yes | Not bound | 2 |
| `calcDynamic` | Yes | Not bound | 0 |
| `PromptPhotonDensity` | Yes | Not bound | 0 |
| `ExtendedEmissionLuminosity` | Yes | Not bound | 1 |
| `PlateauEmissionLuminosity` | Yes | Not bound | 1 |
| `ExternalDensity` | Yes | Not bound | 3 |
| `MomentumInj` | Yes | Not bound | 4 |
| `MomentumMax` | Yes | Not bound | 5 |
| `DecelerationRadius` | Yes | Not bound | 3 |
| `EvolutionThinShellAnalytic` | Yes | Not bound | 1 |
| `LorentzBoost` | Yes | Not bound | 2 |
| `FluxNorm` | Yes | Not bound | 2 |
| `MomentumInjProton` | Not found | Not bound | 0 |
| `MomentumMaxDyn` | Yes | Not bound | 0 |

## Findings to review

- `MomentumInjProton` is declared but has no definition in this .cpp and no
  repository references beyond its declaration were found. It is not bound.
  Confirm intended scope before implementing or removing it.
- `haveAttenuFF` and `haveAttenuGGCosmic` store flags, but the corresponding
  flags have no active reads in the audited .cpp. These switches should not
  be advertised as functioning attenuation features without implementation review.
- `help` and `calcDynamic` have definitions but are not bound to Python.
- `EvolutionThinShellAnalytic`, `InitJet`, `EvolutionThickShell`, radiation
  helpers, and density calculations are implementation machinery even though
  several are declared public. A later refactor may narrow visibility.
- The EIC path exists, but the prompt-density call inside the multi-zone path
  is commented out. Audit external target setup before documenting EIC as ready
  for general use. Its setters and switch are not currently bound to Python.
- Parameter counts and return types in some binding docstrings disagree with
  the implementation. The user guide records actual behavior; fix docstrings
  in a separate code change.
- Array-length validation is incomplete in numerical entry points. MCMC input
  is indexed directly. Add input checks before recommending these to new users.
- Source lifetime is managed by the example object, but the binding currently
  does not declare a keep-alive policy. Retain the Source explicitly.

No C++ functions were deleted or behavior changed for this documentation pass.
This inventory is not a scientific validation of every branch.

## Maintaining and previewing the documentation

From the repository root, preferably in a separate documentation environment:

```sh
python -m pip install -r requirements-docs.txt
python -m mkdocs serve
```

Open the local URL printed by MkDocs. Edit Markdown files under `docs/` and
refresh the preview. Check the static output with:

```sh
python -m mkdocs build --strict
```

This builds HTML under `site/` without compiling or importing AMES.
The MkDocs configuration uses the bundled Read the Docs theme; no custom
website code is required. See the [MkDocs guide](https://www.mkdocs.org/user-guide/).

## Optional publication later

Review the generated pages first. For GitHub Pages, MkDocs provides:

```sh
python -m mkdocs gh-deploy
```

This command commits generated pages and pushes a `gh-pages` branch to your
configured remote. Run it only when ready to publish, then configure the
repository's Pages settings to serve that branch. The current change does not
publish anything or enable automatic deployments. Read the Docs hosting can
also be configured later; it is not required for the local preview.
See [MkDocs deployment documentation](https://www.mkdocs.org/user-guide/deploying-your-docs/).
