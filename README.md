# IDR Explorer

The disordered regions of the reviewed human proteome, read two ways, with the modifications that tune their charge:

- **IDRs**: every MobiDB disordered region of 20 or more residues (33,542 on 12,297 proteins), with its charge, its
  charge patterning, its acidic stretches, RG/RGG clusters and basic patches, and its sites.
- **dtIDRs** (domain-tethered IDRs): the window beside every folded-domain edge (121,444), scored for a **basic pole**
  that arginine and lysine writers act on and an **acidic pole** that CK2 and the glutamylases act on.

The site is in [`docs/`](docs/) (GitHub Pages: `main`, folder `/docs`): an Ideas page that states each claim with its
numbers and figure, the two tables, and a map of every protein. Methylarginine evidence at residue level is in the
[Methylarginine Site Explorer](https://shechterlab.github.io/PTM_Pipeline/).

## idrkit

The sequence definitions the lab's IDR work shares, so each is written once:

| | |
| --- | --- |
| `charge_vector`, `ncpr`, `fcr`, `HIS` | residue charge (histidine under the value each tool uses) |
| `delta`, `delta_z` | Das & Pappu's delta (blobs of 5 and 6, as localCIDER computes it), and its z-score against shuffles of the same residues: local charge blocks |
| `scd`, `scd_z` | sequence charge decoration, and its shuffle z-score: long-range segregation |
| `acidic_stretches`, `acidic_runs`, `basic_patches`, `rg_clusters` | the sequence elements of the two poles |
| `e_tract_count`, `ck2_motif_count`, `rg_motif_count`, `longest_acidic_run` | the grammars the dtIDR scores count |
| `idr_intervals`, `merge_intervals` | IDRs from disorder intervals |

```bash
pip install "idrkit @ git+https://github.com/Shechterlab/idr-explorer"
python -c "import idrkit as K; print(K.delta_z('KKKKKEEEEE' * 6))"
pytest            # in a clone: the tests include localCIDER's delta and a literature SCD value
```

glu_proteomics' `glu.sequence` imports its acidic stretches, RG clusters and charge from here.

## How the data are made

`docs/data/` and `docs/figures/` are generated, never edited, by `Analysis/scripts/build_idr_explorer.py` in the lab's
MicroDomains repository, with [PTM_Pipeline](https://github.com/Shechterlab/PTM_Pipeline), glu_proteomics and this
repository cloned beside it. `docs/data/build.json` records the commit of each repository, the MD5 of every input, the
rules, and every number the Ideas page states. MicroDomains and glu_proteomics are private until publication; the inputs
they read are public and listed below.

| Data | Source | Licence |
| --- | --- | --- |
| reference proteome | UniProt 2026_03 reviewed human | CC BY 4.0 |
| disorder | MobiDB (Piovesan et al. 2025) | CC BY 4.0 |
| domains | InterPro 108.0 | CC0 |
| methylarginines | PTM_Pipeline release v1.4.1 (DeAngelo, Silverstein & Shechter, Cell Chem Biol 2026) | CC BY-NC-SA 4.0 |
| phosphorylation, lysine acylation, ubiquitin/SUMO, lysine methylation | iPTMnet export (Huang et al. 2018) staged in PTM_Pipeline | CC BY-NC-SA 4.0 |
| modified residues with enzymes (CK2), ADP-ribosylation, citrullination, polyglutamylation | UniProt 2026_03 MOD_RES (PTM_Pipeline scripts/53) | CC BY 4.0 |
| acetyl, ubiquitin, SUMO and methyl sites, Gold/Silver/Bronze | PTMeXchange human builds, PRIDE PXD079218, PXD068989, PXD068489, PXD076372 (PTM_Pipeline scripts/54) | CC0 |
| glutamylation | TTLL4 re-analysis of PRIDE PXD065251 (Dey et al., J Cell Sci 2026) | lab analysis |

## Licence

Code (idrkit, the site): MIT. Data in `docs/data/`: CC BY-NC-SA 4.0, because the methylarginine release and iPTMnet
are; the sources keep their own terms (table above). See [LICENSE.md](LICENSE.md).
